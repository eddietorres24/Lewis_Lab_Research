library(dplyr)
library(rtracklayer)
library(ggplot2)
library(tidyr)

#Read in PRC2 targets
Prc2targets <- read.table("./bed_files/K27_genes_trimmed.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

###Calculate ChIP-seq Signal Averages

# Load your gene/promoter BED file as GRanges
colnames(K27_genes)[1:3] <- c("chr", "start", "end")  # set column names if missing
gene_ranges <- GRanges(
  seqnames = K27_genes$chr,
  ranges = IRanges(
    start = K27_genes$start,
    end = K27_genes$end
  ),
  strand = K27_genes$V6
)
mcols(gene_ranges)$gene_id <- K27_genes$V10


```{r, Read in Data and Run Diffbind}
#Read in H3K27me3 ChIP Data (This should contain 2+ replicates of ChIP data & controls with your desired modification & strains. we don't need peak calls here, since we are assigning our own regions)
## reading in my sample sheet w/ all CAF-1 data and just using WT strains
chip_data <- read.csv("../DiffBind_CAF-1_K27_CS_nopeak_ET.csv")
chip_data <- chip_data[1:12,]

# Fix paths and get rid of peaks (dont need to fix paths if in correct directory)
chip_data$bamControl <- paste0("../", chip_data$bamControl)
chip_data$bamReads <- paste0("../", chip_data$bamReads)

#Make a diffbind object
chip_data <- dba(sampleSheet = chip_data)

###NOTE: unfortunately, we will have to run diffibnd on the + & - strand genes separately. This is because if 2 genes overlap, diffbind will combine them into one "peak". ideally, we want to measure the signal over each gene individually to get an accurate picture of K27 coverage. This is the strategy ChatGPT helped me devise to address this. There may be a better way to do this, but this also wasn't too difficult.

plus_r  <- genes_gr[strand(genes_gr) == "+"]
minus_r <- genes_gr[strand(genes_gr) == "-"]

# Run dba.count separately
chip_data_plus  <- dba.count(chip_data, peaks = plus_r, summits = FALSE, filter = 0, minOverlap = 0)
chip_data_minus <- dba.count(chip_data, peaks = minus_r, summits = FALSE, filter = 0, minOverlap = 0)

#dba analyze to obtain differential binding
chip_data_plus <- dba.analyze(chip_data_plus)
chip_data_minus <- dba.analyze(chip_data_minus)

# Extract differentially bound sites for each contrast (automatically defined by dba.analyze)
plus_db  <- dba.report(chip_data_plus, th=1)  # th=1 returns all sites regardless of FDR
minus_db <- dba.report(chip_data_minus, th=1)

# Convert to data.frame
plus_df  <- as.data.frame(plus_db)
minus_df <- as.data.frame(minus_db)

# Add strand info for reference
plus_df$strand  <- "+"
minus_df$strand <- "-"

# Combine
combined_df <- rbind(plus_df, minus_df)

# Optional: Create unique row names
rownames(combined_df) <- paste0(combined_df$Chr, ":", combined_df$Start, "-", combined_df$End, "(", combined_df$strand, ")")







# # Extract counts as GRanges
# counts_plus  <- dba.peakset(chip_data_plus,  bRetrieve = TRUE)
# counts_minus <- dba.peakset(chip_data_minus, bRetrieve = TRUE)
# 
# # Combine into full GRanges object
# full_counts <- c(counts_plus, counts_minus)
# 
# # Extract metadata (sample names, counts matrix)
# count_matrix <- as.data.frame(mcols(full_counts))
# 
# # Extract genomic coordinates
# region_info <- data.frame(
#   Chromosome = as.character(seqnames(full_counts)),
#   Start = start(full_counts),
#   End   = end(full_counts)
# )
# 
# # Combine into final count dataframe with coordinate info
# norm_counts <- cbind(region_info, count_matrix)
```

###############################################################
# Step 1: Split gene_ranges by strand
plus_ranges  <- gene_ranges[strand(gene_ranges) == "+"]
minus_ranges <- gene_ranges[strand(gene_ranges) == "-"]

# Step 2: Run dba.count separately
chip_data_plus  <- dba.count(chip_data, peaks = plus_ranges, summits = FALSE, filter = 0, minOverlap = 0)
chip_data_minus <- dba.count(chip_data, peaks = minus_ranges, summits = FALSE, filter = 0, minOverlap = 0)

# Step 3: Extract counts as GRanges
counts_plus  <- dba.peakset(chip_data_plus,  bRetrieve = TRUE)
counts_minus <- dba.peakset(chip_data_minus, bRetrieve = TRUE)

# Step 4: Combine into full GRanges object
full_counts <- c(counts_plus, counts_minus)

# Extract metadata (sample names, counts matrix)
count_matrix <- as.data.frame(mcols(full_counts))

# Extract genomic coordinates
region_info <- data.frame(
  CHR   = as.character(seqnames(full_counts)),
  START = start(full_counts),
  END   = end(full_counts),
  STRAND = as.character(strand(full_counts))
)

# Combine into final count dataframe
norm_counts <- cbind(region_info, count_matrix)

# Ensure coordinate columns match in name
####################################################################


###DEseq2 for ChIP data

# Step 1: Convert full_counts to a count matrix
count_matrix <- as.data.frame(mcols(full_counts))
rownames(count_matrix) <- paste0(seqnames(full_counts), ":", start(full_counts), "-", end(full_counts))

# Step 2: Create sample information (colData)
sample_info <- data.frame(
  sample = colnames(count_matrix),
  strain = gsub("_CS_K27_Rep_\\d+", "", colnames(count_matrix))  # Extract strain names
)
rownames(sample_info) <- sample_info$sample

# Step 3: Make DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),         # DESeq2 expects integer counts
  colData = sample_info,
  design = ~ strain
)

# Step 4: Run DESeq
dds2 <- DESeq(dds)
plotDispEsts(dds2)

# Step 5: Extract results (mutants vs. WT)
res_cac1 <- results(dds2, contrast = c("strain", "cac.1", "WT"))
res_cac2 <- results(dds2, contrast = c("strain", "cac.2", "WT"))
res_cac3 <- results(dds2, contrast = c("strain", "cac.3", "WT"))

#Step 6: Make new df w/ log2FC values

#Extract coordinates from DESeq2 results
coord_split <- do.call(rbind, strsplit(rownames(res_cac1), "[:-]"))
summary_df <- data.frame(
  CHR = coord_split[, 1],
  START = as.numeric(coord_split[, 2]),
  END = as.numeric(coord_split[, 3]),
  log2FC_cac1_vs_WT = res_cac1$log2FoldChange,
  log2FC_cac2_vs_WT = res_cac2$log2FoldChange,
  log2FC_cac3_vs_WT = res_cac3$log2FoldChange
)

#Merge with known gene IDs from norm_counts
summary_df <- merge(
  summary_df,
  norm_counts[, c("CHR", "START", "END", "gene_id")],
  by = c("CHR", "START", "END"),
  all.x = TRUE
)

#Reorder columns
summary_df <- summary_df[, c("gene_id", "CHR", "START", "END",
                             "log2FC_cac1_vs_WT", "log2FC_cac2_vs_WT", "log2FC_cac3_vs_WT")]

###Incorporate Analyzed RNA-seq data

###NEED TO RUN DEseq TO OBTAIN RNA-seq Log2FC VALUES###
##assign RNA-seq Log2FC for each strain after running DEseq (I am not eliminating non-DE genes in this code because i want both DE & non-DE genes present on my plot since we are correlating DE w/ K27 changes)
cac1_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac1_new_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac2_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac2_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac3_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac3_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac1_cac2_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac1_cac2_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)

#append RNA-seq log2FC from DEseq to the end of the ChIP log2FC df
summary_df$cac1RNA <- cac1_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac1_seq))]
summary_df$cac2RNA <- cac2_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac2_seq))]
summary_df$cac3RNA <- cac3_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac3_seq))]
summary_df$cac1_cac2RNA <- cac1_cac2_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac1_cac2_seq))]

#remove any rows w/ NA
summary_df <- na.omit(summary_df)

#statistics
cor_test <- cor.test(summary_df$log2FC_cac3_vs_WT, summary_df$cac3RNA, method = "pearson")
r2 <- round(cor_test$estimate^2, 2)
pval <- signif(cor_test$p.value, 3)

#Scater Plot
ggplot(summary_df, aes(x = log2FC_cac3_vs_WT, y = cac3RNA)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  annotate("text", x = -7.5, y = 12.5, 
           label = paste0("R² = ", r2, ", p = ", pval), hjust = 0) +
  theme_minimal() +
  labs(
    x = "H3K27me3 Signal (cac-3/WT)",
    y = "log2 Fold Change (cac-3/WT)",
    title = "ChIP-seq vs RNA-seq Scatter Plot"
  )


ggsave(filename = "./K27_v_RNA_cac3_Paper.pdf", plot = plot, dpi=600, height=9, width=12, units = "in")

###TO-DO###
#1. color genes that are located within ectopic CAF-1 peaks
#2. 

