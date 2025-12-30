# antisense_deseq_minimal.R
suppressPackageStartupMessages({
  library(DESeq2)
})
#libraries
library("DESeq2")
library("ggplot2")
library("gplots")
library("ggrepel")
library("dplyr")
library("pheatmap")
library("grid")
library("corrplot")
library("RColorBrewer")
library("scales")
library("readr")
library("tidyr")
library("tidyverse")
library("scater")
library("reshape2")
library("grDevices")
library("tibble")

## --- USER PATHS (edit these two) ---
BASE_DIR   <- "C:/Users/eddie/Research/GitHub/Research"
COLDATA_CSV <- "C:/Users/eddie/Research/GitHub/coldata_anti.csv"

## If your copied counts live in BASE_DIR/counts/..., this is fine.
## Otherwise, point COUNTS_DIR to wherever the .antisense.txt files are.
COUNTS_DIR <- file.path(BASE_DIR, "counts")

## Switch between "antisense" and "sense" by changing this:
COUNT_KIND <- "antisense"

#--- helpers ---------------------------------------------------------------

findCountsFile <- function(accession, kind = c("antisense","sense"),
                           counts_dir = COUNTS_DIR) {
  kind <- match.arg(kind)
  # candidate layouts:
  cand <- c(
    file.path(counts_dir, accession, paste0(accession, ".", kind, ".txt")),
    file.path(counts_dir, paste0(accession, ".", kind, ".txt"))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) {
    stop(sprintf("Counts file not found for %s (%s): tried\n  - %s\n  - %s",
                 accession, kind, cand[1], cand[2]), call. = FALSE)
  }
  hit[1]
}

read_fc_one <- function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   comment.char = "#", check.names = FALSE, quote = "")
  gene  <- df[[1]]                 # "Geneid" column
  count <- df[[ncol(df)]]          # the only sample column in each file
  data.frame(gene = gene, cnt = as.numeric(count), 
             stringsAsFactors = FALSE)
}

#--- load coldata ----------------------------------------------------------

coldata <- read.csv(COLDATA_CSV, row.names = 1, stringsAsFactors = FALSE)
stopifnot(all(c("condition","type","file_id") %in% colnames(coldata)))
coldata$condition <- factor(coldata$condition)
coldata$type      <- factor(coldata$type)

# Basenames to look up files by
samples <- rownames(coldata)
file_ids <- setNames(coldata$file_id, samples)

#--- build counts matrix ---------------------------------------------------

message("Locating count files...")
files <- vapply(file_ids, findCountsFile, FUN.VALUE = character(1),
                kind = COUNT_KIND, counts_dir = COUNTS_DIR)

message("Reading featureCounts tables...")
lst <- lapply(files, read_fc_one)

# merge by gene ID
merged <- Reduce(function(x,y) merge(x,y, by = "gene", all = TRUE), lst)
rownames(merged) <- merged$gene
merged$gene <- NULL
colnames(merged) <- samples  # label columns by sample_id

#subset merged by genes that dont overlap (from previously mad beds)
merged = subset(merged, rownames(merged) %in% plusgene$V10 | rownames(merged) %in% minusgene$V10 )

# replace NAs with 0 and coerce to integer matrix
merged[is.na(merged)] <- 0
cts <- as.matrix(round(merged, 0))

# ensure order matches coldata
cts <- cts[, rownames(coldata), drop = FALSE]

# quick write-out for inspection
write.csv(data.frame(Geneid = rownames(cts), cts, check.names = FALSE),
          file = file.path(BASE_DIR, paste0("combined_", COUNT_KIND, "_counts.csv")),
          row.names = FALSE)

#--- minimal DESeq2 --------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = ~ condition)

# light prefilter: keep genes with at least 10 total counts
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# set reference
if ("WT" %in% levels(dds$condition)) {
  dds$condition <- relevel(dds$condition, ref = "WT")
}

dds <- DESeq(dds)

# example results: set7 vs WT (change the contrast as you like)
if ("set7" %in% levels(dds$condition)) {
  res_set7 <- results(dds, contrast = c("condition","set7","WT"))
  write.csv(as.data.frame(res_set7),
            file = file.path(BASE_DIR, "DE_set7_vs_WT_antisense.csv"))
}

# save the DESeq2 object for downstream work
saveRDS(dds, file = file.path(BASE_DIR, paste0("dds_", COUNT_KIND, ".rds")))

##setting p-value cutoff, starting with 0.5
alpha = 0.1

#Subsetting results based on p-value
cac1_old_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac1_old", "WT"))
cac2_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac2", "WT"))
cac3_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac3", "WT"))
asf1_ant <- results(dds, alpha = 0.1, contrast=c("condition", "asf1", "WT"))
set7_ant <- results(dds, alpha = 0.1, contrast=c("condition", "set7", "WT"))
cac1_new_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac1_new", "WT"))
cac1_cac2_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac1_cac2", "WT"))
cac1_suz12_ant <- results(dds, alpha = 0.1, contrast=c("condition", "cac1_suz12", "WT"))
suz12_ant <- results(dds, alpha = 0.1, contrast=c("condition", "suz12", "WT"))
ash1_ant <- results(dds, alpha = 0.1, contrast=c("condition", "ash1", "WT"))
dpf3_ant <- results(dds, alpha = 0.1, contrast=c("condition", "peri3dpf", "WT"))
dpf4_ant <- results(dds, alpha = 0.1, contrast=c("condition", "peri4dpf", "WT"))
dpf5_ant <- results(dds, alpha = 0.1, contrast=c("condition", "peri5dpf", "WT"))
dpf6_ant <- results(dds, alpha = 0.1, contrast=c("condition", "peri6dpf", "WT"))

#Make csvs
write.csv(as.data.frame(cac1_old_ant), file = "cac1_ant.csv")
write.csv(as.data.frame(cac2_ant), file = "cac2_ant.csv")
write.csv(as.data.frame(cac3_ant), file = "cac3_ant.csv")
write.csv(as.data.frame(asf1_ant), file = "asf1_ant.csv")
write.csv(as.data.frame(set7_ant), file = "set7_ant.csv")
write.csv(as.data.frame(cac1_new_ant), file = "cac1_new_ant.csv")
write.csv(as.data.frame(cac1_cac2_ant), file = "cac1_cac2_ant.csv")
write.csv(as.data.frame(cac1_suz12_ant), file = "cac1_suz12_ant.csv")
write.csv(as.data.frame(suz12_ant), file = "suz12_ant.csv")
write.csv(as.data.frame(ash1_ant), file = "ash1_ant.csv")
write.csv(as.data.frame(dpf3_ant), file = "3dpf_ant.csv")
write.csv(as.data.frame(dpf4_ant), file = "4dpf_ant.csv")
write.csv(as.data.frame(dpf5_ant), file = "5dpf_ant.csv")
write.csv(as.data.frame(dpf6_ant), file = "6dpf_ant.csv")


##############Making IGV Track
#calling list of files
##Note: this will only work if the ONLY csv's in this directory are the ones containing the Data you want on the track (i.e. you need to move any other csv's from this directory, or move the ones you want on the track to a new one and change the path)
list_of_files <- list.files(path = "../Research",
                            recursive = FALSE,
                            pattern = ".csv$")

setwd("C:/Users/eddie/Research/GitHub/Research")

induced_genes <- readr::read_csv(list_of_files, id = "file_name")

##convert to wide format
l2fc_data <- data.frame(pivot_wider(data = induced_genes, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))

#assigning row names and deleting extra column
rownames(l2fc_data) <- l2fc_data[,1]
l2fc_data <- l2fc_data[,-1]

#assigning NA values to 0
l2fc_data[is.na(l2fc_data)] <- 0

#tibble NCU names back to first column
l2fc <- tibble::rownames_to_column(l2fc_data, "VALUE")

#read in file w/ all genes
ngenes = read.csv("../csv_files/neurospora_genes_edit.csv")

#Make IGV Track
##Note: only genes with DE in one or more strains will display data in the track
combine <- transform(merge(ngenes, l2fc, by = "VALUE"))
track = combine %>% select("SequenceID", "FeatureStart", "FeatureEnd", "VALUE", "asf1_ant.csv", "cac1_ant.csv", "cac2_ant.csv", "cac3_ant.csv", "set7_ant.csv", "cac1_new_ant.csv", "cac1_cac2_ant.csv", "cac1_suz12_ant.csv", "suz12_ant.csv", "ash1_ant.csv", "X3dpf_ant.csv", "X4dpf_ant.csv", "X5dpf_ant.csv", "X6dpf_ant.csv")

#write track to file
write.table(track, file="antisense_dev.igv", sep="\t", row.names = FALSE, quote = FALSE)



#####Making heatmap

##Run calculateTPM function - must load "scater" library
Interactors_tpm <- calculateTPM(merged, lengths = countdataInteractors[,5])

##convert to matrix again
allDataTPM  <- data.matrix(Interactors_tpm)

#rename columns with human readable sample names
##Note: THIS NEXT SECTION REQUIRED MANUAL INPUT INTO A SPREADSHEET OR TEXT EDITOR, this code only reads in the file that was already renamed

names <- read.table(file="names.txt" , sep="", row.names = 1, header = TRUE)
colnames(allDataTPM) <- names[,1]

#Make new matrix
Ordered_KO_data <- cbind(allDataTPM[,13:16],allDataTPM[,31:33],allDataTPM[,20:22],allDataTPM[,4:6],allDataTPM[,23:25],allDataTPM[,26:28],allDataTPM[,29:30],allDataTPM[,34:35],allDataTPM[,36:37],allDataTPM[,38:39],allDataTPM[,40:41])
Averaged_Orderd_KO_data <- cbind(rowMeans(allDataTPM[,13:16], na.rm = TRUE),
                                 rowMeans(allDataTPM[,31:33], na.rm = TRUE),
                                 rowMeans(allDataTPM[,20:22], na.rm = TRUE),
                                 rowMeans(allDataTPM[,4:6], na.rm = TRUE),
                                 rowMeans(allDataTPM[,23:25], na.rm = TRUE),
                                 rowMeans(allDataTPM[,26:28], na.rm = TRUE),
                                 rowMeans(allDataTPM[,29:30], na.rm = TRUE),
                                 rowMeans(allDataTPM[,34:35], na.rm = TRUE),
                                 rowMeans(allDataTPM[,36:37], na.rm = TRUE),
                                 rowMeans(allDataTPM[,38:39], na.rm = TRUE),
                                 rowMeans(allDataTPM[,40:41], na.rm = TRUE))
#Rename columns
averageRowIDs=c("WT","suz12","cac-1","cac-2","cac-1_2","cac-1_suz12","ash-1","3dpf","4dpf","5dpf","6dpf")
colnames(Averaged_Orderd_KO_data) <- averageRowIDs

#subset data to sample only PRC2-target genes
##reading in geneIDs of PRC2 target genes (promoters and genes fully covered by K27)
Prc2targets <- read.table("../bed_files/K27_genes_stringent.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

Prc2targetTPM <- subset(allDataTPM, rownames(allDataTPM)%in%Prc2targets[,10])

#This gene gets left out becausei had to rename it in my bed file due to blanks. I need to add it back manually because the subsetting doesn't recognize the new name
lac <- allDataTPM[rownames(allDataTPM) == "laccase precursor", ]
Prc2targetTPM <- rbind(Prc2targetTPM, lac)

#rename rows of added gene back to original name
rownames(Prc2targetTPM)[rownames(Prc2targetTPM) == "lac"] <- "laccase precursor"

#Subset Data (n = 531)
AVERAGE_Prc2targetTPM <- subset(Averaged_Orderd_KO_data, rownames(Averaged_Orderd_KO_data)%in%rownames(Prc2targetTPM))

#Subset data to filter out non-repressed PRC2 targets regions that got through (cutting out any genes over 10 tpm in WT). Choice of 10 TPM cutoff was arbitrary
##These are likely genes whose promoters are not marked by H3K27me3, genes on the edge of K27 regions, and/or bivalent genes
AVERAGE_Prc2targetTPM <- subset(AVERAGE_Prc2targetTPM, (AVERAGE_Prc2targetTPM[,1] < 10))
AVERAGE_AlldataTPM <- subset(Averaged_Orderd_KO_data, (Averaged_Orderd_KO_data[,1] > -0.1))

#resubet PRC2 targets after filtering
Prc2targetTPM <- subset(Prc2targetTPM, rownames(Prc2targetTPM)%in%rownames(AVERAGE_Prc2targetTPM))

#Add sudocount and log transform (if necessary)
AVERAGE_Prc2targetTPM <- AVERAGE_Prc2targetTPM + 1
AVERAGE_Prc2targetTPM <- log2(AVERAGE_Prc2targetTPM)

AVERAGE_AlldataTPM <- AVERAGE_AlldataTPM + 1
AVERAGE_AlldataTPM <- log2(AVERAGE_AlldataTPM)
#melt data to get it into a format ggplot can use (library "reshape2")
meltedAveragePRC2targetData <- reshape2::melt(AVERAGE_Prc2targetTPM, value.name = 'Count',
                                    varnames=c('GeneID', 'Sample'))

meltedAllData <- reshape2::melt(AVERAGE_AlldataTPM, value.name = 'Count',
                      varnames=c('GeneID', 'Sample'))

###Plots
altorder = rev(c("WT","suz12","cac-1","cac-2","cac-1_2","cac-1_suz12","ash-1","3dpf","4dpf","5dpf","6dpf"))

#make factor
meltedAveragePRC2targetData$Sample <- factor(meltedAveragePRC2targetData$Sample)

meltedAllData$Sample <- factor(meltedAllData$Sample)

#histograms
hist <- as.data.frame(AVERAGE_AlldataTPM)

long_hist <- hist %>%
  select(WT, `cac-1_suz12`) %>%
  pivot_longer(cols = everything(), names_to = "Strain", values_to = "TPM")

# Plot overlapping histograms
ggplot(long_hist, aes(x = TPM, fill = Strain)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, color = "black") +
  theme_minimal() +
  labs(title = "TPM Distribution in WT and cac-1",
       x = "TPM",
       y = "Frequency") +
  scale_fill_manual(values = c("WT" = "steelblue", "cac-1_suz12" = "tomato"))


#violin
# Define factor order and formatted x-axis labels
label_order <- c("WT","suz12","cac-1","cac-2","cac-1_2","cac-1_suz12","ash-1","3dpf","4dpf","5dpf","6dpf")
formatted_labels <- c(
  "WT",
  expression(italic("\u0394suz12")),
  expression(italic("\u0394cac-1")),
  expression(italic("\u0394cac-2")),
  expression(italic("\u0394cac-1_2")),
  expression(italic("\u0394cac-1_suz12")),
  expression(italic("\u0394ash-1")),
  "3dpf",
  "4dpf",
  "5dpf",
  "6dpf"
)

#Group and summarize
total2 <- meltedAllData %>%
  group_by(Sample)

total_dist <- meltedAveragePRC2targetData %>%
  group_by(Sample) %>%
  summarise(num = n())

#Set dodge for violin spacing
dodge <- position_dodge(width = 1)

#Build plot
violin <- total2 %>%
  left_join(total_dist) %>%
  arrange(factor(Sample, levels = label_order)) %>%
  mutate(Sample = factor(Sample, levels = label_order)) %>%
  ggplot(aes(x = Sample, y = Count)) +
  geom_violin(position = dodge, scale = "width", trim = FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
  scale_x_discrete(labels = formatted_labels) +
  labs(
    x = "Strain",
    y = expression("Expression Level (log"[2]~"(TPM+1))")) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black")
  )

#check it
print(violin)

#heatmap
GenesWithChanges <- subset(AVERAGE_Prc2targetTPM, (rowSums(AVERAGE_Prc2targetTPM) > 0))

#plot in the desired column order; did this by subsetting the dataset based on sample list 'altorder' above
heatmap <- pheatmap::pheatmap(GenesWithChanges, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                    cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                    legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

#order of heatmap
gene_order <- if (!is.null(heatmap$tree_row)) {
  rownames(GenesWithChanges)[heatmap$tree_row$order]
} else {
  rownames(GenesWithChanges)  # if cluster_rows = FALSE
}

#plot
heatmap_plot <- heatmap[[4]]

#save plot
ggsave(filename = "./CAF-1_antisense_map.pdf", plot = heatmap_plot, dpi=600, height=4, width=3)


###Pull out genes of interest (higher int cac mutants)
# 1) Make the same row-z matrix pheatmap used
mat <- as.matrix(GenesWithChanges)              # the exact matrix you plotted
z   <- t(scale(t(mat), center = TRUE, scale = TRUE))

# (optional) drop rows with zero variance (sd = 0 → NA z-scores)
z <- z[complete.cases(z), , drop = FALSE]

# 2) For each gene (row), find which column has the highest z-score
max_col <- colnames(z)[max.col(z, ties.method = "first")]  # argmax per row

# 3) Keep rows whose max is in your target conditions
targets <- c("cac-1", "cac-2", "cac-1_2", "cac-1_suz12")

# If your column names contain these as substrings (e.g. replicates), map them:
is_target <- Reduce(`|`, lapply(targets, \(p) grepl(p, colnames(z), fixed = TRUE)))
target_cols <- colnames(z)[is_target]

rows_keep <- rownames(z)[max_col %in% target_cols]
subset_df <- GenesWithChanges[rows_keep, , drop = FALSE]    # original values
subset_z  <- z[rows_keep, , drop = FALSE]                   # z-scores (if you prefer)

heatmap <- pheatmap::pheatmap(subset_z, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
                              cellwidth = NA, cellheight = NA, scale = "row", cluster_rows = T, cluster_cols = F, clustering_method="centroid", clustering_distance_cols="euclidean",
                              legend=T, show_rownames=F, show_colnames=T, fontsize_col=10, treeheight_row=0, treeheight_col=5, height = 1.5, width = 2.5)

gene_order <- if (!is.null(heatmap$tree_row)) {
  rownames(GenesWithChanges)[heatmap$tree_row$order]
} else {
  rownames(GenesWithChanges)  # if cluster_rows = FALSE
}

gene_order = as.matrix(gene_order)
