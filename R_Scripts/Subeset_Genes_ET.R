#Code to subset genes that overlap with H3K27me3 marked regions

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# Load gene and peak data (assuming tab-delimited BED-like structure)
genes_df <- read.table("./bed_files/K27_genes_trimmed.bed", header = FALSE, stringsAsFactors = FALSE)
peaks_df <- read.table("./bed_files/WT_H3K27me3_Rep2_peaks.bed", header = FALSE, stringsAsFactors = FALSE)

# Create GRanges objects
genes_gr <- GRanges(seqnames = genes_promoter$V1,
                    ranges = IRanges(start = genes_promoter$V2 + 1, end = genes_promoter$V3),
                    strand = "*",
                    gene_id = genes_promoter$V11)

peaks_gr <- GRanges(seqnames = peaks_df$V1,
                    ranges = IRanges(start = peaks_df$V2 + 1, end = peaks_df$V3),
                    strand = "*")

# Find overlaps and calculate % of gene overlapped
hits <- findOverlaps(genes_gr, peaks_gr)
gene_hits <- genes_gr[queryHits(hits)]
peak_hits <- peaks_gr[subjectHits(hits)]
intersections <- pintersect(gene_hits, peak_hits)

# Summarize total overlap per gene
df_overlap <- data.frame(
  gene_id = mcols(gene_hits)$gene_id,
  gene_length = width(gene_hits),
  intersect_length = width(intersections)
) %>%
  group_by(gene_id, gene_length) %>%
  summarise(total_overlap = sum(intersect_length), .groups = "drop") %>%
  mutate(overlap_percent = (total_overlap / gene_length) * 100)

# Merge back with original gene info
genes_promoter$gene_id <- genes_promoter$V11
genes_annotated <- left_join(genes_promoter, df_overlap, by = "gene_id") %>%
  mutate(overlap_percent = ifelse(is.na(overlap_percent), 0, overlap_percent))

# Output BED files by tier
tiers <- seq(80, 100, by = 5)
for (tier in tiers) {
  subset_df <- genes_annotated %>%
    filter(overlap_percent >= tier) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11)
  
  filename <- paste0("genes_and_promoters_overlap_", tier, "pct.bed")
  write.table(subset_df, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#Make bed file with promoters included
genes_promoter <- genes_df

# Subtract 1000 from V2 if strand is "+"
genes_promoter$V2 <- ifelse(genes_promoter$V6 == "+", 
                            pmax(0, genes_promoter$V2 - 1000),  # avoid negative coords
                            genes_promoter$V2)

# Add 1000 to V3 if strand is "-"
genes_promoter$V3 <- ifelse(genes_promoter$V6 == "-", 
                            genes_promoter$V3 + 1000,
                            genes_promoter$V3)

###Make file w only promoters
# K27_prom = genes_promoter
# 
# K27_prom$V3 <- ifelse(K27_prom$V6 == "+", 
#                       pmax(0, K27_prom$V2 + 1000),  # avoid negative coords
#                       K27_prom$V3)
# 
# # Add 1000 to V3 if strand is "-"
# K27_prom$V2 <- ifelse(K27_prom$V6 == "-", 
#                       pmax(0, K27_prom$V3 - 1000),
#                       K27_prom$V2)
     