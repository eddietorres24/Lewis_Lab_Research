#!/bin/bash
#SBATCH --job-name=rnaseq_qc_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --output=rnaseq_qc_all.%j.out
#SBATCH --error=rnaseq_qc_all.%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# ---------- user-editable inputs (relative to this folder) ----------
export COUNTS_FILE="readcounts_All_CAF1paper.merged.txt"
export COLDATA_FILE="coldata.csv"
export GENE_MAP_FILE="gene_annotation.csv"                 # optional
export K27_BED_FILE="H3K27me3_methylated_genes_FINAL.bed"  # optional
export OUTDIR="RNAseq_QC_out"
export REFERENCE_LEVEL="WT"
export MIN_COUNT_FILTER=10
# -------------------------------------------------------------------

# ---- Load Miniconda and ensure an environment with Bioc pkgs exists (Sapelo2 supports conda) ----
module load Miniconda3 || module load Anaconda3 || true   # Sapelo2 docs use Miniconda3
# Make env name stable
ENV_NAME="rnaseq-qc"

# Create env once if missing (R 4.4 works well with DESeq2 on conda)
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    r-base=4.4 \
    bioconductor-deseq2 \
    bioconductor-apeglm \
    bioconductor-summarizedexperiment \
    bioconductor-delayedarray \
    bioconductor-s4arrays \
    bioconductor-sparsearray \
    r-pheatmap r-data.table r-tidyverse r-ggrepel
fi

# Activate env
source activate "$ENV_NAME" 2>/dev/null || conda activate "$ENV_NAME"

# ------------------------------- R ---------------------------------
Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(dplyr); library(tibble)
  library(ggplot2); library(pheatmap); library(readr); library(stringr); library(tidyr)
})
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
have_apeglm  <- requireNamespace("apeglm",  quietly = TRUE)

counts_file      <- Sys.getenv("COUNTS_FILE")
coldata_file     <- Sys.getenv("COLDATA_FILE")
gene_map_file    <- Sys.getenv("GENE_MAP_FILE")
k27_bed_file     <- Sys.getenv("K27_BED_FILE")
outdir           <- Sys.getenv("OUTDIR")
reference_level  <- Sys.getenv("REFERENCE_LEVEL")
min_count_filter <- as.integer(Sys.getenv("MIN_COUNT_FILTER"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
figdir <- file.path(outdir, "figures"); dir.create(figdir, FALSE, TRUE)
tbldir <- file.path(outdir, "tables");  dir.create(tbldir, FALSE, TRUE)
qcdir  <- file.path(outdir, "qc");      dir.create(qcdir,  FALSE, TRUE)

message("\n=== Loading counts ===")
counts_dt <- data.table::fread(counts_file)
stopifnot(ncol(counts_dt) >= 2)
gene_id_col <- names(counts_dt)[1]
counts_df <- as.data.frame(counts_dt)
rownames(counts_df) <- counts_df[[gene_id_col]]
counts_df[[gene_id_col]] <- NULL
counts_df[] <- lapply(counts_df, function(x) as.integer(round(as.numeric(x))))
stopifnot(all(is.finite(as.matrix(counts_df))))

message("=== Loading coldata ===")
coldata <- read.csv(coldata_file, stringsAsFactors = FALSE)
stopifnot("sample_id" %in% names(coldata))
rownames(coldata) <- coldata$sample_id
miss_counts <- setdiff(coldata$sample_id, colnames(counts_df))
miss_col    <- setdiff(colnames(counts_df), coldata$sample_id)
if (length(miss_counts)) stop("Samples in coldata but not counts: ", paste(miss_counts, collapse=", "))
if (length(miss_col))    stop("Samples in counts but not coldata: ", paste(miss_col, collapse=", "))
coldata <- coldata[colnames(counts_df), , drop=FALSE]

if ("condition" %in% names(coldata)) {
  coldata$condition <- factor(coldata$condition)
  if (reference_level %in% levels(coldata$condition)) {
    coldata$condition <- stats::relevel(coldata$condition, ref = reference_level)
  }
} else warning("'condition' missing; using ~1 design (no contrasts).")

gene_map <- NULL
if (file.exists(gene_map_file)) {
  gm <- read.csv(gene_map_file, stringsAsFactors = FALSE)
  if (all(c("gene_id","gene_symbol") %in% names(gm))) {
    gene_map <- gm[!duplicated(gm$gene_id), c("gene_id","gene_symbol")]
  }
}

k27_ids <- NULL
if (file.exists(k27_bed_file)) {
  bed <- tryCatch(read.delim(k27_bed_file, header = FALSE, stringsAsFactors = FALSE), error=function(e) NULL)
  if (!is.null(bed)) {
    k27_ids <- unique(if (ncol(bed) >= 10) bed[[10]] else bed[[ncol(bed)]])
    k27_ids <- k27_ids[!is.na(k27_ids)]
  }
}

message("=== Build DESeq2 & QC ===")
keep <- rowSums(counts_df) >= min_count_filter
counts_f <- counts_df[keep, , drop=FALSE]
design_formula <- if ("condition" %in% names(coldata)) ~ condition else ~ 1
dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata, design = design_formula)
dds <- estimateSizeFactors(dds)

# sample-level QC
lib_sizes <- colSums(counts_df)
detected  <- colSums(counts_df > 0)
pct_zero  <- colSums(counts_df == 0) / nrow(counts_df) * 100
qc_tbl <- tibble(sample_id = names(lib_sizes),
                 library_size = as.numeric(lib_sizes),
                 detected_genes = as.numeric(detected),
                 pct_zero = round(pct_zero, 2),
                 size_factor = as.numeric(sizeFactors(dds)))
write.csv(qc_tbl, file.path(qcdir, "sample_qc_metrics.csv"), row.names = FALSE)

p1 <- ggplot(qc_tbl, aes(reorder(sample_id, library_size), library_size, fill = coldata$condition))+
  geom_col()+coord_flip()+labs(x="Sample", y="Library size", fill="condition")+theme_minimal(11)
ggsave(file.path(figdir,"library_sizes.pdf"), p1, width=7.5, height=6, useDingbats=FALSE)

p2 <- ggplot(qc_tbl, aes(reorder(sample_id, detected_genes), detected_genes, fill = coldata$condition))+
  geom_col()+coord_flip()+labs(x="Sample", y="Detected genes (>0)", fill="condition")+theme_minimal(11)
ggsave(file.path(figdir,"detected_genes.pdf"), p2, width=7.5, height=6, useDingbats=FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(tbldir, "normalized_counts.csv"))
write.csv(as.data.frame(counts(dds)), file.path(tbldir, "raw_counts_filtered.csv"))
write.csv(tibble(sample_id = names(sizeFactors(dds)), size_factor = sizeFactors(dds)),
          file.path(tbldir, "size_factors.csv"), row.names = FALSE)

log_norm <- log10(norm_counts + 1)
ln_df <- as.data.frame(log_norm) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to="sample_id", values_to="log10_norm")
ln_df$condition <- coldata[ln_df$sample_id, "condition", drop=TRUE]
p3 <- ggplot(ln_df, aes(sample_id, log10_norm, fill = condition))+
  geom_violin(scale="width", trim=TRUE)+coord_flip()+
  labs(x="Sample", y="log10(normalized counts + 1)", fill="condition")+theme_minimal(11)
ggsave(file.path(figdir,"normalized_count_violin.pdf"), p3, width=7.5, height=6.5, useDingbats=FALSE)

message("=== VST, PCA, correlation ===")
vsd <- vst(dds, blind = FALSE)
write.csv(as.data.frame(assay(vsd)), file.path(tbldir, "vst_matrix.csv"))

pca_df <- DESeq2::plotPCA(vsd, intgroup = intersect(c("condition","type","batch"), colnames(coldata)), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = name)) + geom_point(size = 3)
if (have_ggrepel) p_pca <- p_pca + ggrepel::geom_text_repel(size = 3, max.overlaps = 50)
p_pca <- p_pca + labs(x=paste0("PC1: ", percentVar[1], "%"), y=paste0("PC2: ", percentVar[2], "%"), color="condition") + theme_minimal(12)
ggsave(file.path(figdir,"PCA_vst.pdf"), p_pca, width=6.5, height=5.5, useDingbats=FALSE)

cm <- cor(assay(vsd), method="pearson")
ann <- NULL
if ("condition" %in% names(coldata)) { ann <- data.frame(condition=coldata$condition); rownames(ann) <- rownames(coldata) }
pheatmap(cm, annotation_col=ann, annotation_row=ann, display_numbers=TRUE, number_format="%.2f",
         fontsize_number=7, main="Sample–sample correlation (Pearson, VST)",
         filename=file.path(figdir,"correlation_heatmap_vst.pdf"), width=7.5, height=7)

message("=== DESeq2 fit (for dispersion, MA, p-hist) ===")
dds <- DESeq(dds)

pdf(file.path(figdir,"dispersion_plot.pdf"), width=6, height=5, useDingbats=TRUE); plotDispEsts(dds); dev.off()

ref <- if ("condition" %in% names(coldata)) levels(coldata$condition)[1] else NA_character_
if (!is.na(ref)) {
  for (lvl in setdiff(levels(coldata$condition), ref)) {
    res <- results(dds, contrast = c("condition", lvl, ref), alpha = 0.05)
    res_df <- as.data.frame(res) |> rownames_to_column("gene_id")

    if (have_apeglm) {
      coef_name <- paste0("condition_", lvl, "_vs_", ref)
      res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")
      res_shr_df <- as.data.frame(res_shr) |> rownames_to_column("gene_id")
      names(res_shr_df)[names(res_shr_df)=="log2FoldChange"] <- "log2FoldChange_shrunk"
      res_df <- dplyr::left_join(res_df, res_shr_df[, c("gene_id","log2FoldChange_shrunk")], by="gene_id")
    } else res_df$log2FoldChange_shrunk <- NA_real_

    if (!is.null(gene_map)) res_df <- dplyr::left_join(res_df, gene_map, by="gene_id")
    write.csv(res_df, file.path(tbldir, paste0("DE_results_", lvl, "_vs_", ref, ".csv")), row.names = FALSE)

    pdf(file.path(figdir, paste0("MA_", lvl, "_vs_", ref, ".pdf")), width=6.2, height=5.2, useDingbats=TRUE)
    plotMA(res, ylim=c(-6,6), main=paste(lvl,"vs",ref)); dev.off()
    if (have_apeglm) {
      pdf(file.path(figdir, paste0("MA_shrunk_", lvl, "_vs_", ref, ".pdf")), width=6.2, height=5.2, useDingbats=TRUE)
      plotMA(res_shr, ylim=c(-6,6), main=paste(lvl,"vs",ref,"(shrunk)")); dev.off()
    }
    pdf(file.path(figdir, paste0("pvalue_hist_", lvl, "_vs_", ref, ".pdf")), width=5.5, height=4.5, useDingbats=TRUE)
    hist(res$pvalue, breaks=50, col="grey", main=paste("p-value distribution:", lvl, "vs", ref), xlab="p-value"); dev.off()
  }
}

writeLines(c(capture.output(sessionInfo())), con = file.path(outdir, "sessionInfo.txt"))
message("\nDone. Outputs in: ", normalizePath(outdir))
RSCRIPT
