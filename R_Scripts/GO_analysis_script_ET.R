#GO Analysis
```{r, my script}
# GO BP enrichment + lollipop plot for Δcac-1 & Δcac-2 DE genes
# Requires:
#   - Nc_GO_mapping.txt  (your FungiDB export)
#   - cac1_cac2__alpha0.05__DE_intersection.csv  (shared DE genes)

## --- packages ---
# cran_pkgs <- c("readr","dplyr","tidyr","ggplot2")
# miss <- setdiff(cran_pkgs, rownames(installed.packages()))
# if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", repos = "https://cloud.r-project.org")
# BiocManager::install(c("clusterProfiler","GO.db","AnnotationDbi"), ask = FALSE, update = FALSE)
# 
# suppressPackageStartupMessages({
#   library(readr); library(dplyr); library(tidyr); library(ggplot2)
#   library(clusterProfiler); library(GO.db); library(AnnotationDbi)
# })

## --- paths / outputs ---
fungidb_go_path <- "Nc_GO_mapping.txt"
de_csv_path     <- "../CAF-1_RNA-seq_Analysis/csv_files/cac1_cac2__alpha0.05__DE_intersection.csv"
out_tag  <- "CAF1_shared_GO"
out_res  <- paste0(out_tag, "_BP_enrichment.tsv")
out_plot <- paste0(out_tag, "_BP_lollipop.pdf")
top_n_terms <- 20

## --- read GO mapping (FungiDB tab with separate ID/name cols) ---
stopifnot(file.exists(fungidb_go_path))
go_raw <- readr::read_tsv(fungidb_go_path, show_col_types = FALSE)

gene_col <- intersect(names(go_raw), c("Gene ID","gene_id","Genes","ID"))[1]
if (is.na(gene_col)) stop("Could not find a 'Gene ID' column in Nc_GO_mapping.txt")

pairs <- list(
  list(id="Curated GO Process IDs",   nm="Curated GO Processes",   ns="BP"),
  list(id="Computed GO Process IDs",  nm="Computed GO Processes",   ns="BP"),
  list(id="Curated GO Component IDs", nm="Curated GO Components",   ns="CC"),
  list(id="Computed GO Component IDs",nm="Computed GO Components",  ns="CC"),
  list(id="Curated GO Function IDs",  nm="Curated GO Functions",    ns="MF"),
  list(id="Computed GO Function IDs", nm="Computed GO Functions",   ns="MF")
)

make_long <- function(idcol, nmcol, ns){
  if (!(idcol %in% names(go_raw))) return(NULL)
  tmp <- go_raw %>%
    dplyr::transmute(
      gene_id = .data[[gene_col]],
      go_id   = as.character(.data[[idcol]]),
      go_name = if (nmcol %in% names(go_raw)) as.character(.data[[nmcol]]) else NA_character_
    ) %>%
    dplyr::filter(!is.na(go_id), nzchar(go_id)) %>%
    tidyr::separate_rows(go_id,  sep = "\\s*;\\s*") %>%
    tidyr::separate_rows(go_name, sep = "\\s*;\\s*") %>%
    dplyr::mutate(namespace = ns) %>%
    dplyr::distinct()
  tmp
}

go_map <- do.call(
  dplyr::bind_rows,
  lapply(pairs, function(p) make_long(p$id, p$nm, p$ns))
) %>%
  dplyr::filter(grepl("^GO:", go_id)) %>%
  dplyr::distinct()

# fill missing names from GO.db
need_names <- setdiff(unique(go_map$go_id), unique(go_map$go_id[!is.na(go_map$go_name) & nzchar(go_map$go_name)]))
if (length(need_names)) {
  dbn <- AnnotationDbi::select(GO.db, keys = need_names, keytype = "GOID", columns = "TERM")
  dbn <- dplyr::rename(dbn, go_id = GOID, go_name_db = TERM)
  go_map <- go_map %>%
    dplyr::left_join(dbn, by = "go_id") %>%
    dplyr::mutate(go_name = dplyr::coalesce(go_name, go_name_db)) %>%
    dplyr::select(-go_name_db)
}

## --- read shared DE genes, auto-detect ID column ---
stopifnot(file.exists(de_csv_path))
de_df <- readr::read_csv(de_csv_path, show_col_types = FALSE)

guess_id_col <- function(df){
  # prefer obvious columns
  cand <- intersect(names(df), c("gene_id","Gene ID","ID","Genes","NCU","NCU_ID","NCUid","Attributes","VALUE"))
  if (length(cand)) return(cand[1])
  # else pick column with most NCU-like entries
  score <- sapply(df, function(x) sum(grepl("^NCU\\d{5}$", as.character(x))))
  names(which.max(score))
}
de_col <- guess_id_col(de_df)
genes_shared <- unique(as.character(de_df[[de_col]]))

## --- universe: all genes that have any GO term (quick + conservative) ---
universe_genes <- unique(as.character(go_map$gene_id))

## --- TERM2GENE / TERM2NAME for Biological Process only ---
TERM2GENE <- go_map %>%
  dplyr::filter(tolower(namespace) %in% c("bp","biological process","biological_process","p")) %>%
  dplyr::select(go_id, gene_id) %>% dplyr::distinct()
TERM2NAME <- go_map %>%
  dplyr::select(go_id, go_name) %>% dplyr::distinct()

## --- enrichment (BH FDR; we’ll filter to FDR <= 0.05 for plotting) ---
enr_bp <- clusterProfiler::enricher(
  gene          = genes_shared,
  universe      = universe_genes,
  TERM2GENE     = TERM2GENE,
  TERM2NAME     = TERM2NAME,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1.0,
  qvalueCutoff  = 1.0
)
if (is.null(enr_bp) || nrow(enr_bp@result) == 0)
  stop("No BP terms returned. Check IDs vs GO mapping.")

## --- compute fold-enrichment and tidy ---
res <- enr_bp@result %>%
  dplyr::mutate(
    GeneRatioNum = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    BgRatioNum   = sapply(strsplit(BgRatio,   "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    FoldEnrich   = GeneRatioNum / BgRatioNum,
    FDR          = p.adjust(pvalue, "BH"),
    minusLog10FDR = -log10(FDR),
    Term = dplyr::coalesce(Description, ID)   # Description is the term name from TERM2NAME
  ) %>%
  dplyr::arrange(dplyr::desc(FoldEnrich))

readr::write_tsv(res, out_res)

## --- lollipop plot: names on y, bigger size spread, vivid colors; FDR-filtered ---
max_chars <- 37
shorten <- function(x, n = max_chars) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)

plot_df <- res %>%
  dplyr::filter(FDR <= 0.05) %>%                # significant terms only
  dplyr::slice_head(n = top_n_terms) %>%
  dplyr::mutate(
    LabelShow = shorten(Term),                  # what you SEE (no GO IDs)
    LabelKey  = paste0(Term, "\u200B", dplyr::row_number()), # invisible uniqueness
    Label     = factor(LabelKey, levels = rev(LabelKey))
  )

lab_map <- setNames(plot_df$LabelShow, plot_df$LabelKey)

p <- ggplot(plot_df, aes(FoldEnrich, Label)) +
  geom_segment(aes(x = 0, xend = FoldEnrich, yend = Label), linewidth = 0.8) +
  geom_point(aes(size = Count, color = minusLog10FDR)) +
  scale_size(range = c(3, 10), name = "N. of Genes") +
  scale_color_viridis_c(option = "plasma", name = expression(-log[10]("FDR"))) +
  scale_y_discrete(labels = lab_map) +          # show shortened names only
  labs(x = "Fold Enrichment", y = NULL,
       title = expression("GO Biological Process Enrichment (cac-1 cac-2 DE)")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 9),
    legend.text  = element_text(size = 9),
    legend.title = element_text(size = 10),
    plot.title   = element_text(size = 12, hjust = 0.5),
    panel.grid.major.y = element_blank()
  )

ggplot2::ggsave(out_plot, p, width = 7, height = 5)
p
