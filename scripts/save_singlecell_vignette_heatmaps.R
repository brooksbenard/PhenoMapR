#!/usr/bin/env Rscript
# Reproduce single-cell vignette steps through marker heatmaps and save PDFs to Desktop.
# Run from package root: Rscript scripts/save_singlecell_vignette_heatmaps.R
# Requires: vignettes/PAAD_CRA001160_expression.h5 and ...CellMetainfo_table.tsv
#   (or googledrive downloads). PhenoMap scoring twice on full matrix can take many minutes.
#
# Usage:
#   cd /path/to/PhenoMapR && Rscript scripts/save_singlecell_vignette_heatmaps.R
#   Rscript scripts/save_singlecell_vignette_heatmaps.R /path/to/PhenoMapR

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[1L])
} else {
  normalizePath(".")
}

# Load package from source tree when run from repo (pkgload)
if (file.exists(file.path(root, "DESCRIPTION"))) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload or run after devtools::install()")
  }
  pkgload::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(PhenoMapR))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
})

vdir <- file.path(root, "vignettes")
h5 <- file.path(vdir, "PAAD_CRA001160_expression.h5")
tsv <- file.path(vdir, "PAAD_CRA001160_CellMetainfo_table.tsv")

if (!file.exists(h5) || !file.exists(tsv)) {
  message("Local data not found; attempting Google Drive download to vignettes/ ...")
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    stop("Install googledrive or place PAAD_CRA001160_expression.h5 and PAAD_CRA001160_CellMetainfo_table.tsv in vignettes/")
  }
  options(googledrive_quiet = TRUE)
  googledrive::drive_deauth()
  googledrive::drive_download(googledrive::as_id("1PolTXggREz8XmhutCLTQJGCfKxFAzqMl"), h5, overwrite = TRUE)
  googledrive::drive_download(googledrive::as_id("17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu"), tsv, overwrite = TRUE)
}

message("Reading expression matrix ...")
expr_mat <- Seurat::Read10X_h5(h5)

message("Reading metadata ...")
meta <- read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(`Celltype (original)` = gsub(" cell", "", `Celltype (original)`)) %>%
  dplyr::rename(
    celltype_malignancy = `Celltype (malignancy)`,
    celltype_original = `Celltype (original)`
  )

pal_cells <- PhenoMapR::get_celltype_palette(meta$celltype_original)

message("PhenoMap (TCGA PAAD) ...")
scores_tcga <- PhenoMap(expression = expr_mat, reference = "tcga", cancer_type = "PAAD", verbose = TRUE)
message("PhenoMap (PRECOG Pancreatic) ...")
scores_precog <- PhenoMap(expression = expr_mat, reference = "precog", cancer_type = "Pancreatic", verbose = TRUE)

for (col in names(scores_tcga)) meta[[col]] <- scores_tcga[meta$Cell, col]
for (col in names(scores_precog)) meta[[col]] <- scores_precog[meta$Cell, col]

score_tcga_col <- grep("weighted_sum_score.*PAAD", names(meta), value = TRUE, ignore.case = TRUE)[1]
score_precog_col <- grep("weighted_sum_score.*Pancreatic", names(meta), value = TRUE, ignore.case = TRUE)[1]
if (is.na(score_tcga_col)) score_tcga_col <- names(scores_tcga)[1]
if (is.na(score_precog_col)) score_precog_col <- names(scores_precog)[1]

scores_df <- meta[, c(score_tcga_col, score_precog_col), drop = FALSE]
groups <- define_phenotype_groups(scores_df, percentile = 0.05, score_columns = score_precog_col)
group_col <- grep("phenotype_group", names(groups), value = TRUE)[1]
meta[[group_col]] <- groups[rownames(meta), group_col]

message("find_phenotype_markers (global) ...")
markers <- find_phenotype_markers(
  expr_mat,
  group_labels = meta[[group_col]],
  pval_threshold = 0.05,
  max_cells_per_ident = 5000L
)

out_dir <- file.path(Sys.getenv("HOME"), "Desktop")
if (!dir.exists(out_dir)) {
  out_dir <- root
  message("Desktop not found; writing to: ", out_dir)
}

pdf1 <- file.path(out_dir, "single_cell_vignette_heatmap_global_markers.pdf")
message("Drawing global marker heatmap -> ", pdf1)
grDevices::pdf(pdf1, width = 12, height = 8)
plot_phenotype_markers(
  markers = markers,
  expr_mat = expr_mat,
  meta = meta,
  cell_id_col = "Cell",
  group_col = group_col,
  score_col = score_precog_col,
  celltype_col = "celltype_original",
  celltype_palette = pal_cells,
  heatmap_type = "global",
  top_n_markers = 20L,
  n_mark_labels = 5L,
  p_adj_threshold = 0.05
)
grDevices::dev.off()

message("find_phenotype_markers (cell_type_specific) ...")
group_df <- data.frame(
  cell_id = meta$Cell,
  phenotype_group = meta[[group_col]],
  cell_type = meta$celltype_original,
  stringsAsFactors = FALSE
)
markers_ct <- find_phenotype_markers(
  expr_mat,
  group_labels = group_df,
  group_column = "phenotype_group",
  cell_id_column = "cell_id",
  cell_type_column = "cell_type",
  marker_scope = "cell_type_specific",
  pval_threshold = 0.05,
  max_cells_per_ident = 5000L,
  verbose = FALSE
)

pdf2 <- file.path(out_dir, "single_cell_vignette_heatmap_celltype_specific_markers.pdf")
message("Drawing cell-type-specific marker heatmap -> ", pdf2)
grDevices::pdf(pdf2, width = 14, height = 10)
plot_phenotype_markers(
  markers = markers_ct,
  expr_mat = expr_mat,
  meta = meta,
  cell_id_col = "Cell",
  group_col = group_col,
  score_col = score_precog_col,
  celltype_col = "celltype_original",
  celltype_palette = pal_cells,
  heatmap_type = "cell_type_specific",
  top_n_markers = 20L,
  n_mark_labels = 5L,
  p_adj_threshold = 0.05
)
grDevices::dev.off()

message("Done.")
message("  ", pdf1)
message("  ", pdf2)
