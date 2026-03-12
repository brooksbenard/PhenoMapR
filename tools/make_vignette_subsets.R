#!/usr/bin/env Rscript
# Regenerate vignette subset RDS files (max 5000 cells per cell type, no types lost).
# Run from package root: Rscript tools/make_vignette_subsets.R
# Requires: googledrive, Seurat. Full datasets are downloaded from Google Drive.

args <- commandArgs(trailingOnly = TRUE)
dest_dir <- if (length(args) >= 1) args[[1]] else file.path("inst", "extdata", "vignette_subsets")

dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

download_from_drive <- function(file_id, dest) {
  if (file.exists(dest)) return(invisible(dest))
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    stop("googledrive is required. Install with install.packages('googledrive').")
  }
  googledrive::drive_deauth()
  googledrive::drive_download(googledrive::as_id(file_id), path = dest, overwrite = TRUE)
  invisible(dest)
}

find_celltype_col <- function(obj) {
  md <- obj@meta.data
  candidates <- c(
    "Celltype..major.lineage.",
    "cell_type", "CellType", "celltype",
    "annotation", "Annotation",
    "cluster", "seurat_clusters"
  )
  for (col in candidates) {
    if (col %in% names(md) && (is.character(md[[col]]) || is.factor(md[[col]])) &&
        length(unique(md[[col]])) >= 1) return(col)
  }
  NULL
}

subset_cells_by_type <- function(obj, max_per_type = 5000L, celltype_col = NULL) {
  if (is.null(celltype_col)) celltype_col <- find_celltype_col(obj)
  if (is.null(celltype_col)) {
    n <- min(5000L, ncol(obj))
    return(colnames(obj)[seq_len(n)])
  }
  grp <- as.character(obj@meta.data[[celltype_col]])
  names(grp) <- colnames(obj)
  keep <- character(0)
  set.seed(1)
  for (g in unique(grp)) {
    cells <- names(grp)[grp == g]
    n <- min(as.integer(max_per_type), length(cells))
    keep <- c(keep, if (length(cells) <= n) cells else sample(cells, n))
  }
  keep
}

save_seurat_subset <- function(src_rds, dest_rds, max_per_type = 5000L, n_genes = 1500L) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  obj <- readRDS(src_rds)
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object in: ", src_rds)
  obj <- tryCatch(Seurat::UpdateSeuratObject(obj), error = function(e) obj)
  keep_cells <- subset_cells_by_type(obj, max_per_type = max_per_type)
  obj <- subset(obj, cells = keep_cells)
  assay <- if ("RNA" %in% names(obj@assays)) "RNA" else names(obj@assays)[1]
  Seurat::DefaultAssay(obj) <- assay
  feats <- Seurat::VariableFeatures(obj)
  if (length(feats) < 10) feats <- rownames(obj[[assay]])
  n_genes <- min(as.integer(n_genes), length(feats))
  obj <- subset(obj, features = feats[seq_len(n_genes)])
  obj <- Seurat::DietSeurat(obj, assays = assay, counts = TRUE, data = TRUE,
    scale.data = FALSE, dimreducs = intersect(c("pca", "umap", "tsne"), names(obj@reductions)),
    graphs = character(0), misc = FALSE)
  saveRDS(obj, dest_rds, compress = "xz")
  invisible(dest_rds)
}

save_spatial_subset <- function(src_rds, dest_rds, max_per_type = 5000L, n_genes = 1500L) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  obj <- readRDS(src_rds)
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object in: ", src_rds)
  obj <- tryCatch(Seurat::UpdateSeuratObject(obj), error = function(e) obj)
  assay <- if ("Spatial" %in% names(obj@assays)) "Spatial" else if ("RNA" %in% names(obj@assays)) "RNA" else names(obj@assays)[1]
  Seurat::DefaultAssay(obj) <- assay
  keep_cells <- subset_cells_by_type(obj, max_per_type = max_per_type)
  if (length(keep_cells) == 0) keep_cells <- colnames(obj)[seq_len(min(3000L, ncol(obj)))]
  obj <- subset(obj, cells = keep_cells)
  feats <- Seurat::VariableFeatures(obj)
  if (length(feats) < 10) feats <- rownames(obj[[assay]])
  n_genes <- min(as.integer(n_genes), length(feats))
  obj <- subset(obj, features = feats[seq_len(n_genes)])
  saveRDS(obj, dest_rds, compress = "xz")
  invisible(dest_rds)
}

message("Writing vignette subsets to: ", dest_dir)
ids <- list(gse111672 = "1vJxIlW_kvqFXOPw9qn1L7D6QbMpe-P0c",
  cra001160 = "14p_fYIFeuuRdXBF3J-5ZsXElq_mduSzb",
  spatial = "1HM0dBrQnaNsdm5mnq23aaQ2ILofJ0_vj")
tmp_dir <- tempdir()
gse_full <- file.path(tmp_dir, "PAAD_GSE111672_seurat.rds")
cra_full <- file.path(tmp_dir, "PAAD_CRA001160_seurat.rds")
spatial_full <- file.path(tmp_dir, "HT270P1_processed.rds")
download_from_drive(ids$gse111672, gse_full)
download_from_drive(ids$cra001160, cra_full)
download_from_drive(ids$spatial, spatial_full)
save_seurat_subset(gse_full, file.path(dest_dir, "PAAD_GSE111672_seurat_subset.rds"), max_per_type = 5000, n_genes = 1500)
save_seurat_subset(cra_full, file.path(dest_dir, "PAAD_CRA001160_seurat_subset.rds"), max_per_type = 5000, n_genes = 1500)
save_spatial_subset(spatial_full, file.path(dest_dir, "HT270P1_processed_subset.rds"), max_per_type = 5000, n_genes = 1500)
message("Done.")
