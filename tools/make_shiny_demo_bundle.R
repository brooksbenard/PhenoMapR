#!/usr/bin/env Rscript
# Build inst/extdata/shiny/PAAD_CRA001160_demo_5000.rds for the Shiny demo.
# Run from package root: Rscript tools/make_shiny_demo_bundle.R
#
# Loads the FULL CRA001160 cohort (TISCH2 H5 + metadata) from Google Drive
# (or reuses local copies under vignettes/), then writes a 5,000-cell subset
# for the Shiny app. The single-cell vignette uses the full cohort separately.

args <- commandArgs(trailingOnly = TRUE)
dest <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path("inst", "extdata", "shiny", "PAAD_CRA001160_demo_5000.rds")
}

helpers <- file.path("inst", "shiny", "helpers.R")
if (!file.exists(helpers)) {
  stop("Run from the PhenoMapR package root (inst/shiny/helpers.R not found).")
}
sys.source(helpers, envir = .GlobalEnv)

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is required to build the Shiny demo bundle.")
}
if (!requireNamespace("googledrive", quietly = TRUE)) {
  stop("googledrive is required to download CRA001160 from Google Drive.")
}

dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
vignette_dir <- if (dir.exists("vignettes")) "vignettes" else "."
h5_path <- file.path(vignette_dir, "PAAD_CRA001160_expression.h5")
meta_path <- file.path(vignette_dir, "PAAD_CRA001160_CellMetainfo_table.tsv")

options(googledrive_quiet = TRUE)
googledrive::drive_deauth()
drive_files <- list(
  list(path = h5_path, id = "1iFWJa13s5UClrP362CQtAAE7KYEo8iBc"),
  list(path = meta_path, id = "1yC7Vw3oQ2APB6ZlK7BUl-2DLKGJhFbGN")
)
for (f in drive_files) {
  if (!file.exists(f$path)) {
    message("Downloading ", basename(f$path), " from Google Drive...")
    googledrive::drive_download(
      googledrive::as_id(f$id), path = f$path, overwrite = TRUE
    )
  }
}

pool <- .load_shiny_demo_pool_from_h5_tsv(h5_path, meta_path)
if (is.null(pool)) {
  stop("Failed to load CRA001160 H5 + metadata for demo bundle generation.")
}
pool$source_info <- shiny_demo_source_info()

n_cells <- 5000L
n_genes <- 500L
bundle_seed <- 20250625L

set.seed(bundle_seed)
cells <- colnames(pool$expression)
keep <- sample(cells, min(n_cells, length(cells)), replace = FALSE)

expr_sub <- pool$expression[, keep, drop = FALSE]
gene_var <- .shiny_row_var(expr_sub)
top_genes <- order(gene_var, decreasing = TRUE)[seq_len(min(n_genes, length(gene_var)))]
expr <- if (inherits(expr_sub, "Matrix")) {
  as.matrix(expr_sub[top_genes, , drop = FALSE])
} else {
  expr_sub[top_genes, , drop = FALSE]
}

md <- pool$metadata[match(keep, pool$metadata$.cell_id), , drop = FALSE]
rownames(md) <- NULL

# Minor lineage is the default cell-type label in the Shiny app.
minor <- md$cell_type_minor
if (is.null(minor) || all(is.na(minor) | !nzchar(minor))) {
  minor <- md$cell_type_original
}
md$cell_type <- as.character(minor)

expr <- .attach_demo_matrix_obsm(expr, md)

info <- pool$source_info
info$n_cells_pool <- length(cells)
info$n_cells_sampled <- length(keep)
info$n_genes_sampled <- nrow(expr)
info$bundle_seed <- bundle_seed
info$is_presampled <- TRUE

bundle <- list(
  expression = expr,
  metadata = md,
  source_info = info,
  from = pool$from
)

saveRDS(bundle, dest, compress = "xz")
message("Wrote ", dest, " (", nrow(expr), " genes x ", ncol(expr), " cells).")
