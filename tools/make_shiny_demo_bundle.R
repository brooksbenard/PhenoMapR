#!/usr/bin/env Rscript
# Build inst/extdata/shiny/PAAD_CRA001160_demo_5000.rds for the Shiny demo.
# Run from package root: Rscript tools/make_shiny_demo_bundle.R
#
# Requires CRA001160 vignette files or a full Seurat RDS (see load_shiny_demo_pool).

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

dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

n_cells <- 5000L
n_genes <- 500L
bundle_seed <- 20250625L

pool <- load_shiny_demo_pool()
if (is.null(pool)) {
  stop(
    "CRA001160 source data not found. Place vignette files under vignettes/ ",
    "or set PHENOMAPR_CRA001160_H5 / PHENOMAPR_CRA001160_META."
  )
}

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
