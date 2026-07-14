#!/usr/bin/env Rscript
# Build the FULL CRA001160 vignette bundle used by pkgdown / CI.
# Requires local vignettes/PAAD_CRA001160_seurat.rds + CellMetainfo TSV.
suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload", call. = FALSE)
  }
  pkgload::load_all(quiet = TRUE, helpers = FALSE, export_all = FALSE)
  library(Seurat)
  library(Matrix)
})

seu_path <- if (file.exists("vignettes/PAAD_CRA001160_seurat.rds")) {
  "vignettes/PAAD_CRA001160_seurat.rds"
} else {
  "PAAD_CRA001160_seurat.rds"
}
meta_path <- if (file.exists("vignettes/PAAD_CRA001160_CellMetainfo_table.tsv")) {
  "vignettes/PAAD_CRA001160_CellMetainfo_table.tsv"
} else {
  "PAAD_CRA001160_CellMetainfo_table.tsv"
}
stopifnot(file.exists(seu_path), file.exists(meta_path))

seu <- readRDS(seu_path)
expr <- PhenoMapR:::.get_assay_data_compat(seu, assay = "RNA", slot = "counts")
meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
cells <- intersect(colnames(expr), meta$Cell)
expr <- expr[, cells, drop = FALSE]
meta <- meta[match(cells, meta$Cell), , drop = FALSE]
obj <- list(expression = expr, metadata = meta, from = "CRA001160_full")

outdir <- "/tmp/phenomapr_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, "PAAD_CRA001160_full.rds")
tmp <- paste0(out, ".tmp")
saveRDS(obj, tmp, compress = "gzip")
file.rename(tmp, out)

message("Wrote ", out, " (", round(file.info(out)$size / 1e6, 1), " MB)")
message("dims ", paste(dim(expr), collapse = "x"))
