#!/usr/bin/env Rscript
## Smoke tests for the Python-free h5ad load path:
##
##   1. inst/shiny/helpers.R::.read_h5ad() prefers the pure-R hdf5r
##      reader over reticulate. Previously reticulate had priority,
##      so users with both packages installed landed in the AnnData-
##      Python branch and find_phenotype_markers() later errored
##      with "'expression' must be a matrix, Matrix, data.frame,
##      Seurat, or SingleCellExperiment".
##
##   2. .read_h5ad_hdf5r()'s "missing hdf5r" error message no longer
##      points users at Python (no anndata / py_install mention).
##
##   3. .read_h5ad_hdf5r() now reads /obs columns via hdf5r and
##      stashes them on the returned matrix's `phenomapr_obs`
##      attribute, so users get cell-level metadata without ever
##      touching Python. extract_object_metadata() recognises that
##      attribute and surfaces it as `.cell_id` + obs columns.
##
##   4. R/prognostic_analysis.R::process_expression_for_markers()
##      gained a `python.builtin.object` branch that pipes AnnData
##      input through .anndata_X_to_genes_cells(), so users who pass
##      a Python AnnData directly (e.g. from a reticulate session)
##      can find markers without the "must be a matrix..." error.
##
## We assert source-level shape rather than spinning up an actual
## h5ad fixture, so this smoke runs anywhere the package source is
## present without depending on hdf5r or reticulate.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

helpers_R <- system.file("shiny", "helpers.R", package = "PhenoMapR")
dev_dir   <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  helpers_R <- file.path(dev_dir, "helpers.R")
}
stopifnot(file.exists(helpers_R))
helpers_src <- paste(readLines(helpers_R, warn = FALSE), collapse = "\n")

prog_R <- file.path(dirname(helpers_R), "..", "..", "R", "prognostic_analysis.R")
if (!file.exists(prog_R)) {
  prog_R <- system.file("R", "prognostic_analysis.R", package = "PhenoMapR")
}
stopifnot(file.exists(prog_R))
prog_src <- paste(readLines(prog_R, warn = FALSE), collapse = "\n")

## --- 1. .read_h5ad() prefers hdf5r BEFORE reticulate ----------------------
## The dispatcher and the hdf5r helper share the same prefix, so anchor on
## a line-start `.read_h5ad <- function` (multiline mode) and a trailing
## word boundary so `.read_h5ad_hdf5r <- function` doesn't accidentally
## match.
read_pos <- regexpr(
  "(?m)^\\.read_h5ad\\b\\s*<-\\s*function",
  helpers_src, perl = TRUE
)
stopifnot_msg(read_pos > 0L, "found .read_h5ad() definition")
## Slice from the start of .read_h5ad through .read_h5ad_hdf5r so we
## measure ordering inside the dispatcher only (not the helper below it).
end_marker <- regexpr(
  "\\.read_h5ad_hdf5r\\s*<-\\s*function",
  helpers_src, perl = TRUE
)
stopifnot_msg(end_marker > 0L, "found .read_h5ad_hdf5r() definition")
dispatch_window <- substr(helpers_src, read_pos, end_marker - 1L)

hdf5r_check <- regexpr(
  'requireNamespace\\(\\s*"hdf5r"',
  dispatch_window, perl = TRUE
)
retic_check <- regexpr(
  'requireNamespace\\(\\s*"reticulate"',
  dispatch_window, perl = TRUE
)
stopifnot_msg(hdf5r_check > 0L,
              "dispatcher checks requireNamespace(\"hdf5r\")")
stopifnot_msg(retic_check > 0L,
              "dispatcher checks requireNamespace(\"reticulate\")")
stopifnot_msg(hdf5r_check < retic_check,
              "hdf5r check is BEFORE reticulate check (Python-free first)")

## --- 2. .read_h5ad_hdf5r() error message is Python-free ------------------
hdf5r_def_pos <- end_marker
hdf5r_window <- substr(
  helpers_src, hdf5r_def_pos,
  min(nchar(helpers_src), hdf5r_def_pos + 6000L)
)
stopifnot_msg(
  grepl("Reading \\.h5ad files requires the `hdf5r` R package",
        hdf5r_window, perl = TRUE),
  "missing-hdf5r error advertises a pure-R install path"
)
stopifnot_msg(
  !grepl("py_install", hdf5r_window, fixed = TRUE),
  "missing-hdf5r error no longer suggests reticulate::py_install()"
)
stopifnot_msg(
  !grepl("Python `anndata`", hdf5r_window, fixed = TRUE),
  "missing-hdf5r error no longer suggests Python anndata"
)

## --- 3. .read_h5ad_hdf5r() reads /obs and stashes phenomapr_obs ----------
stopifnot_msg(
  grepl('\\.read_h5ad_obs_hdf5r\\s*<-\\s*function',
        helpers_src, perl = TRUE),
  ".read_h5ad_obs_hdf5r() helper exists"
)
stopifnot_msg(
  grepl('attr\\(\\s*X_genes_x_cells\\s*,\\s*"phenomapr_obs"\\s*\\)\\s*<-',
        helpers_src, perl = TRUE),
  ".read_h5ad_hdf5r() stamps the matrix with attr(\"phenomapr_obs\")"
)
stopifnot_msg(
  grepl('attr\\(\\s*obj\\s*,\\s*"phenomapr_obs"\\s*\\)',
        helpers_src, perl = TRUE),
  "extract_object_metadata() reads attr(obj, \"phenomapr_obs\")"
)

## --- 4. process_expression_for_markers() has a python.builtin.object br --
stopifnot_msg(
  grepl('inherits\\(\\s*expression\\s*,\\s*"python\\.builtin\\.object"\\s*\\)',
        prog_src, perl = TRUE),
  'process_expression_for_markers() handles python.builtin.object'
)
stopifnot_msg(
  grepl('\\.anndata_X_to_genes_cells\\(\\s*expression\\s*\\)',
        prog_src, perl = TRUE),
  "AnnData branch pipes through .anndata_X_to_genes_cells()"
)
stopifnot_msg(
  grepl(
    "(?s)must be a matrix, Matrix, data\\.frame, Seurat,.{0,80}SingleCellExperiment, or AnnData",
    prog_src, perl = TRUE
  ),
  "catch-all stop() now mentions AnnData as a supported type"
)

cat("\nAll h5ad-Python-free-load smoke tests passed.\n")
