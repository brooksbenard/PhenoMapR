#!/usr/bin/env Rscript
## Smoke tests for the spatial-segmentation embedding option.
##
## The Visualization tab's Reduction dropdown now exposes
## "segmentation" as a sibling of "spatial" whenever the loaded
## object carries cell-segmentation boundaries:
##
##   * Seurat 5 FOV objects -- Xenium / CosMx / MERSCOPE / Visium HD
##     -- where any image's `Boundaries()` is non-empty.
##   * AnnData with `obsm["segmentation"]` (or one of the conventional
##     X_-prefixed variants).
##
## This file asserts source-level structure so it runs anywhere the
## package source is present, without needing a Seurat 5 install or
## a Python anndata environment in CI.

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

app_R <- file.path(dirname(helpers_R), "app.R")
if (!file.exists(app_R)) {
  app_R <- system.file("shiny", "app.R", package = "PhenoMapR")
}
stopifnot(file.exists(app_R))
app_src <- paste(readLines(app_R, warn = FALSE), collapse = "\n")

## --- 1. helpers.R declares the AnnData segmentation key whitelist --------
stopifnot_msg(
  grepl(
    "(?s)\\.anndata_segmentation_obsm_keys\\s*<-\\s*c\\([^)]*\"segmentation\"",
    helpers_src, perl = TRUE
  ),
  "AnnData segmentation obsm key whitelist exists"
)
for (k in c("X_segmentation", "X_segmentation_centroid")) {
  stopifnot_msg(
    grepl(sprintf("\"%s\"", k), helpers_src, fixed = TRUE),
    sprintf("whitelist includes \"%s\"", k)
  )
}

## --- 2. Seurat 5 boundaries detection helper exists ----------------------
stopifnot_msg(
  grepl("\\.has_segmentation_boundaries\\s*<-\\s*function",
        helpers_src, perl = TRUE),
  ".has_segmentation_boundaries() helper exists"
)
stopifnot_msg(
  grepl("Seurat::Boundaries\\(\\s*img\\s*\\)",
        helpers_src, perl = TRUE),
  "boundary detection delegates to Seurat::Boundaries()"
)

## --- 3. list_available_embeddings() advertises "segmentation" ------------
le_pos <- regexpr(
  "(?m)^list_available_embeddings\\s*<-\\s*function",
  helpers_src, perl = TRUE
)
stopifnot_msg(le_pos > 0L, "found list_available_embeddings() definition")
le_window <- substr(
  helpers_src, le_pos,
  min(nchar(helpers_src), le_pos + 4000L)
)
stopifnot_msg(
  grepl(
    "(?s)\\.has_segmentation_boundaries.{0,400}nms\\s*<-\\s*c\\(\\s*nms\\s*,\\s*\"segmentation\"\\s*\\)",
    le_window, perl = TRUE
  ),
  "Seurat path appends \"segmentation\" when boundaries are present"
)
stopifnot_msg(
  grepl(
    "(?s)seg_hit\\s*<-\\s*intersect\\(\\s*\\.anndata_segmentation_obsm_keys\\s*,\\s*keys\\s*\\)",
    le_window, perl = TRUE
  ),
  "AnnData path looks up segmentation obsm keys via the whitelist"
)
stopifnot_msg(
  grepl(
    "(?s)setdiff\\(\\s*keys\\s*,\\s*\\.anndata_segmentation_obsm_keys\\s*\\).{0,100}\"segmentation\"",
    le_window, perl = TRUE
  ),
  "AnnData path collapses raw obsm keys into a single \"segmentation\" entry"
)

## --- 4. extract_embedding() handles "segmentation" sibling branch --------
ee_pos <- regexpr(
  "(?m)^extract_embedding\\s*<-\\s*function",
  helpers_src, perl = TRUE
)
stopifnot_msg(ee_pos > 0L, "found extract_embedding() definition")
ee_window <- substr(
  helpers_src, ee_pos,
  min(nchar(helpers_src), ee_pos + 12000L)
)
stopifnot_msg(
  grepl(
    "is_segmentation\\s*<-\\s*identical\\(\\s*name\\s*,\\s*\"segmentation\"\\s*\\)",
    ee_window, perl = TRUE
  ),
  "extract_embedding() introduces an is_segmentation flag"
)
stopifnot_msg(
  grepl(
    'GetTissueCoordinates\\([^)]*which\\s*=\\s*"centroids"',
    ee_window, perl = TRUE
  ),
  "Seurat segmentation branch reads `which = \"centroids\"`"
)
stopifnot_msg(
  grepl(
    "is_spatial\\s*\\|\\|\\s*is_segmentation",
    ee_window, perl = TRUE
  ),
  "segmentation rides on the same tissue-frame layout as spatial"
)

## --- 5. UI copy in app.R mentions the new "segmentation" reduction ------
## The R source stores the embedded quotes as `\"...\"`, so when we read
## the file as plain text the bytes include the backslashes. We allow
## either escaped or unescaped variants here.
stopifnot_msg(
  grepl("(\\\\\"|\")segmentation(\\\\\"|\") reduction",
        app_src, perl = TRUE),
  "Embedding help <details> mentions the \"segmentation\" reduction"
)
stopifnot_msg(
  grepl("Xenium", app_src, fixed = TRUE) &&
    grepl("CosMx", app_src, fixed = TRUE),
  "help copy lists imaging-based platforms (Xenium / CosMx)"
)
stopifnot_msg(
  grepl("obsm\\['segmentation'\\]", app_src, perl = TRUE),
  "help copy advertises AnnData obsm['segmentation'] convention"
)

cat("\nAll spatial-segmentation smoke tests passed.\n")
