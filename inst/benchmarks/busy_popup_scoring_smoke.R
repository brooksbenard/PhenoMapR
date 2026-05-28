#!/usr/bin/env Rscript
## Smoke tests for the four "post-upload" fixes:
##
##   1. JS handleShow() resets dismissedThisRun so a previously
##      X-dismissed popup does NOT silently suppress a fresh
##      R-driven show on the next upload / scoring run.
##   2. The expression-upload observer in inst/shiny/app.R calls
##      phenomapr_busy_hide() BEFORE touching state, then defers
##      the heavy state propagation via later::later() so the
##      busy-hide message ships to the client ahead of the
##      renderUI / renderPlot cascade.
##   3. The Scoring summary "Input data" row appends the
##      predicted single-cell vs bulk classification when the
##      uploaded expression is a plain matrix.
##   4. score_show_slot_block / score_have_matrix output flags
##      gate the slot/assay control block via conditionalPanel
##      (hidden for matrix uploads, visible otherwise).
##
## Mocks the relevant fragments of state, input, expr_diagnostics
## etc.  Avoids spinning up an actual Shiny session.

suppressPackageStartupMessages({
  library(shiny)
})

stopifnot(requireNamespace("PhenoMapR", quietly = TRUE))

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R     <- system.file("shiny", "app.R",     package = "PhenoMapR")
helpers_R <- system.file("shiny", "helpers.R", package = "PhenoMapR")
dev_dir   <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R     <- file.path(dev_dir, "app.R")
  helpers_R <- file.path(dev_dir, "helpers.R")
}
stopifnot(file.exists(app_R), file.exists(helpers_R))

app_src     <- paste(readLines(app_R, warn = FALSE),   collapse = "\n")
helpers_src <- paste(readLines(helpers_R, warn = FALSE), collapse = "\n")

## --- 1. JS reset of dismissedThisRun in handleShow -------------------------
stopifnot_msg(
  grepl("dismissedThisRun = false;\\s*currentMessage", helpers_src),
  "JS handleShow resets dismissedThisRun before rendering"
)
## And the *blocking* `if (dismissedThisRun) ... return;` branch must be gone
## from handleShow (it survived in earlier iterations and was the root cause
## of the "popup never reappears on the next upload" report).
hs_block <- regmatches(
  helpers_src,
  regexpr("function handleShow\\(p\\) \\{[\\s\\S]*?\\n  \\}", helpers_src,
          perl = TRUE)
)
stopifnot_msg(length(hs_block) == 1L,
              "located the handleShow function body in helpers.R")
stopifnot_msg(
  !grepl("if \\(dismissedThisRun\\)", hs_block),
  "handleShow no longer suppresses renderShow when dismissedThisRun is set"
)

## --- 2. Expression observer defers state assignments via later::later ------
## Locate the slice of app.R between phenomapr_busy_show("Reading expression
## file...") and the matching `}, ignoreNULL = FALSE, ignoreInit = TRUE)`
## terminator so we are inspecting *only* the expression-upload observer.
slice_match <- regexpr(
  "phenomapr_busy_show\\(\"Reading expression file[\\s\\S]*?ignoreNULL = FALSE",
  app_src, perl = TRUE
)
stopifnot_msg(slice_match > 0,
              "located the expression-upload observer body in app.R")
slice <- substr(app_src, slice_match,
                slice_match + attr(slice_match, "match.length") - 1L)

stopifnot_msg(
  grepl("phenomapr_busy_hide\\(\\)", slice),
  "observer calls phenomapr_busy_hide() explicitly after parse"
)
stopifnot_msg(
  grepl("later::later\\(function\\(\\)", slice),
  "observer defers state assignments via later::later(function() ...)"
)
stopifnot_msg(
  grepl("state\\$expression <- res_obj\\$object", slice),
  "deferred body assigns state$expression from captured res_obj"
)
stopifnot_msg(
  grepl("session = sess", slice),
  "deferred body passes session = sess to showNotification calls"
)
## The old `on.exit(phenomapr_busy_hide(), add = TRUE)` pattern must NOT be
## present in the expr-upload observer any more -- it would defeat the whole
## point of the defer (state assignments would re-queue the hide message
## behind the downstream cascade in the same flush).
stopifnot_msg(
  !grepl("on\\.exit\\(phenomapr_busy_hide", slice),
  "expr-upload observer no longer uses on.exit(phenomapr_busy_hide())"
)

## Same checks for the metadata-upload observer.
meta_match <- regexpr(
  "phenomapr_busy_show\\(\"Loading cell metadata[\\s\\S]*?ignoreNULL = FALSE",
  app_src, perl = TRUE
)
stopifnot_msg(meta_match > 0,
              "located the metadata-upload observer body in app.R")
meta_slice <- substr(app_src, meta_match,
                     meta_match + attr(meta_match, "match.length") - 1L)
stopifnot_msg(
  grepl("later::later\\(function\\(\\)", meta_slice),
  "metadata observer also defers via later::later"
)
stopifnot_msg(
  !grepl("on\\.exit\\(phenomapr_busy_hide", meta_slice),
  "metadata observer no longer uses on.exit(phenomapr_busy_hide())"
)

## --- 3. score_data_status surfaces sc/bulk for matrix inputs ---------------
## Re-create the matrix branch logic and confirm both the "predicted
## single-cell" and "predicted bulk" labels reach the rendered HTML.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

build_data_row <- function(s, diag = NULL) {
  kind_label <- switch(s$kind %||% "loaded",
    seurat  = "Seurat object", sce = "SingleCellExperiment",
    spatial = "Spatial dataset", matrix = "Expression matrix",
    anndata = "AnnData", sprintf("%s object", s$kind %||% "loaded"))
  bits <- character(0)
  if (!is.na(s$default_assay %||% NA_character_) && nzchar(s$default_assay)) {
    bits <- c(bits, sprintf("assay: %s", s$default_assay))
  }
  if (length(s$layers_avail %||% character(0))) {
    bits <- c(bits, sprintf("layers: %s",
                            paste(s$layers_avail, collapse = ", ")))
  }
  if (identical(s$kind %||% "", "matrix") && !is.null(diag) &&
      length(diag$sc_or_bulk %||% "")) {
    sb_label <- switch(diag$sc_or_bulk,
      single_cell = "predicted single-cell",
      bulk        = "predicted bulk",
      unclear     = "single-cell vs bulk unclear",
      diag$sc_or_bulk
    )
    if (nzchar(sb_label)) bits <- c(bits, sb_label)
  }
  tags$div(class = "sds-row sds-row-ok",
    tags$div(class = "sds-icon", icon("check-circle")),
    tags$div(class = "sds-content",
             tags$div(class = "sds-label", "Input data"),
             tags$div(class = "sds-value", kind_label),
             if (length(bits))
               tags$div(class = "sds-detail",
                        paste(bits, collapse = "  \u00B7  "))))
}

## Matrix + single-cell prediction.
h <- as.character(build_data_row(
  list(kind = "matrix"),
  list(sc_or_bulk = "single_cell")
))
stopifnot_msg(grepl("Expression matrix", h),
              "matrix data-row labels the input correctly")
stopifnot_msg(grepl("predicted single-cell", h, fixed = TRUE),
              "matrix data-row surfaces 'predicted single-cell'")

## Matrix + bulk prediction.
h <- as.character(build_data_row(
  list(kind = "matrix"),
  list(sc_or_bulk = "bulk")
))
stopifnot_msg(grepl("predicted bulk", h, fixed = TRUE),
              "matrix data-row surfaces 'predicted bulk'")

## Matrix + unclear prediction.
h <- as.character(build_data_row(
  list(kind = "matrix"),
  list(sc_or_bulk = "unclear")
))
stopifnot_msg(grepl("single-cell vs bulk unclear", h, fixed = TRUE),
              "matrix data-row surfaces 'unclear' classification verbatim")

## Object input must NOT include any predicted bulk/sc label (that's for
## matrices only -- objects already report assay + layers).
h <- as.character(build_data_row(
  list(kind = "seurat", default_assay = "RNA",
       layers_avail = c("data", "counts")),
  list(sc_or_bulk = "single_cell")
))
stopifnot_msg(!grepl("predicted single-cell", h, fixed = TRUE),
              "Seurat/object data-row does NOT inject the sc/bulk label")
stopifnot_msg(grepl("assay: RNA", h, fixed = TRUE) &&
              grepl("layers: data, counts", h, fixed = TRUE),
              "object data-row still surfaces assay + layers")

## --- 4. score_show_slot_block / score_have_matrix flags + conditionalPanel -
stopifnot_msg(
  grepl("output\\$score_show_slot_block <- reactive", app_src),
  "score_show_slot_block reactive is defined server-side"
)
stopifnot_msg(
  grepl("output\\$score_have_matrix <- reactive", app_src),
  "score_have_matrix reactive is defined server-side"
)
stopifnot_msg(
  grepl('outputOptions(output, "score_show_slot_block", suspendWhenHidden = FALSE)',
        app_src, fixed = TRUE),
  "score_show_slot_block is kept live via suspendWhenHidden = FALSE"
)
stopifnot_msg(
  grepl('outputOptions(output, "score_have_matrix", suspendWhenHidden = FALSE)',
        app_src, fixed = TRUE),
  "score_have_matrix is kept live via suspendWhenHidden = FALSE"
)
stopifnot_msg(
  grepl("conditionalPanel\\(\\s*condition = \"output\\.score_show_slot_block\"",
        app_src, perl = TRUE),
  "slot/assay block is wrapped in conditionalPanel('output.score_show_slot_block')"
)
stopifnot_msg(
  grepl("output.score_have_matrix", app_src, fixed = TRUE),
  "matrix-only note panel is keyed on output.score_have_matrix"
)
stopifnot_msg(
  grepl("Plain matrices are scored directly", app_src),
  "matrix-only note panel surfaces the 'matrices scored directly' helpText"
)

## Exercise the *logic* of the two flags directly.
expr_summaries <- list(
  no_upload = NULL,
  matrix    = list(kind = "matrix"),
  seurat    = list(kind = "seurat",  default_assay = "RNA",
                   layers_avail = c("data", "counts")),
  sce       = list(kind = "sce",     default_assay = "logcounts",
                   layers_avail = c("counts", "logcounts")),
  spatial   = list(kind = "spatial", default_assay = "Spatial",
                   layers_avail = c("data", "counts")),
  anndata   = list(kind = "anndata", default_assay = NA_character_,
                   layers_avail = "X")
)

flag_show_slot   <- function(s) !is.null(s) && !identical(s$kind %||% "", "matrix")
flag_have_matrix <- function(s) !is.null(s) &&  identical(s$kind %||% "", "matrix")

stopifnot_msg(identical(flag_show_slot(expr_summaries$no_upload), FALSE),
              "no-upload  -> slot block hidden")
stopifnot_msg(identical(flag_have_matrix(expr_summaries$no_upload), FALSE),
              "no-upload  -> matrix note hidden")

stopifnot_msg(identical(flag_show_slot(expr_summaries$matrix), FALSE),
              "matrix     -> slot block hidden")
stopifnot_msg(identical(flag_have_matrix(expr_summaries$matrix), TRUE),
              "matrix     -> matrix note shown")

for (kn in c("seurat", "sce", "spatial", "anndata")) {
  stopifnot_msg(identical(flag_show_slot(expr_summaries[[kn]]), TRUE),
                sprintf("%-9s -> slot block shown", kn))
  stopifnot_msg(identical(flag_have_matrix(expr_summaries[[kn]]), FALSE),
                sprintf("%-9s -> matrix note hidden", kn))
}

cat("\nAll busy-popup + scoring-summary + slot-block smoke tests passed.\n")
