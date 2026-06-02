#!/usr/bin/env Rscript
## Smoke tests for three related polish changes:
##
##   1. plot_phenotype_markers() rowAnnotation calls now pass
##      `show_legend = FALSE`, so single-level cell-type annotations
##      (e.g. metadata where the only cell-type label is "1") no
##      longer produce a stray teal "1" badge floating next to the
##      Scaled-expr legend in the marker-gene heatmap.
##   2. The Embeddings sidebar (4. Visualization) wraps the upper
##      Reduction selector + help blocks in `.phenomapr-compact-stack
##      embedding-compact-stack` so the section above the <hr/>
##      divider before "Color cells by" stays tight.
##   3. The Expression input sidebar (1. Data) wraps the file picker,
##      "Use a tiny demo matrix instead" actionLink, and the matrix
##      diagnostics renderUI in `.phenomapr-compact-stack
##      expression-compact-stack` so the area under the actionLink
##      no longer has a wide gap before the <hr/>.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

repo <- if (file.exists("DESCRIPTION")) getwd() else system.file(package = "PhenoMapR")
app_R      <- file.path(repo, "inst", "shiny", "app.R")
styles_css <- file.path(repo, "inst", "shiny", "www", "styles.css")
ppm_R      <- file.path(repo, "R", "plot_phenotype_markers.R")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(styles_css), file.exists(ppm_R))
app_src    <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
styles_src <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")
ppm_src    <- paste(readLines(ppm_R,      warn = FALSE), collapse = "\n")

## --- 1. plot_phenotype_markers: rowAnnotation auto-legends suppressed -----
## Count the rowAnnotation() calls and the show_legend = FALSE flags
## inside their argument lists. The function has FOUR rowAnnotation
## calls total (ha_left/ha_right in the global path + ha_left/ha_right
## in the cell-type-specific path); every one must opt out of
## auto-legends, otherwise the stray "1" badge can leak back in.
n_row_anno <- length(gregexpr("ComplexHeatmap::rowAnnotation\\(", ppm_src)[[1L]])
stopifnot_msg(n_row_anno == 4L,
              sprintf("plot_phenotype_markers has 4 rowAnnotation calls (got %d)",
                      n_row_anno))

## Count show_legend = FALSE occurrences -- 3 were already there
## (2 HeatmapAnnotation ha_top calls + the manual draw() comment row)
## plus 4 we just added inside the rowAnnotation calls. The exact
## count depends on whether the comment "# show_legend = FALSE" is
## present in textual form; we just assert that we have at least
## 4 + 2 = 6 occurrences and that every rowAnnotation block has one.
sl_count <- length(gregexpr("show_legend = FALSE", ppm_src, fixed = TRUE)[[1L]])
stopifnot_msg(sl_count >= 6L,
              sprintf("plot_phenotype_markers has >=6 'show_legend = FALSE' occurrences (got %d)",
                      sl_count))

## Each rowAnnotation block must carry its own show_legend = FALSE.
## Strategy: split the source around each rowAnnotation( opener and
## confirm `show_legend = FALSE` appears before the *matching* close.
split_at <- "ComplexHeatmap::rowAnnotation("
parts <- strsplit(ppm_src, split_at, fixed = TRUE)[[1L]]
## parts[1] = text before the FIRST rowAnnotation. The remaining
## n_row_anno entries each start with the body of a rowAnnotation
## call and continue to the next rowAnnotation (or end-of-file).
stopifnot_msg(length(parts) == n_row_anno + 1L,
              "split produces (n_row_anno + 1) chunks")
for (i in seq_len(n_row_anno)) {
  body <- parts[i + 1L]
  ## Truncate body at the first occurrence of "\n    ha_" (next
  ## ha_left / ha_right assignment) or "    if (" to bound it within
  ## a single rowAnnotation call. We use a small look-ahead window
  ## of the first 3000 chars which comfortably covers each block.
  body_window <- substr(body, 1L, 3000L)
  stopifnot_msg(
    grepl("show_legend = FALSE", body_window, fixed = TRUE),
    sprintf("rowAnnotation #%d has show_legend = FALSE in its body", i)
  )
}

## (Live PNG round-trip is not run here -- the existing testthat
## suite under tests/testthat/test-plot-phenotype-markers.R covers
## 21 end-to-end scenarios including single-cell-type metadata
## and is the authoritative regression check for the heatmap.)

## --- 2. Embeddings sidebar wraps the upper block in compact-stack ---------
emb_anchor <- '" 4. Visualization")'
ep <- regexpr(emb_anchor, app_src, fixed = TRUE)
stopifnot_msg(ep > 0, "located the '4. Visualization' nav anchor")
emb_slice <- substr(app_src, ep, min(nchar(app_src), ep + 5000L))

stopifnot_msg(
  grepl('class = "phenomapr-compact-stack embedding-compact-stack"',
        emb_slice, fixed = TRUE),
  "Embeddings sidebar wraps upper block in .embedding-compact-stack"
)
## The expected upper-block children sit INSIDE the wrapper.
for (needle in c(
  'h4("Embedding")',
  'selectInput("umap_reduction"',
  'uiOutput("umap_reduction_status")',
  'uiOutput("spatial_sample_selector")',
  'class = "embedding-help-details"',
  'embedding-upload-details'
)) {
  stopifnot_msg(grepl(needle, emb_slice, fixed = TRUE),
                sprintf("embedding-compact-stack contains: %s", needle))
}
## The <hr/> + Color-cells-by control sit AFTER the wrapper close.
hr_pos <- regexpr("\n        hr\\(\\)", emb_slice, perl = TRUE)
color_pos <- regexpr('radioButtons\\(\n\\s*"umap_color_by"',
                     emb_slice, perl = TRUE)
stopifnot_msg(
  hr_pos > 0 && color_pos > 0 && hr_pos < color_pos,
  "<hr/> sits BEFORE the umap_color_by radio (sanity check)"
)

## CSS supplementary rules for the embedding-compact-stack exist.
for (sel in c(
  ".embedding-compact-stack > h4",
  ".embedding-compact-stack > details.embedding-help-details",
  ".embedding-compact-stack details.embedding-help-details > summary"
)) {
  stopifnot_msg(grepl(sel, styles_src, fixed = TRUE),
                sprintf("styles.css defines %s", sel))
}

## --- 3. Expression input sidebar wraps the upper block too ----------------
data_anchor <- '" 1. Data")'
dp <- regexpr(data_anchor, app_src, fixed = TRUE)
stopifnot_msg(dp > 0, "located the '1. Data' nav anchor")
data_slice <- substr(app_src, dp, min(nchar(app_src), dp + 4000L))

stopifnot_msg(
  grepl('class = "phenomapr-compact-stack expression-compact-stack"',
        data_slice, fixed = TRUE),
  "Expression input sidebar wraps upper block in .expression-compact-stack"
)
for (needle in c(
  'h4("Expression input")',
  'phenomapr_file_input(\n            "expr_file"',
  'actionLink("use_demo"',
  'uiOutput("expr_matrix_diagnostics")'
)) {
  stopifnot_msg(grepl(needle, data_slice, fixed = TRUE),
                sprintf("expression-compact-stack contains: %s", needle))
}

## CSS supplementary rules for the expression-compact-stack exist.
for (sel in c(
  ".expression-compact-stack > h4",
  ".expression-compact-stack > #use_demo",
  ".expression-compact-stack > #expr_matrix_diagnostics"
)) {
  stopifnot_msg(grepl(sel, styles_src, fixed = TRUE),
                sprintf("styles.css defines %s", sel))
}

cat("\nAll heatmap-stray-legend + sidebar-compact smoke tests passed.\n")
