#!/usr/bin/env Rscript
## Smoke tests for four UI/UX polish changes:
##
##   1. The metadata file picker has been moved INSIDE the
##      `<details id="metadata_upload_details">` collapsible (used to be
##      a sibling rendered statically below the renderUI). The wrapper
##      auto-opens (and gets a `mu-prompt` modifier class) when an
##      expression file is loaded but no metadata was detected. The
##      `<summary>` content is rendered via uiOutput with a
##      `tags$summary(...)` container so HTML5's disclosure semantics
##      stay intact. Open / prompt state is pushed via a new custom
##      message ("phenomapr-set-details-state") handled in helpers.R.
##
##   2. The "PhenoMap() parameters" sub-section in the 3. Score
##      sidebar is wrapped in a `phenomapr-compact-stack
##      score-params-compact-stack` div, with new CSS rules trimming
##      heading / form-group / helpText / button margins.
##
##   3. The Embedding (umap_plot) cell_type and group color branches
##      gain a `guides(color = guide_legend(override.aes =
##      list(size = 4, alpha = 1)))` so their legend point sizes
##      match those of the score_rank_plot.
##
##   4. The "Per-cell-type group enrichment" panel:
##        a. drops the "Only shown when a cell-type column has been
##           selected on the Data tab." helpText.
##        b. orders cell types by ascending median PhenoMapR score
##           (pulled from cell_table()) so the x-axis matches the
##           "Score by cell type and source" boxplot above.
##
## All assertions are made against the source files; we don't spin up
## the live app here.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R      <- system.file("shiny", "app.R",            package = "PhenoMapR")
helpers_R  <- system.file("shiny", "helpers.R",        package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  helpers_R  <- file.path(dev_dir, "helpers.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(helpers_R), file.exists(styles_css))
app_src     <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
helpers_src <- paste(readLines(helpers_R,  warn = FALSE), collapse = "\n")
styles_src  <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## ---- 1. metadata_upload_details: file picker INSIDE <details> -----------
stopifnot_msg(
  grepl('id    = "metadata_upload_details"', app_src, fixed = TRUE),
  "static <details id='metadata_upload_details'> wrapper exists"
)

## Locate the wrapper definition and verify both the summary uiOutput
## AND the file picker are children of the same details() call.
det_open <- regexpr('id    = "metadata_upload_details"', app_src, fixed = TRUE)
stopifnot_msg(det_open > 0L, "located the metadata-upload <details> block")
det_window <- substr(app_src, det_open,
                     min(nchar(app_src), det_open + 2500L))

stopifnot_msg(
  grepl('uiOutput\\(\\s*"metadata_upload_summary"', det_window),
  "summary uiOutput sits inside the <details> block"
)
stopifnot_msg(
  grepl('container = function\\(\\.\\.\\.\\)\\s*\\n?\\s*tags\\$summary',
        det_window),
  "summary uiOutput container is tags$summary"
)
stopifnot_msg(
  grepl('class = "metadata-upload-help"', det_window, fixed = TRUE),
  "helpText sits inside a static .metadata-upload-help div"
)
stopifnot_msg(
  grepl('class = "metadata-upload-file-host"', det_window, fixed = TRUE) &&
    grepl('phenomapr_file_input\\(\\s*\\n?\\s*"meta_file"', det_window),
  "file picker is inside the <details> block (not a sibling)"
)

## The file picker should NOT also exist as a sibling outside of
## metadata-upload-details. We check this by counting the number of
## phenomapr_file_input("meta_file", ...) calls; there should be
## exactly one.
n_meta_pick <- length(gregexpr('phenomapr_file_input\\(\\s*\\n?\\s*"meta_file"',
                               app_src)[[1L]])
stopifnot_msg(
  n_meta_pick == 1L,
  sprintf("phenomapr_file_input('meta_file') appears exactly once (got %d)",
          n_meta_pick)
)

## The renderUI for metadata_upload_summary returns CHILDREN ONLY (a
## tagList or a plain string) -- never a top-level tags$summary().
sum_open <- regexpr("output\\$metadata_upload_summary <- renderUI\\(\\{",
                    app_src)
stopifnot_msg(sum_open > 0L,
              "metadata_upload_summary renderUI block exists")
sum_end  <- regexpr("\\n\\s*\\}\\)\\s*\\n", substr(app_src, sum_open,
                                                  nchar(app_src)))
sum_window <- substr(app_src, sum_open, sum_open + sum_end + 5L)
stopifnot_msg(
  !grepl("tags\\$summary\\(", sum_window),
  "renderUI body emits children only (no nested tags$summary)"
)

## Server pushes both `open` AND `prompt` flags via the new custom
## message.
stopifnot_msg(
  grepl('"phenomapr-set-details-state"', app_src, fixed = TRUE),
  "server-side observe sends phenomapr-set-details-state message"
)
stopifnot_msg(
  grepl('id     = "metadata_upload_details"', app_src, fixed = TRUE) ||
    grepl('id = "metadata_upload_details"',     app_src, fixed = TRUE),
  "message payload references metadata_upload_details element"
)
stopifnot_msg(
  grepl("open\\s*=\\s*prompt", app_src) &&
    grepl("prompt\\s*=\\s*prompt", app_src),
  "message payload carries both open and prompt flags"
)

## helpers.R registers the new handler and toggles `open` + `mu-prompt`.
stopifnot_msg(
  grepl('"phenomapr-set-details-state"', helpers_src, fixed = TRUE),
  "helpers.R defines the phenomapr-set-details-state handler"
)
stopifnot_msg(
  grepl('el\\.classList\\.add\\("mu-prompt"\\)', helpers_src) &&
    grepl('el\\.classList\\.remove\\("mu-prompt"\\)', helpers_src),
  "JS handler toggles the mu-prompt class"
)
stopifnot_msg(
  grepl('el\\.setAttribute\\("open"', helpers_src) &&
    grepl('el\\.removeAttribute\\("open"\\)', helpers_src),
  "JS handler sets/removes the open attribute"
)

## CSS: prompt look is keyed off `.metadata-upload-details.mu-prompt
## > summary`, not the standalone `.metadata-upload-details-prompt`.
stopifnot_msg(
  grepl(".metadata-upload-details.mu-prompt > summary.metadata-upload-details-summary",
        styles_src, fixed = TRUE),
  "styles.css scopes the prompt look to .mu-prompt > summary"
)
stopifnot_msg(
  grepl(".metadata-upload-help", styles_src, fixed = TRUE),
  "styles.css trims margins on .metadata-upload-help"
)

## ---- 2. PhenoMap() parameters compact stack ----------------------------
stopifnot_msg(
  grepl('"phenomapr-compact-stack score-params-compact-stack"',
        app_src, fixed = TRUE),
  "PhenoMap() parameters block wraps in score-params-compact-stack"
)
stopifnot_msg(
  grepl(".score-params-compact-stack > h4", styles_src, fixed = TRUE),
  "styles.css trims margins on .score-params-compact-stack > h4"
)
stopifnot_msg(
  grepl(".score-params-compact-stack .help-block", styles_src, fixed = TRUE),
  "styles.css trims margins on .score-params-compact-stack helpText"
)
stopifnot_msg(
  grepl(".score-params-compact-stack > .btn", styles_src, fixed = TRUE),
  "styles.css trims margins on .score-params-compact-stack buttons"
)

## h4("PhenoMap() parameters") still appears (we didn't lose the
## heading in the wrap).
stopifnot_msg(
  grepl('h4\\("PhenoMap\\(\\) parameters"\\)', app_src),
  "h4('PhenoMap() parameters') heading is still present"
)

## ---- 3. Embedding legend point size override ---------------------------
## Grab the umap_plot renderPlot block and check that BOTH the
## cell_type and group branches now apply the size = 4 / alpha = 1
## guide_legend override (matching the score_rank_plot legend).
umap_open <- regexpr("output\\$umap_plot <- renderPlot\\(\\{", app_src)
stopifnot_msg(umap_open > 0L, "located output$umap_plot renderPlot block")
umap_end <- regexpr("\\n\\s*\\}\\)\\s*\\n\\s*phenomapr_register_plot_download\\(output, \"umap_plot\"",
                    substr(app_src, umap_open, nchar(app_src)))
stopifnot_msg(umap_end > 0L, "located the end of the umap_plot block")
umap_window <- substr(app_src, umap_open, umap_open + umap_end + 200L)

## The override.aes = list(size = 4, alpha = 1) line should appear
## at least 3 times (cell_type, source, group branches).
n_override <- length(gregexpr("override\\.aes = list\\(size = 4, alpha = 1\\)",
                              umap_window)[[1L]])
stopifnot_msg(
  n_override >= 3L,
  sprintf("umap_plot has at least 3 size=4 legend overrides (got %d)",
          n_override)
)

## Make sure the cell_type branch contains the override (the regex is
## anchored to the cell_type branch's labs(color = "Cell type") line).
ct_match <- regexpr('labs\\(color = "Cell type"\\)', umap_window)
stopifnot_msg(ct_match > 0L,
              "located the cell_type branch in umap_plot")
ct_slice <- substr(umap_window, ct_match,
                   min(nchar(umap_window), ct_match + 600L))
stopifnot_msg(
  grepl("override\\.aes = list\\(size = 4, alpha = 1\\)", ct_slice),
  "cell_type branch now has the size=4 legend override"
)

## And the group branch contains it (anchored to the
## scale_color_manual call's `name = "Group"` line).
grp_match <- regexpr('name = "Group"', umap_window, fixed = TRUE)
stopifnot_msg(grp_match > 0L,
              "located the group branch in umap_plot")
grp_slice <- substr(umap_window, grp_match,
                    min(nchar(umap_window), grp_match + 800L))
stopifnot_msg(
  grepl("override\\.aes = list\\(size = 4, alpha = 1\\)", grp_slice),
  "group branch now has the size=4 legend override"
)

## ---- 4. Per-cell-type group enrichment panel: helpText removed +
##         cell-type ordering matches score boxplot ----------------------
## a. helpText removed from the panel UI.
stopifnot_msg(
  !grepl("Only shown when a cell-type column has been selected on the",
         app_src, fixed = TRUE),
  "'Only shown when a cell-type column...' helpText has been removed"
)

## b. group_by_celltype_plot pulls the cell-type ordering from
##    cell_table() and applies it as a factor.
gbc_open <- regexpr("output\\$group_by_celltype_plot <- renderPlot\\(\\{",
                    app_src)
stopifnot_msg(gbc_open > 0L,
              "located output$group_by_celltype_plot block")
gbc_end <- regexpr("\\n\\s*\\}\\)\\s*\\n\\s*phenomapr_register_plot_download\\(output, \"group_by_celltype_plot\"",
                   substr(app_src, gbc_open, nchar(app_src)))
stopifnot_msg(gbc_end > 0L,
              "located end of group_by_celltype_plot block")
gbc_window <- substr(app_src, gbc_open, gbc_open + gbc_end + 200L)

stopifnot_msg(
  grepl("cell_table\\(\\)", gbc_window),
  "group_by_celltype_plot consults cell_table() for ordering"
)
stopifnot_msg(
  grepl("stats::aggregate\\(\\s*\\n?\\s*score ~ cell_type", gbc_window),
  "group_by_celltype_plot computes median score per cell_type"
)
stopifnot_msg(
  grepl("med <- med\\[order\\(med\\$score\\)", gbc_window),
  "ordering is ascending by median score"
)
stopifnot_msg(
  grepl("df_count\\$cell_type <- factor\\(", gbc_window),
  "cell_type is forced to a factor with the computed levels"
)

cat("\nAll metadata-details + sidebar-polish smoke tests passed.\n")
