#!/usr/bin/env Rscript
## Smoke tests for the two new "Score by cell type and source" plot
## controls:
##
##   1. radioButtons("score_box_source_palette", ...) -- a Color /
##      Greyscale toggle that swaps the Source-outline palette
##      between the brand colours and a black-to-light-grey ramp.
##
##   2. numericInput("score_box_celltype_legend_ncol", ...) -- the
##      number of columns the cell-type fill legend is laid out in.
##      Defaults to 1; bumping to >= 2 keeps the legend on-canvas
##      when there are many cell types.
##
## Both controls live inside the Visualization-tab SIDEBAR (below
## the embedding download button), keyed off a
## `.phenomapr-score-box-controls` div, and feed input$* into the
## renderPlot() that builds the boxplot. We assert:
##
##   * UI: both inputs are declared, sit inside
##     `.phenomapr-score-box-controls`, and that controls div lives
##     in the sidebar -- BEFORE both `card(`-blocks for the
##     embedding plotOutput and the score-box plotOutput.
##   * Server: the renderPlot consults both inputs when building the
##     plot, the greyscale branch uses grDevices::gray.colors(), the
##     colour branch keeps using pm_brand_palette(), and the cell-
##     type legend honours guide_legend(ncol = ...).
##   * CSS: a companion `.phenomapr-score-box-controls` rule exists.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R <- system.file("shiny", "app.R", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R <- file.path(dev_dir, "app.R")
}
stopifnot(file.exists(app_R))
app_src <- paste(readLines(app_R, warn = FALSE), collapse = "\n")

styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(styles_css))
css_src <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## --- 1. UI: both controls declared --------------------------------------
pal_pos <- regexpr(
  'radioButtons\\(\\s*"score_box_source_palette"',
  app_src, perl = TRUE
)
stopifnot_msg(pal_pos > 0L,
              'radioButtons("score_box_source_palette") declared')

ncol_pos <- regexpr(
  'numericInput\\(\\s*"score_box_celltype_legend_ncol"',
  app_src, perl = TRUE
)
stopifnot_msg(ncol_pos > 0L,
              'numericInput("score_box_celltype_legend_ncol") declared')

## Both choices on the palette toggle are present.
pal_window <- substr(app_src, pal_pos,
                     min(nchar(app_src), pal_pos + 600L))
stopifnot_msg(
  grepl('"Color"\\s*=\\s*"color"', pal_window, perl = TRUE),
  'palette toggle exposes a "Color" choice -> "color"'
)
stopifnot_msg(
  grepl('"Greyscale"\\s*=\\s*"greyscale"', pal_window, perl = TRUE),
  'palette toggle exposes a "Greyscale" choice -> "greyscale"'
)
stopifnot_msg(
  grepl('selected\\s*=\\s*"color"', pal_window, perl = TRUE),
  'palette toggle defaults to "color" (preserves prior plot look)'
)

## --- 2. UI: controls sit in the Visualization-tab SIDEBAR ----------------
##
## There are several `sidebar(width = 360,` calls in app.R (one per
## tab that has a sidebar). We need the one on the Visualization
## tab specifically: anchor on its nav_panel value (`value = "umap"`,
## kept stable for backward-compat with bookmarks/tests).
viz_anchor <- regexpr(
  'value\\s*=\\s*"umap"',
  app_src, perl = TRUE
)
stopifnot_msg(viz_anchor > 0L,
              'found Visualization-tab anchor (value = "umap")')

post_viz <- substr(app_src, viz_anchor, nchar(app_src))
sb_offset <- regexpr(
  'sidebar\\s*=\\s*sidebar\\(\\s*width\\s*=\\s*360',
  post_viz, perl = TRUE
)
stopifnot_msg(sb_offset > 0L,
              "Visualization-tab sidebar(width = 360) declared")
sidebar_pos <- viz_anchor + sb_offset - 1L

## First `card(` after the sidebar opens marks the end of the
## sidebar's UI block (everything between sidebar(...) and the
## first card() lives inside the sidebar).
post_sidebar <- substr(app_src, sidebar_pos, nchar(app_src))
first_card_offset <- regexpr("\\bcard\\s*\\(", post_sidebar, perl = TRUE)
stopifnot_msg(first_card_offset > 0L,
              "first card() after the sidebar found")
first_card_pos <- sidebar_pos + first_card_offset - 1L

stopifnot_msg(
  pal_pos > sidebar_pos && pal_pos < first_card_pos,
  '"score_box_source_palette" radio lives inside the Visualization sidebar'
)
stopifnot_msg(
  ncol_pos > sidebar_pos && ncol_pos < first_card_pos,
  '"score_box_celltype_legend_ncol" numeric lives inside the Visualization sidebar'
)

## They live inside the dedicated `.phenomapr-score-box-controls` div,
## which itself opens inside the sidebar.
controls_div_pos <- regexpr(
  'class\\s*=\\s*"phenomapr-score-box-controls"',
  app_src, perl = TRUE
)
stopifnot_msg(controls_div_pos > 0L,
              'controls wrapped in .phenomapr-score-box-controls div')
stopifnot_msg(
  controls_div_pos > sidebar_pos && controls_div_pos < first_card_pos,
  ".phenomapr-score-box-controls div opens inside the sidebar"
)
stopifnot_msg(
  controls_div_pos < pal_pos && controls_div_pos < ncol_pos,
  "controls div opens BEFORE both inputs"
)

## And -- the boxplot's plotOutput still exists, just downstream
## of the sidebar (inside one of the cards on the same tab).
plotout_pos <- regexpr(
  'plotOutput\\(\\s*"score_box_source_plot"',
  app_src, perl = TRUE
)
stopifnot_msg(plotout_pos > 0L,
              'plotOutput("score_box_source_plot") found')
stopifnot_msg(
  plotout_pos > first_card_pos,
  "plotOutput sits in a card downstream of the sidebar"
)

## --- 3. Server: renderPlot consults both inputs --------------------------
render_pos <- regexpr(
  'output\\$score_box_source_plot\\s*<-\\s*renderPlot',
  app_src, perl = TRUE
)
stopifnot_msg(render_pos > 0L,
              'output$score_box_source_plot renderPlot block found')
render_window <- substr(
  app_src, render_pos,
  min(nchar(app_src), render_pos + 12000L)
)

stopifnot_msg(
  grepl('input\\$score_box_source_palette', render_window, perl = TRUE),
  'renderPlot reads input$score_box_source_palette'
)
stopifnot_msg(
  grepl('input\\$score_box_celltype_legend_ncol', render_window, perl = TRUE),
  'renderPlot reads input$score_box_celltype_legend_ncol'
)

## Greyscale branch uses grDevices::gray.colors(); colour branch keeps
## pm_brand_palette(). Both must be present in the renderPlot block.
stopifnot_msg(
  grepl('grDevices::gray\\.colors', render_window, perl = TRUE),
  'greyscale branch uses grDevices::gray.colors()'
)
stopifnot_msg(
  grepl('pm_brand_palette', render_window, perl = TRUE),
  'colour branch still uses pm_brand_palette()'
)

## Cell-type fill legend is split into N columns via
## guide_legend(ncol = legend_ncol, ...). Match across whitespace /
## newlines because the call is multi-line wrapped.
stopifnot_msg(
  grepl(
    'guide_legend\\([^)]*ncol\\s*=\\s*legend_ncol',
    render_window, perl = TRUE
  ),
  'cell-type guide_legend honours ncol = legend_ncol'
)

## --- 4. CSS companion rule ----------------------------------------------
stopifnot_msg(
  grepl('\\.phenomapr-score-box-controls\\b', css_src, perl = TRUE),
  '.phenomapr-score-box-controls CSS rule exists'
)

cat("\nAll score-box-source plot-control smoke tests passed.\n")
