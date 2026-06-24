#!/usr/bin/env Rscript
## Smoke tests for the redesigned "About markers" panel on the Markers
## tab.
##
## The previous implementation was a static `card_header` + `card_body`
## that always consumed ~30 % of the column height with a wall of
## markdown text. The new layout:
##
##   1. Wraps the panel in <details>/<summary> so users can collapse
##      it (collapsed by default; one click opens it).
##   2. Replaces the long prose with a compact numbered list (1/2/3)
##      that maps onto the new schematic figure.
##   3. Splits the body into a two-column layout: simplified text on
##      the left, the schematic image on the right. The image lives at
##      `inst/figures/PhenoMapR_marker_schematic.png` and is copied
##      into `inst/shiny/www/` so Shiny serves it without any
##      addResourcePath() plumbing.
##
## We assert the source structure and the image asset placement.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (!nzchar(dev_dir) || !dir.exists(dev_dir)) {
  dev_dir <- system.file("shiny", package = "PhenoMapR")
}
stopifnot(dir.exists(dev_dir))
app_R <- file.path(dev_dir, "app.R")
css   <- file.path(dev_dir, "www", "styles.css")
img   <- file.path(dev_dir, "www", "PhenoMapR_marker_schematic.png")
stopifnot(file.exists(app_R), file.exists(css))

app_src <- paste(readLines(app_R, warn = FALSE), collapse = "\n")
css_src <- paste(readLines(css,   warn = FALSE), collapse = "\n")

## --- 1. The image asset is in inst/shiny/www/ ----------------------------
stopifnot_msg(
  file.exists(img),
  "PhenoMapR_marker_schematic.png is shipped under inst/shiny/www/"
)

## --- 2. About-markers panel uses <details>/<summary> ---------------------
stopifnot_msg(
  grepl(
    '(?s)tags\\$details\\(\\s*class\\s*=\\s*"phenomapr-about-markers-details"',
    app_src, perl = TRUE
  ),
  "About-markers panel is wrapped in tags$details with the styled class"
)
stopifnot_msg(
  grepl(
    '(?s)class\\s*=\\s*"phenomapr-about-markers-details".{0,1500}tags\\$summary\\(',
    app_src, perl = TRUE
  ),
  "details block declares a <summary> for the disclosure header"
)
stopifnot_msg(
  grepl('phenomapr-about-markers-summary-label', app_src, fixed = TRUE),
  "summary label uses the styled-summary class"
)
stopifnot_msg(
  grepl('phenomapr-about-markers-summary-hint', app_src, fixed = TRUE) &&
    grepl("click to expand", app_src, fixed = TRUE),
  "summary advertises a \"click to expand\" hint"
)

## --- 3. Two-column layout: text left, schematic right -------------------
stopifnot_msg(
  grepl(
    '(?s)class\\s*=\\s*"phenomapr-about-markers-details".{0,4000}layout_columns\\(',
    app_src, perl = TRUE
  ),
  "details body uses layout_columns() for text + figure split"
)
stopifnot_msg(
  grepl(
    '(?s)layout_columns\\(\\s*col_widths\\s*=\\s*c\\(\\s*7\\s*,\\s*5\\s*\\)',
    app_src, perl = TRUE
  ),
  "layout_columns gives the text column 7 and the figure column 5"
)
stopifnot_msg(
  grepl('class\\s*=\\s*"phenomapr-about-markers-text"',
        app_src, perl = TRUE),
  "left column is the styled text container"
)
stopifnot_msg(
  grepl('class\\s*=\\s*"phenomapr-about-markers-figure"',
        app_src, perl = TRUE),
  "right column is the styled figure container"
)
stopifnot_msg(
  grepl(
    'src\\s*=\\s*"PhenoMapR_marker_schematic\\.png"',
    app_src, perl = TRUE
  ),
  "<img> sources from PhenoMapR_marker_schematic.png in www/"
)

## --- 4. Simplified copy: numbered list mapped to the schematic ----------
for (lbl in c("\\(1\\) Cohort-wide", "\\(2\\)", "\\(3\\)")) {
  stopifnot_msg(
    grepl(lbl, app_src, perl = TRUE),
    sprintf("body text includes the schematic-numbered scope label `%s`", lbl)
  )
}
## And the legacy long-form sentences are gone (regression guard).
stopifnot_msg(
  !grepl("which is exactly\n\\s*the situation where the next mode returns empty",
         app_src, perl = TRUE),
  "long-form 'mid-stringency' prose has been replaced with the schematic-aligned summary"
)

## --- 5. CSS rules for the new layout exist ------------------------------
for (sel in c(
  ".phenomapr-about-markers-details",
  ".phenomapr-about-markers-details > summary",
  ".phenomapr-about-markers-details[open] > summary",
  ".phenomapr-about-markers-figure",
  ".phenomapr-about-markers-img"
)) {
  stopifnot_msg(
    grepl(sel, css_src, fixed = TRUE),
    sprintf("styles.css declares `%s`", sel)
  )
}

cat("\nAll about-markers-collapsible smoke tests passed.\n")
