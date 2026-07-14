#!/usr/bin/env Rscript
## Smoke tests for two polish changes on the 3. Score tab:
##
##   1. The "PhenoMapR computes a weighted sum of expression x reference
##      z-score across signature genes..." paragraph that used to live
##      directly under the slot/assay conditionalPanel header is now
##      collapsed inside a `tags$details(class = "score-slot-details",
##      tags$summary("Details"), helpText(...))` chip.
##   2. The "About scoring", "About phenotype tails", and "Group sizes"
##      cards on the 3. Score tab now use `fill = FALSE` so they size
##      to their (short) prose / table bodies instead of stretching to
##      whatever the taller sibling card pushes the row to. The
##      previous default behaviour produced a wide blank band beneath
##      each of these cards.
##
## Both changes are static structural edits in inst/shiny/app.R +
## inst/shiny/www/styles.css; we assert against the source.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R     <- system.file("shiny", "app.R",            package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir   <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(styles_css))
app_src    <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
styles_src <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## --- 1. <details> wrapper around the slot rationale -----------------------
stopifnot_msg(
  grepl('class = "score-slot-details"', app_src, fixed = TRUE),
  "Score sidebar uses <details class='score-slot-details'> wrapper"
)
## The summary label is "Details".
det_anchor <- 'class = "score-slot-details"'
dp <- regexpr(det_anchor, app_src, fixed = TRUE)
stopifnot_msg(dp > 0, "located the score-slot-details wrapper")
det_slice <- substr(app_src, dp, min(nchar(app_src), dp + 1500L))
stopifnot_msg(
  grepl('tags\\$summary\\("Details"\\)', det_slice),
  "score-slot-details summary is labeled 'Details'"
)
## The original helpText paragraph is INSIDE the <details> body.
stopifnot_msg(
  grepl("PhenoMapR computes a weighted sum of expression", det_slice, fixed = TRUE),
  "the explanatory paragraph lives INSIDE the <details> wrapper"
)
## And it's NOT also rendered as a top-level helpText elsewhere in the
## sidebar (the move should be exclusive, not a duplicate).
n_para <- length(gregexpr(
  "PhenoMapR computes a weighted sum of expression",
  app_src, fixed = TRUE
)[[1L]])
stopifnot_msg(
  n_para == 1L,
  sprintf("paragraph appears exactly ONCE in app.R (got %d)", n_para)
)

## CSS supplementary rule for the new wrapper.
stopifnot_msg(
  grepl(".score-slot-details", styles_src, fixed = TRUE),
  "styles.css defines .score-slot-details"
)
stopifnot_msg(
  grepl(".score-slot-details > summary", styles_src, fixed = TRUE),
  "styles.css styles the .score-slot-details summary"
)

## --- 2. fill = FALSE on the three short-body cards ------------------------
## Locate each card by its header label and confirm the card() call
## that wraps it carries `fill = FALSE` in its arg list.

## Helper: extract a window of N chars starting at the position of an
## anchor we know lives in a card_header() child. We then walk
## backwards a small fixed distance to capture the `card(` opener that
## introduces the card.
card_window <- function(src, header_anchor) {
  pos <- regexpr(header_anchor, src, fixed = TRUE)
  if (pos < 1L) return(NULL)
  start <- max(1L, pos - 600L)
  end   <- min(nchar(src), pos + 50L)
  substr(src, start, end)
}

for (anchor in c(
  "card_header(icon(\"circle-info\"), \" About scoring\")",
  "card_header(icon(\"circle-info\"), \" About phenotype tails\")",
  "card_header(tags$strong(\"Group sizes\"))"
)) {
  win <- card_window(app_src, anchor)
  stopifnot_msg(!is.null(win),
                sprintf("located card with header %s", sQuote(anchor)))
  stopifnot_msg(
    grepl("fill = FALSE", win, fixed = TRUE),
    sprintf("card with header %s now sets fill = FALSE",
            sQuote(anchor))
  )
}

cat("\nAll Score-tab sidebar polish smoke tests passed.\n")
