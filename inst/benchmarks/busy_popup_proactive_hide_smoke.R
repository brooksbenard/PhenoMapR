#!/usr/bin/env Rscript
## Smoke tests for the popup-lingering fix.
##
## When R bundles many invalidatedOutputValues into a single websocket
## frame after a long heavy op, the centered busy popup could remain on
## screen for several seconds AFTER the new content was already painted
## behind it -- the explicit phenomapr_busy_hide() custom message was
## queued in the same frame as the heavy output payloads and the
## browser sometimes processed the output values before the hide
## message in that batch.
##
## The fix has two prongs:
##
##   A. JS-side proactive hide: helpers.R registers a shiny:value /
##      shiny:visualchange listener that hides the popup as soon as
##      ANY output update lands at the browser, provided the popup has
##      been visible for at least SHOWN_MIN_MS (a small grace window
##      so fast ops whose show+hide land in the same flush are not
##      pre-empted before they ever paint).
##
##   B. R-side deferred state assignment: the heavy state-mutating
##      observers (run_score, run_markers, diag_run_clean) now call
##      phenomapr_busy_hide() BEFORE the state assignment, then defer
##      the state mutation via later::later(..., delay = 0) so the
##      hide message lands in a flush ahead of the cascading reactive
##      output updates.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R     <- system.file("shiny", "app.R",     package = "PhenoMapR")
helpers_R <- system.file("shiny", "helpers.R", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R     <- file.path(dev_dir, "app.R")
  helpers_R <- file.path(dev_dir, "helpers.R")
}
stopifnot(file.exists(app_R), file.exists(helpers_R))
app_src     <- paste(readLines(app_R,     warn = FALSE), collapse = "\n")
helpers_src <- paste(readLines(helpers_R, warn = FALSE), collapse = "\n")

## ---- A. JS-side proactive hide --------------------------------------------
stopifnot_msg(
  grepl("var shownAt", helpers_src, fixed = TRUE),
  "helpers.R declares the shownAt timestamp variable"
)
stopifnot_msg(
  grepl("var SHOWN_MIN_MS", helpers_src, fixed = TRUE),
  "helpers.R declares the SHOWN_MIN_MS grace-window constant"
)
stopifnot_msg(
  grepl("function maybeProactiveHide", helpers_src, fixed = TRUE),
  "helpers.R defines maybeProactiveHide()"
)

## maybeProactiveHide() must guard with the SHOWN_MIN_MS grace window
## (so a fast op's show+hide doesn't get pre-empted before it paints)
## and use setTimeout(0) so the hide runs after the current handler.
proactive_anchor <- "function maybeProactiveHide"
pa <- regexpr(proactive_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(pa > 0L, "located maybeProactiveHide() body")
pa_window <- substr(helpers_src, pa, min(nchar(helpers_src), pa + 800L))
stopifnot_msg(
  grepl("SHOWN_MIN_MS", pa_window, fixed = TRUE),
  "maybeProactiveHide() respects the SHOWN_MIN_MS grace window"
)
stopifnot_msg(
  grepl("setTimeout(function", pa_window, fixed = TRUE) ||
    grepl("setTimeout\\(function", pa_window),
  "maybeProactiveHide() defers via setTimeout(..., 0)"
)
stopifnot_msg(
  grepl("renderHide()", pa_window, fixed = TRUE),
  "maybeProactiveHide() ultimately calls renderHide()"
)

## Both onShinyValue and onShinyVisualChange hooks must call
## maybeProactiveHide.
stopifnot_msg(
  grepl('function onShinyValue\\(\\) \\{ maybeProactiveHide\\("shiny:value"\\)',
        helpers_src),
  "onShinyValue forwards to maybeProactiveHide"
)
stopifnot_msg(
  grepl('function onShinyVisualChange\\(\\) \\{ maybeProactiveHide\\("shiny:visualchange"\\)',
        helpers_src),
  "onShinyVisualChange forwards to maybeProactiveHide"
)

## Both event listeners must be wired up via BOTH the native and
## jQuery layers so we are robust against either event-dispatch
## strategy Shiny uses on a given platform.
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:value",\\s*onShinyValue',
        helpers_src),
  "shiny:value is wired via native addEventListener"
)
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:visualchange",\\s*onShinyVisualChange',
        helpers_src),
  "shiny:visualchange is wired via native addEventListener"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:value",\\s*onShinyValue\\)',
        helpers_src),
  "shiny:value is wired via jQuery (.on)"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:visualchange",\\s*onShinyVisualChange\\)',
        helpers_src),
  "shiny:visualchange is wired via jQuery (.on)"
)

## renderShow() records shownAt; renderHide() resets it. Without this
## bookkeeping the grace-window guard cannot work.
rs_anchor <- "function renderShow"
rs <- regexpr(rs_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(rs > 0L, "located renderShow() body")
rs_window <- substr(helpers_src, rs, min(nchar(helpers_src), rs + 600L))
stopifnot_msg(
  grepl("shownAt = ", rs_window, fixed = TRUE),
  "renderShow() records the shownAt timestamp"
)
rh_anchor <- "function renderHide"
rh <- regexpr(rh_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(rh > 0L, "located renderHide() body")
rh_window <- substr(helpers_src, rh, min(nchar(helpers_src), rh + 400L))
stopifnot_msg(
  grepl("shownAt = 0", rh_window, fixed = TRUE),
  "renderHide() resets the shownAt timestamp"
)

## ---- B. R-side deferred state assignments ---------------------------------
## Helper: pull a window of N chars starting at a fixed anchor.
window_after <- function(src, anchor, n = 1500L) {
  pos <- regexpr(anchor, src, fixed = TRUE)
  if (pos < 1L) return(NULL)
  substr(src, pos, min(nchar(src), pos + n))
}

## Run-score observer: should call hide BEFORE state$scores assignment
## and the assignment must live inside a later::later() call.
score_window <- window_after(
  app_src,
  'phenomapr_busy_show("Computing PhenoMap scores...", ref_label)',
  3000L
)
stopifnot_msg(!is.null(score_window),
              "located the run_score observer body")
## hide call sits between the heavy work and the state assignment.
hide_pos <- regexpr("phenomapr_busy_hide()", score_window, fixed = TRUE)
state_pos <- regexpr("state$scores <- scores", score_window, fixed = TRUE)
later_pos <- regexpr("later::later(function()", score_window, fixed = TRUE)
stopifnot_msg(
  hide_pos > 0L && state_pos > 0L && later_pos > 0L,
  "run_score has hide, state$scores assignment, and later::later wrapper"
)
stopifnot_msg(
  hide_pos < later_pos && later_pos < state_pos,
  "run_score: hide() comes before later::later() which wraps state$scores"
)
## And there is no on.exit(phenomapr_busy_hide()) inside the
## observer (we replaced it with the inline hide + deferred state).
stopifnot_msg(
  !grepl("on\\.exit\\(phenomapr_busy_hide\\(\\), add = TRUE\\)",
         score_window) ||
    regexpr("on\\.exit\\(phenomapr_busy_hide\\(\\), add = TRUE\\)",
            score_window) > later_pos,
  "run_score does not rely on on.exit for the popup hide"
)

## Run-markers observer: same pattern.
markers_window <- window_after(
  app_src,
  'phenomapr_busy_show(\n      "Finding marker genes..."',
  3000L
)
stopifnot_msg(!is.null(markers_window),
              "located the run_markers observer body")
hide_pos2  <- regexpr("phenomapr_busy_hide()", markers_window, fixed = TRUE)
state_pos2 <- regexpr("state$markers <- markers",
                      markers_window, fixed = TRUE)
later_pos2 <- regexpr("later::later(function()", markers_window, fixed = TRUE)
stopifnot_msg(
  hide_pos2 > 0L && state_pos2 > 0L && later_pos2 > 0L,
  "run_markers has hide, state$markers assignment, and later::later wrapper"
)
stopifnot_msg(
  hide_pos2 < later_pos2 && later_pos2 < state_pos2,
  "run_markers: hide() comes before later::later() wrapping state$markers"
)

## Matrix-cleanup observer: same pattern.
clean_window <- window_after(
  app_src,
  'phenomapr_busy_show("Cleaning matrix..."',
  3000L
)
stopifnot_msg(!is.null(clean_window),
              "located the diag_run_clean observer body")
hide_pos3  <- regexpr("phenomapr_busy_hide()", clean_window, fixed = TRUE)
state_pos3 <- regexpr("state$expression <- res$matrix",
                      clean_window, fixed = TRUE)
later_pos3 <- regexpr("later::later(function()", clean_window, fixed = TRUE)
stopifnot_msg(
  hide_pos3 > 0L && state_pos3 > 0L && later_pos3 > 0L,
  "diag_run_clean has hide, state$expression assignment, and later wrapper"
)
stopifnot_msg(
  hide_pos3 < later_pos3 && later_pos3 < state_pos3,
  "diag_run_clean: hide() comes before later::later() wrapping the state assignment"
)

cat("\nAll busy-popup proactive-hide smoke tests passed.\n")
