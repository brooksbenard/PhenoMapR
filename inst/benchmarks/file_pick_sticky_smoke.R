#!/usr/bin/env Rscript
## Smoke tests for the "sticky file pick" behaviour added to
## phenomapr_file_pick() in inst/shiny/helpers.R.
##
## Background: prior to this change, every time the user re-opened
## the shinyFiles browse dialog (or shinyFiles emitted an empty
## selection for any other reason), the picker observer reset
## pick_state(NULL), silently unloading the currently-loaded file.
## That's a footgun -- a user with a 2GB anndata file open would
## click "Browse..." to verify the path and lose all their loaded
## metadata + scores + groups + markers in the process.
##
## The fix is to drop pick_state ONLY in response to:
##   * the explicit "x" remove button next to the chosen filename
##     (input[[<picker_id>_clear]]), OR
##   * a successful new pick whose parseFilePaths() yields a real
##     file row (which naturally replaces the previous pick via the
##     success branch).
##
## Empty / NULL / integer(0) selections from the picker -- which
## happen when the user opens and cancels the dialog, or when the
## input is re-armed mid-session -- must be no-ops.
##
## This test pattern-matches helpers.R for both the comment hint
## (so future maintainers know the constraint exists) and the
## structural change in the observer body.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

helpers_R <- system.file("shiny", "helpers.R", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  helpers_R <- file.path(dev_dir, "helpers.R")
}
stopifnot(file.exists(helpers_R))
helpers_src <- paste(readLines(helpers_R, warn = FALSE), collapse = "\n")

## --- 1. The phenomapr_file_pick() observer carries the sticky-pick comment.
##         (Doubles as documentation for future maintainers.)
stopifnot_msg(
  grepl("sticky picks", helpers_src, fixed = TRUE),
  "phenomapr_file_pick() declares the 'sticky picks' contract"
)

## --- 2. The empty-selection branch must be a NO-OP, not a clear.
##         In R syntax: `if (is.null(sel) || identical(sel, integer(0))) {`
##         followed by an early `return()` -- and *no* `pick_state(NULL)`
##         on that branch.
##
## We grep for the early return and the absence of a clearing call
## on the empty-selection branch.
empty_sel_pattern <- paste0(
  "if \\(is\\.null\\(sel\\) \\|\\| identical\\(sel, integer\\(0\\)\\)\\) \\{",
  "[^}]*return\\(\\)"
)
m <- regmatches(helpers_src,
                regexpr(empty_sel_pattern, helpers_src, perl = TRUE))
stopifnot_msg(
  length(m) == 1L && nzchar(m),
  "empty picker selection branch returns early"
)
stopifnot_msg(
  !grepl("pick_state\\(NULL\\)", m),
  "empty picker selection branch does NOT clear pick_state"
)

## --- 3. The parse-failure branch is also a no-op (don't drop the
##         previously-loaded file just because shinyFiles handed us
##         a selection we couldn't materialise into a path).
parse_fail_pattern <- paste0(
  "if \\(is\\.null\\(path_df\\) \\|\\| nrow\\(path_df\\) == 0L\\) \\{",
  "[^}]*return\\(\\)"
)
m2 <- regmatches(helpers_src,
                 regexpr(parse_fail_pattern, helpers_src, perl = TRUE))
stopifnot_msg(
  length(m2) == 1L && nzchar(m2),
  "parse-failure branch in phenomapr_file_pick() returns early"
)
stopifnot_msg(
  !grepl("pick_state\\(NULL\\)", m2),
  "parse-failure branch does NOT clear pick_state"
)

## --- 4. The remove-button branch (input[[clear_id]]) IS still a
##         pick_state(NULL) -- this is the only user gesture that drops
##         a loaded file (apart from selecting a new one, which is
##         handled by the success branch overwriting pick_state).
##
## Locate the clear_id observer and confirm it sets pick_state(NULL).
clear_obs <- regmatches(
  helpers_src,
  regexpr(
    "shiny::observeEvent\\(input\\[\\[clear_id\\]\\][^}]*pick_state\\(NULL\\)",
    helpers_src,
    perl = TRUE
  )
)
stopifnot_msg(
  length(clear_obs) == 1L && nzchar(clear_obs),
  "explicit 'x' remove button still clears pick_state to NULL"
)

cat("\nAll file-pick sticky smoke tests passed.\n")
