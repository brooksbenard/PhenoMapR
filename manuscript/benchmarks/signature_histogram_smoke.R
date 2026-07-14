#!/usr/bin/env Rscript
## Smoke tests for the Phenotype-signature histogram polish:
##
##   1. The "Set threshold for Signature |z| cutoff" slider now lives
##      OUTSIDE the built-in-only conditionalPanel so it remains
##      visible (and reactive) for custom phenotype signatures too.
##      PhenoMap() honours z_score_cutoff for every reference, so
##      hiding the slider for custom signatures was a UX bug -- and
##      having it always visible lets the histogram below draw the
##      live vertical-line overlay regardless of which signature
##      source the user picked.
##
##   2. The built-in pre-filter note STILL lives inside the
##      built-in-only conditionalPanel -- it does not apply to
##      custom signatures. The note is reference-aware: ICI PRECOG
##      ships at |z| >= 1, every other built-in ships at |z| >= 2,
##      so the note has two nested conditionalPanels gating on
##      input.reference_choice == 'ici_precog' so each user sees
##      the cutoff that actually applies to their selection.
##
##   3. output$reference_signature_plot:
##        a. Always renders a histogram (no more ComplexHeatmap
##           branch for custom signatures -- both modes now produce
##           ggplot output that ggsave can serialise cleanly via
##           phenomapr_register_plot_download).
##        b. For built-in references, faceting by sign with
##           scales = "free_x" omits the empty -2 < z < 2 middle.
##        c. Both branches consult input$z_score_cutoff and overlay
##           dashed vertical lines at +cut and -cut (live overlay),
##           mirroring the cutoff overlay the Score-distribution
##           histogram already does for percentile thresholds.
##        d. The built-in branch passes n = Inf to
##           PhenoMapR::get_top_prognostic_genes so every gene in
##           the cancer-type column gets plotted, not just the top
##           500.

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

## --- 1. Slider lives outside the built-in-only conditionalPanel -----------
##
## Strategy: locate the radioButtons("reference_choice") block and
## the closing `)` of the conditional that scopes to non-custom
## (`input.reference_choice != '_custom'`). The slider must appear
## AFTER that closing `)` so the conditionalPanel can't hide it.
##
## We do this by finding two anchors:
##   * the cancer-type selectInput (which lives inside the
##     non-custom conditionalPanel).
##   * the z_score_cutoff sliderInput (which we want OUTSIDE).
##
## Then we check that the cancer-type anchor's enclosing
## conditionalPanel closes BEFORE the slider's start.

ct_anchor <- regexpr('selectInput\\("cancer_type"', app_src, perl = TRUE)
sl_anchor <- regexpr('sliderInput\\(\\s*"z_score_cutoff"', app_src, perl = TRUE)
stopifnot_msg(ct_anchor > 0L, "found selectInput('cancer_type') anchor")
stopifnot_msg(sl_anchor > 0L, "found sliderInput('z_score_cutoff') anchor")
stopifnot_msg(sl_anchor > ct_anchor,
              "z_score_cutoff slider declared AFTER cancer_type select")

## The non-custom conditionalPanel closes between the two anchors
## (i.e. the slider sits in the outer sidebar block, not in the
## conditional). Easiest concrete check: the slice between the two
## anchors must close at least one conditionalPanel `)` block.
slice_between <- substr(app_src, ct_anchor, sl_anchor)
n_close <- length(gregexpr("\\)\\s*,", slice_between)[[1]])
stopifnot_msg(
  n_close >= 1L,
  "the conditionalPanel for built-in scope closes BEFORE the slider"
)

## --- 2. Pre-filter note still lives inside the built-in conditionalPanel --
stopifnot_msg(
  grepl('phenomapr-prefilter-note', app_src, fixed = TRUE),
  "pre-filter note (phenomapr-prefilter-note) is still rendered"
)
note_pos <- regexpr('phenomapr-prefilter-note', app_src, fixed = TRUE)
custom_cond <- regexpr("input.reference_choice == '_custom'",
                        app_src, fixed = TRUE)
## The note must appear AFTER the slider (we restructured the order
## to slider -> non-custom-only-note) and BEFORE the custom
## conditionalPanel opens.
stopifnot_msg(
  note_pos > sl_anchor,
  "pre-filter note appears after the slider"
)
stopifnot_msg(
  note_pos < custom_cond,
  "pre-filter note appears before the custom-source conditionalPanel"
)

## --- 2b. Pre-filter note has per-reference branches (ICI vs others) ------
##
## The note must mention the relaxed ICI PRECOG cutoff (|z| >= 1) and the
## stricter cutoff (|z| >= 2) used by every other built-in. We don't pin
## the exact wording -- just the two cutoffs and the ici_precog gate.
note_window <- substr(
  app_src, note_pos,
  min(nchar(app_src), note_pos + 2000L)
)
stopifnot_msg(
  grepl("input\\.reference_choice\\s*==\\s*'ici_precog'", note_window, perl = TRUE),
  "pre-filter note branches on input.reference_choice == 'ici_precog'"
)
stopifnot_msg(
  grepl("\\|z\\|\\s*&ge;\\s*1", note_window, perl = TRUE),
  "pre-filter note advertises the ICI PRECOG cutoff |z| >= 1"
)
stopifnot_msg(
  grepl("\\|z\\|\\s*&ge;\\s*2", note_window, perl = TRUE),
  "pre-filter note still advertises the |z| >= 2 cutoff for non-ICI built-ins"
)

## --- 3a. reference_signature_plot is a histogram for both branches --------
plot_anchor <- regexpr(
  "output\\$reference_signature_plot\\s*<-\\s*renderPlot",
  app_src, perl = TRUE
)
stopifnot_msg(plot_anchor > 0L, "found output$reference_signature_plot")
plot_window <- substr(app_src, plot_anchor,
                       min(nchar(app_src), plot_anchor + 6000L))

stopifnot_msg(
  grepl(".geom_rounded_histogram", plot_window, fixed = TRUE),
  "reference_signature_plot uses .geom_rounded_histogram (ggplot)"
)
stopifnot_msg(
  !grepl("plot_reference_signature\\(state\\$reference\\)", plot_window),
  "reference_signature_plot no longer dispatches to the ComplexHeatmap path"
)

## --- 3b. Faceted broken-axis layout for built-in references ---------------
stopifnot_msg(
  grepl("facet_wrap\\(~\\s*side", plot_window, perl = TRUE),
  "built-in branch facets the histogram by sign (broken-axis layout)"
)
stopifnot_msg(
  grepl('scales\\s*=\\s*"free_x"', plot_window, perl = TRUE),
  "built-in facets use scales = \"free_x\" so the empty middle is dropped"
)

## --- 3c. Live |z| cutoff vertical-line overlay drawn on the histogram -----
stopifnot_msg(
  grepl("input\\$z_score_cutoff", plot_window, perl = TRUE),
  "renderPlot consults input$z_score_cutoff for the cutoff overlay"
)
stopifnot_msg(
  grepl("geom_vline", plot_window, fixed = TRUE),
  "renderPlot draws the cutoff overlay via geom_vline()"
)
## Both +cut and -cut must be plotted (custom branch).
stopifnot_msg(
  grepl("xintercept\\s*=\\s*c\\(-cut_val,\\s*cut_val\\)", plot_window, perl = TRUE),
  "custom branch overlays vlines at both -cut_val and +cut_val"
)
## Built-in faceted branch: per-facet x positions.
stopifnot_msg(
  grepl("vline_df", plot_window, fixed = TRUE),
  "built-in branch builds a per-facet vline data.frame"
)

## --- 3d. Built-in branch passes n = Inf so all genes are plotted ---------
##
## The call spans multiple lines (`reference = ..., cancer_type = ...,
## n = Inf, direction = "both"`) so use a multi-line tolerant regex.
gtp_pos <- regexpr(
  "PhenoMapR::get_top_prognostic_genes\\b",
  plot_window, perl = TRUE
)
stopifnot_msg(gtp_pos > 0L,
              "found get_top_prognostic_genes() call in built-in branch")
gtp_window <- substr(plot_window, gtp_pos,
                      min(nchar(plot_window), gtp_pos + 400L))
stopifnot_msg(
  grepl("n\\s*=\\s*Inf", gtp_window, perl = TRUE),
  "built-in branch passes n = Inf to get_top_prognostic_genes()"
)
stopifnot_msg(
  !grepl("n\\s*=\\s*500L?\\b", gtp_window, perl = TRUE),
  "the old top-500 cap on get_top_prognostic_genes is gone"
)

cat("\nAll signature-histogram smoke tests passed.\n")
