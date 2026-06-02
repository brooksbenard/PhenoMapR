#!/usr/bin/env Rscript
## Smoke tests for two follow-up fixes:
##
##   A. Aggressive popup-hide fallback. The previous proactive-hide
##      relied on Shiny dispatching `shiny:value` / `shiny:visualchange`
##      events to the document level, but in some Shiny version /
##      browser combinations those events do NOT bubble to document
##      reliably. Adding more event listeners (`shiny:flushed`,
##      `shiny:message`) helps, but the bullet-proof fix is a
##      requestAnimationFrame poll that watches the body's
##      `.shiny-busy` class -- which Shiny.js maintains itself
##      regardless of event-dispatch quirks. When the popup has been
##      visible for >= SHOWN_MIN_MS AND the body has been non-busy
##      for >= IDLE_GRACE_MS, the poll hides the popup. Plus an
##      ABSOLUTE_TIMEOUT (3 min) hard cap as a final safety net.
##
##   B. Score-by-cell-type-and-source plot star labels. Previously
##      drawn via ggsignif::geom_signif(manual = TRUE), which in some
##      ggsignif/ggplot2 combos drops short labels ("*", "**") while
##      keeping long ones ("***"). We now hand-draw each bracket from
##      pval_df via geom_segment + geom_text so every bracket
##      reliably gets its label, AND we expand the y axis so labels
##      are never clipped at the panel edge. A plot subtitle now
##      explains the significance levels.

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

## ---- A. rAF body-class polling fallback ---------------------------------
stopifnot_msg(
  grepl("var IDLE_GRACE_MS", helpers_src, fixed = TRUE),
  "helpers.R declares IDLE_GRACE_MS"
)
stopifnot_msg(
  grepl("var ABSOLUTE_TIMEOUT", helpers_src, fixed = TRUE),
  "helpers.R declares ABSOLUTE_TIMEOUT"
)
stopifnot_msg(
  grepl("var bodyIdleSince", helpers_src, fixed = TRUE),
  "helpers.R tracks bodyIdleSince"
)
stopifnot_msg(
  grepl("function isShinyBusy", helpers_src, fixed = TRUE),
  "helpers.R defines isShinyBusy()"
)
stopifnot_msg(
  grepl('classList\\.contains\\("shiny-busy"\\)', helpers_src),
  "isShinyBusy reads the body's .shiny-busy class"
)
stopifnot_msg(
  grepl("function busyPollTick", helpers_src, fixed = TRUE),
  "helpers.R defines busyPollTick()"
)

## busyPollTick body must trigger renderHide once shiny has been
## idle for >= IDLE_GRACE_MS and the popup itself has been visible
## for >= SHOWN_MIN_MS.
poll_anchor <- "function busyPollTick"
pp <- regexpr(poll_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(pp > 0L, "located busyPollTick body")
poll_window <- substr(helpers_src, pp, min(nchar(helpers_src), pp + 1500L))
stopifnot_msg(
  grepl("idleFor >= IDLE_GRACE_MS", poll_window, fixed = TRUE),
  "busyPollTick gates the hide on idleFor >= IDLE_GRACE_MS"
)
stopifnot_msg(
  grepl("visibleFor >= SHOWN_MIN_MS", poll_window, fixed = TRUE),
  "busyPollTick gates the hide on visibleFor >= SHOWN_MIN_MS"
)
stopifnot_msg(
  grepl("visibleFor >= ABSOLUTE_TIMEOUT", poll_window, fixed = TRUE),
  "busyPollTick has an ABSOLUTE_TIMEOUT escape hatch"
)
stopifnot_msg(
  grepl("requestAnimationFrame(busyPollTick)", poll_window, fixed = TRUE),
  "busyPollTick reschedules via requestAnimationFrame"
)

## bind() must kick off the rAF chain.
bind_anchor <- "function bind() {"
ba <- regexpr(bind_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(ba > 0L, "located bind() body")
bind_window <- substr(helpers_src, ba, min(nchar(helpers_src), ba + 3000L))
stopifnot_msg(
  grepl("requestAnimationFrame(busyPollTick)", bind_window, fixed = TRUE),
  "bind() starts the busyPollTick rAF chain"
)

## Additional event listeners are wired up.
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:flushed",\\s*onShinyFlushed',
        helpers_src),
  "shiny:flushed wired via native addEventListener"
)
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:message",\\s*onShinyMessage',
        helpers_src),
  "shiny:message wired via native addEventListener"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:flushed",\\s*onShinyFlushed\\)',
        helpers_src),
  "shiny:flushed wired via jQuery .on()"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:message",\\s*onShinyMessage\\)',
        helpers_src),
  "shiny:message wired via jQuery .on()"
)

## ---- B. Manual bracket / star rendering on score boxplot ----------------

## Locate the score_box_source_plot block.
sbsp_open <- regexpr("output\\$score_box_source_plot <- renderPlot\\(\\{",
                     app_src)
stopifnot_msg(sbsp_open > 0L,
              "located output$score_box_source_plot block")
sbsp_end <- regexpr("\\n\\s*\\}\\)\\s*\\n\\s*phenomapr_register_plot_download\\(output, \"score_box_source_plot\"",
                    substr(app_src, sbsp_open, nchar(app_src)))
stopifnot_msg(sbsp_end > 0L,
              "located the end of the score_box_source_plot block")
sbsp_window <- substr(app_src, sbsp_open, sbsp_open + sbsp_end + 200L)

## ggsignif::geom_signif() must NOT be used in the bracket path.
## We allow the substring to appear in a comment (the rationale block)
## but the actual `p <- p + ggsignif::geom_signif(...)` call must be
## gone. We strip out any line whose first non-space character is `#`
## before scanning.
sbsp_no_comments <- paste(
  Filter(function(L) !grepl("^\\s*#", L),
         strsplit(sbsp_window, "\n", fixed = TRUE)[[1L]]),
  collapse = "\n"
)
stopifnot_msg(
  !grepl("ggsignif::geom_signif", sbsp_no_comments, fixed = TRUE),
  "score_box_source_plot no longer calls ggsignif::geom_signif() (only mentioned in comments)"
)

## Manual brackets: geom_segment for spine + tips, geom_text for label.
stopifnot_msg(
  grepl("geom_segment\\(\\s*\\n?\\s*data = seg_df", sbsp_window) ||
    grepl("data = seg_df", sbsp_window, fixed = TRUE),
  "score_box_source_plot draws the bracket spine via geom_segment"
)
stopifnot_msg(
  grepl("data = tick_df", sbsp_window, fixed = TRUE),
  "score_box_source_plot draws bracket tips via geom_segment"
)
stopifnot_msg(
  grepl("data = text_df", sbsp_window, fixed = TRUE),
  "score_box_source_plot draws the bracket label via geom_text"
)
stopifnot_msg(
  grepl('label = as\\.character\\(pval_df\\$label\\)', sbsp_window),
  "label data uses as.character(pval_df$label) so '*' / '**' / '***' all stringify"
)

## Manual draw uses fontface = 'bold' and size 5 so star labels read
## clearly (the ggsignif default of textsize=4 looked anemic).
stopifnot_msg(
  grepl('fontface = "bold"', sbsp_window, fixed = TRUE),
  "geom_text label uses fontface = 'bold'"
)

## Y axis must be expanded so the topmost star labels never clip.
stopifnot_msg(
  grepl("ggplot2::expand_limits\\(y\\s*=\\s*y_max\\)", sbsp_window),
  "score_box_source_plot expands y limits to fit star labels"
)

## Subtitle explaining significance levels.
stopifnot_msg(
  grepl('"Significance: \\*\\*\\* p < 0.001, \\*\\* p < 0.01, \\* p < 0.05"',
        sbsp_window),
  "score_box_source_plot now carries a significance-key subtitle"
)
stopifnot_msg(
  grepl("plot.subtitle\\s*=\\s*element_text", sbsp_window),
  "score_box_source_plot styles plot.subtitle"
)

## render_kind == 'annotation' branch is preserved (the ANOVA-only case
## is still drawn as a single centred label, no stars).
stopifnot_msg(
  grepl('identical\\(render_kind, "annotation"\\)', sbsp_window),
  "annotation render_kind branch is still routed separately"
)

cat("\nAll popup-poll + significance-stars smoke tests passed.\n")
