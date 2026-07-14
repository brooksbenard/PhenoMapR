#!/usr/bin/env Rscript
## Smoke tests for two follow-up fixes:
##
##   A. Custom-message-driven popup with debounce + watchdog. The
##      centered "busy" popup is now driven by R-side custom messages
##      (phenomapr-busy-show / phenomapr-busy-hide) with a JS-side
##      SHOW_DELAY_MS debounce on the show side. R observers MUST
##      defer their heavy synchronous work into later::later() so the
##      busy_show message reaches the browser BEFORE libuv is blocked
##      by the work. A 3-minute absolute watchdog is still wired as a
##      last-resort hard cap; shiny:idle is an advisory fallback hide.
##      (See busy_popup_proactive_hide_smoke.R for the full
##      structural-shape assertions.)
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

## ---- A. Custom-message-driven popup with debounce + watchdog -----------
stopifnot_msg(
  grepl("var SHOW_DELAY_MS", helpers_src, fixed = TRUE),
  "helpers.R declares SHOW_DELAY_MS (debounce before showing)"
)
stopifnot_msg(
  grepl("var IDLE_SETTLE_MS", helpers_src, fixed = TRUE),
  "helpers.R declares IDLE_SETTLE_MS (advisory shiny:idle hide delay)"
)
stopifnot_msg(
  grepl("var ABSOLUTE_TIMEOUT", helpers_src, fixed = TRUE),
  "helpers.R declares ABSOLUTE_TIMEOUT (3-min watchdog cap)"
)

## handleShow() drives the visible popup via SHOW_DELAY_MS debounce.
hs_anchor <- "function handleShow"
hs_pos <- regexpr(hs_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(hs_pos > 0L, "located handleShow body")
hs_window <- substr(helpers_src, hs_pos, min(nchar(helpers_src), hs_pos + 1500L))
stopifnot_msg(
  grepl("SHOW_DELAY_MS", hs_window, fixed = TRUE) &&
    grepl("setTimeout", hs_window, fixed = TRUE),
  "handleShow applies SHOW_DELAY_MS debounce via setTimeout"
)
stopifnot_msg(
  grepl("renderShow()", hs_window, fixed = TRUE),
  "handleShow ultimately renders the popup once debounce elapses"
)

## handleHide() cancels any pending show timer and hides the popup.
hh_anchor <- "function handleHide"
hh_pos <- regexpr(hh_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(hh_pos > 0L, "located handleHide body")
hh_window <- substr(helpers_src, hh_pos, min(nchar(helpers_src), hh_pos + 800L))
stopifnot_msg(
  grepl("clearPendingShow|clearTimeout", hh_window),
  "handleHide cancels any pending show timer"
)
## handleHide() now goes through forceHide(), the unified DOM-probing
## hide path that cannot be defeated by a drifted JS isVisible flag.
## renderHide() is kept as a thin forwarder so legacy callers still
## work.
stopifnot_msg(
  grepl("forceHide(", hh_window, fixed = TRUE) ||
    grepl("renderHide()", hh_window, fixed = TRUE),
  "handleHide hides the popup if it was visible (forceHide / renderHide)"
)

## Watchdog (rAF) enforces the absolute hard cap.
wd_anchor <- "function watchdogTick"
wd_pos <- regexpr(wd_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(wd_pos > 0L, "located watchdogTick body")
wd_window <- substr(helpers_src, wd_pos, min(nchar(helpers_src), wd_pos + 1500L))
stopifnot_msg(
  grepl("ABSOLUTE_TIMEOUT", wd_window, fixed = TRUE),
  "watchdogTick enforces ABSOLUTE_TIMEOUT"
)
stopifnot_msg(
  grepl("requestAnimationFrame(watchdogTick)", wd_window, fixed = TRUE),
  "watchdogTick reschedules via requestAnimationFrame"
)

## bind() must register the custom-message handlers and start the watchdog.
bind_anchor <- "function bind() {"
ba <- regexpr(bind_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(ba > 0L, "located bind() body")
bind_window <- substr(helpers_src, ba, min(nchar(helpers_src), ba + 3000L))
stopifnot_msg(
  grepl('Shiny.addCustomMessageHandler("phenomapr-busy-show", handleShow)',
        bind_window, fixed = TRUE),
  "bind() registers phenomapr-busy-show -> handleShow"
)
stopifnot_msg(
  grepl('Shiny.addCustomMessageHandler("phenomapr-busy-hide", handleHide)',
        bind_window, fixed = TRUE),
  "bind() registers phenomapr-busy-hide -> handleHide"
)
stopifnot_msg(
  grepl("requestAnimationFrame(watchdogTick)", bind_window, fixed = TRUE),
  "bind() starts the watchdogTick rAF chain"
)

## shiny:idle / shiny:disconnected listeners are wired up (advisory hides).
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:idle",\\s*onShinyIdle',
        helpers_src),
  "shiny:idle wired via native addEventListener (advisory)"
)
stopifnot_msg(
  grepl('document\\.addEventListener\\("shiny:disconnected",\\s*onShinyDisconnected',
        helpers_src),
  "shiny:disconnected wired via native addEventListener"
)
stopifnot_msg(
  grepl('\\$j\\(document\\)\\.on\\("shiny:idle",\\s*onShinyIdle\\)',
        helpers_src),
  "shiny:idle also wired via jQuery .on()"
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
