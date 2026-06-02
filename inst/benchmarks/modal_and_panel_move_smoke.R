#!/usr/bin/env Rscript
## Smoke tests for four polish changes:
##
##   1. Sidebar h4 headings ("Expression input", "Cell metadata",
##      "Embedding", "Marker-gene heatmap") no longer carry an
##      explicit `font-size: 1rem` override -- they inherit the
##      default h4 size so every sidebar reads at a consistent
##      larger weight.
##   2. phenomapr_plot_download_modal() now puts Cancel + Download
##      in the title row (top-right) and clears the footer.
##   3. The modal exposes a new `plot_dl_units` radio (inches / cm)
##      and the server-side preview + downloadHandler convert cm
##      to inches before passing widths/heights to ggsave / device.
##   4. The "Score by cell type and source" and "Per-cell-type group
##      enrichment" cards have moved from the 3. Score tab into the
##      4. Visualization tab (UI placement only -- the renderPlot /
##      register_plot_download wiring is unchanged).

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R     <- system.file("shiny", "app.R",            package = "PhenoMapR")
helpers_R <- system.file("shiny", "helpers.R",        package = "PhenoMapR")
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

extract_rule <- function(css, selector) {
  pattern <- paste0("(?m)^", selector, "\\s*\\{[^}]*\\}")
  regmatches(css, regexpr(pattern, css, perl = TRUE))
}

## --- 1. Sidebar h4 headings inherit default size --------------------------
## None of the four compact-stack h4 rules should still carry an
## explicit font-size declaration -- only margin trims. Removing the
## explicit `font-size: 1rem` is the regression-fix for the "headings
## look small" complaint.
for (sel in c("\\.metadata-compact-stack > h4",
              "\\.expression-compact-stack > h4",
              "\\.embedding-compact-stack > h4",
              "\\.marker-heatmap-compact-stack > h4")) {
  rule <- extract_rule(styles_src, sel)
  stopifnot_msg(length(rule) == 1L,
                sprintf("located rule %s", gsub("\\\\", "", sel)))
  stopifnot_msg(
    !grepl("font-size:", rule),
    sprintf("%s no longer carries an explicit font-size declaration",
            gsub("\\\\", "", sel))
  )
  ## Margin trims should still be present.
  stopifnot_msg(
    grepl("margin-top: 0", rule) || grepl("margin-bottom:", rule),
    sprintf("%s still trims margins", gsub("\\\\", "", sel))
  )
}

## --- 2. Plot-download modal: Cancel + Download in title row ---------------
stopifnot_msg(
  grepl('phenomapr-plot-dl-modal-title-row', helpers_src, fixed = TRUE),
  "modal title is now a flex .phenomapr-plot-dl-modal-title-row"
)
stopifnot_msg(
  grepl('phenomapr-plot-dl-modal-actions', helpers_src, fixed = TRUE),
  "modal title row contains a right-aligned actions cluster"
)
stopifnot_msg(
  grepl("footer = NULL", helpers_src, fixed = TRUE),
  "modal explicitly disables the default footer (footer = NULL)"
)
## Cancel + Download buttons must live inside the actions cluster,
## not as a separate `footer = tagList(...)` argument.
title_anchor <- "phenomapr-plot-dl-modal-actions"
ta <- regexpr(title_anchor, helpers_src, fixed = TRUE)
stopifnot_msg(ta > 0, "located the modal actions cluster")
title_slice <- substr(helpers_src, ta,
                      min(nchar(helpers_src), ta + 1500L))
stopifnot_msg(
  grepl('shiny::modalButton\\("Cancel"\\)', title_slice),
  "actions cluster includes modalButton('Cancel')"
)
stopifnot_msg(
  grepl('shiny::downloadButton\\(\\s*"plot_dl_action"', title_slice),
  "actions cluster includes the plot_dl_action downloadButton"
)
## CSS rules for the title row + actions cluster exist.
for (sel in c(".phenomapr-plot-dl-modal-title-row",
              ".phenomapr-plot-dl-modal-title-text",
              ".phenomapr-plot-dl-modal-actions")) {
  stopifnot_msg(grepl(sel, styles_src, fixed = TRUE),
                sprintf("styles.css defines %s", sel))
}

## --- 3. Inches / cm unit selector ----------------------------------------
stopifnot_msg(
  grepl('"plot_dl_units"', helpers_src, fixed = TRUE),
  "modal exposes radioButtons('plot_dl_units')"
)
stopifnot_msg(
  grepl('"inches" = "in"', helpers_src, fixed = TRUE) &&
  grepl('"cm" = "cm"', helpers_src, fixed = TRUE),
  "units selector has inches + cm choices"
)
## The modal label helper uses the current unit.
stopifnot_msg(
  grepl('unit_lbl <- function\\(side\\)', helpers_src),
  "helper builds unit-aware Width/Height labels"
)
## Server-side wiring: prev_dl_units reactiveVal + cm rescale.
stopifnot_msg(
  grepl("prev_dl_units <- reactiveVal", app_src, fixed = TRUE),
  "server tracks previous unit via prev_dl_units reactiveVal"
)
stopifnot_msg(
  grepl(".CM_PER_INCH <- 2.54", app_src, fixed = TRUE),
  "server exposes a 2.54 cm/inch constant"
)
stopifnot_msg(
  grepl(".to_inches <- function", app_src, fixed = TRUE),
  "server defines .to_inches() helper"
)
## The downloadHandler reads units AND converts width/height to inches.
dh_anchor <- 'output$plot_dl_action <- downloadHandler'
dh_pos <- regexpr(dh_anchor, app_src, fixed = TRUE)
stopifnot_msg(dh_pos > 0, "located the plot_dl_action downloadHandler")
dh_slice <- substr(app_src, dh_pos,
                   min(nchar(app_src), dh_pos + 3000L))
stopifnot_msg(
  grepl('isolate\\(input\\$plot_dl_units\\)', dh_slice),
  "downloadHandler reads input$plot_dl_units"
)
stopifnot_msg(
  grepl("w <- .to_inches\\(w, units\\)", dh_slice, fixed = FALSE) &&
  grepl("h <- .to_inches\\(h, units\\)", dh_slice, fixed = FALSE),
  "downloadHandler converts w/h to inches before passing to ggsave/devices"
)
## The preview helper also converts.
pp_anchor <- ".plot_dl_preview_dims <- function"
pp_pos <- regexpr(pp_anchor, app_src, fixed = TRUE)
stopifnot_msg(pp_pos > 0, "located the .plot_dl_preview_dims helper")
pp_slice <- substr(app_src, pp_pos,
                   min(nchar(app_src), pp_pos + 1500L))
stopifnot_msg(
  grepl("w_in <- .to_inches", pp_slice, fixed = TRUE) &&
  grepl("h_in <- .to_inches", pp_slice, fixed = TRUE),
  "preview helper converts w/h to inches before computing pixel box"
)

## --- 4. "Score by cell type and source" + "Per-cell-type group
##        enrichment" panels moved to Visualization tab --------------------
## They must NOT appear before the "4. Visualization" anchor.
v_anchor <- '" 4. Visualization"'
vp <- regexpr(v_anchor, app_src, fixed = TRUE)
stopifnot_msg(vp > 0, "located the '4. Visualization' nav anchor")
pre_viz  <- substr(app_src, 1L, vp - 1L)
post_viz <- substr(app_src, vp, nchar(app_src))

## The plotOutput("score_box_source_plot", ...) and
## plotOutput("group_by_celltype_plot", ...) must each appear EXACTLY
## ONCE total in the whole source (no duplication during the move)
## AND they must live POST the Visualization anchor.
score_box_n <- length(gregexpr('plotOutput\\("score_box_source_plot"',
                               app_src)[[1L]])
group_ct_n  <- length(gregexpr('plotOutput\\("group_by_celltype_plot"',
                               app_src)[[1L]])
stopifnot_msg(score_box_n == 1L,
              sprintf("plotOutput('score_box_source_plot') appears exactly once (got %d)",
                      score_box_n))
stopifnot_msg(group_ct_n == 1L,
              sprintf("plotOutput('group_by_celltype_plot') appears exactly once (got %d)",
                      group_ct_n))
stopifnot_msg(
  grepl('plotOutput\\("score_box_source_plot"', post_viz),
  "score_box_source_plot now lives AFTER the Visualization anchor"
)
stopifnot_msg(
  grepl('plotOutput\\("group_by_celltype_plot"', post_viz),
  "group_by_celltype_plot now lives AFTER the Visualization anchor"
)
stopifnot_msg(
  !grepl('plotOutput\\("score_box_source_plot"', pre_viz),
  "score_box_source_plot no longer appears BEFORE the Visualization anchor"
)
stopifnot_msg(
  !grepl('plotOutput\\("group_by_celltype_plot"', pre_viz),
  "group_by_celltype_plot no longer appears BEFORE the Visualization anchor"
)
## Sanity: both card_headers also live in the Visualization slice.
stopifnot_msg(
  grepl('"Score by cell type and source"', post_viz, fixed = TRUE),
  "Visualization slice carries the 'Score by cell type and source' header"
)
stopifnot_msg(
  grepl('"Per-cell-type group enrichment"', post_viz, fixed = TRUE),
  "Visualization slice carries the 'Per-cell-type group enrichment' header"
)
## The migrated boxplot still has its data-download + plot buttons.
stopifnot_msg(
  grepl('data_download_id = "score_table_download"', post_viz, fixed = TRUE),
  "boxplot's data-download wiring survived the move"
)

cat("\nAll modal + panel-move smoke tests passed.\n")
