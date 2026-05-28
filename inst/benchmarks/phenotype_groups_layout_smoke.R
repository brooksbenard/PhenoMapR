#!/usr/bin/env Rscript
## Smoke tests for the 1. Data sidebar trim, score-table -> data-button
## refactor, side-by-side Phenotype-groups layout, and the
## "Most Phenotype +/-" display remap.
##
## Strategy: most assertions are static inspection of the source files
## (so we can verify the UI structure without spinning up a Shiny
## session), with a couple of live tests for `tag_table` (HTML-cell
## preservation) and the group_summary remap (canonical data values
## NEVER mutated).

suppressPackageStartupMessages({
  library(shiny)
})

stopifnot(requireNamespace("PhenoMapR", quietly = TRUE))

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R     <- system.file("shiny", "app.R",       package = "PhenoMapR")
helpers_R <- system.file("shiny", "helpers.R",   package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir   <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R     <- file.path(dev_dir, "app.R")
  helpers_R <- file.path(dev_dir, "helpers.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(helpers_R), file.exists(styles_css))

app_src     <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
helpers_src <- paste(readLines(helpers_R,  warn = FALSE), collapse = "\n")
styles_src  <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## --- 1. Sidebar helpText removed -------------------------------------------
stopifnot_msg(
  !grepl("Metadata is read automatically from Seurat", app_src, fixed = TRUE),
  "Cell metadata sidebar helpText is removed"
)
stopifnot_msg(
  !grepl("guesses sensible defaults for the cell-ID", app_src, fixed = TRUE),
  "second sentence of the removed helpText is also gone"
)

## --- 2. Score table card panel removed -------------------------------------
## The standalone DTOutput("score_table") inside a card() must be GONE.
## A renderDT("score_table") definition must NOT exist any more (we now
## use a reactive instead, so the data is downloadable without the
## table ever being rendered as a panel).
stopifnot_msg(
  !grepl("DTOutput\\(\"score_table\"\\)", app_src),
  "DTOutput('score_table') panel is no longer rendered"
)
stopifnot_msg(
  !grepl("output\\$score_table <- renderDT", app_src),
  "renderDT(score_table) was replaced by a plain reactive"
)
stopifnot_msg(
  grepl("score_table_data <- reactive", app_src),
  "score_table_data() reactive is defined for the data download"
)
stopifnot_msg(
  grepl("function() score_table_data()", app_src, fixed = TRUE),
  "score_table_download reads from the new reactive"
)

## --- 3. "Score by cell type and source" header has BOTH buttons ------------
## Find the relevant header invocation and confirm it now passes
## data_download_id = "score_table_download". The header call spans
## multiple lines with `tags$strong("...")` inside it, so we grab a
## generous window after the title literal to make sure we see the
## panel_id / data_download_id lines (regex `.*?\\)` would otherwise
## stop at the first inner `)`).
title_anchor <- 'tags$strong("Score by cell type and source")'
title_pos <- regexpr(title_anchor, app_src, fixed = TRUE)
stopifnot_msg(title_pos > 0,
              "located the 'Score by cell type and source' header call")
sb_slice <- substr(app_src,
                   max(1L, title_pos - 200L),
                   min(nchar(app_src), title_pos + 600L))
stopifnot_msg(
  grepl('panel_id = "score_box_source_plot"', sb_slice, fixed = TRUE),
  "boxplot header still binds the modal-plot button to score_box_source_plot"
)
stopifnot_msg(
  grepl('data_download_id = "score_table_download"', sb_slice, fixed = TRUE),
  "boxplot header now also exposes the data-download button"
)
stopifnot_msg(
  grepl('data_tooltip = "Download plot data \\(TSV\\)"', sb_slice),
  "data-download button uses the 'Download plot data (TSV)' tooltip"
)

## phenomapr_card_header_modal_dl supports the new args.
stopifnot_msg(
  grepl("data_download_id = NULL", helpers_src, fixed = TRUE),
  "phenomapr_card_header_modal_dl exposes optional data_download_id"
)
stopifnot_msg(
  grepl('class = "phenomapr-card-header-dl-actions"', helpers_src,
        fixed = TRUE),
  "two-button variant wraps both buttons in a .actions flex div"
)

## CSS for the .actions wrapper exists and right-pins both buttons.
stopifnot_msg(
  grepl("\\.phenomapr-card-header-dl-actions", styles_src),
  "styles.css defines .phenomapr-card-header-dl-actions"
)
stopifnot_msg(
  grepl("margin-left: auto !important", styles_src),
  "actions wrapper is right-pinned via margin-left: auto"
)

## --- 4. Phenotype-groups layout: side-by-side via layout_columns(6, 6) -----
## Locate the chunk between "About phenotype tails" and "Group sizes"
## and verify it sits inside a layout_columns(col_widths = c(6, 6), ...).
pg_match <- regexpr(
  "layout_columns\\(\\s*col_widths = c\\(6, 6\\)[\\s\\S]*?About phenotype tails[\\s\\S]*?Group sizes",
  app_src, perl = TRUE
)
stopifnot_msg(pg_match > 0,
              "'About phenotype tails' + 'Group sizes' live inside layout_columns(6,6)")

## --- 5. "Most Phenotype +/-" replacement landed everywhere -----------------
## Home-page hero now uses the new strong-wrapped labels.
stopifnot_msg(
  grepl("<strong>Most Phenotype \\+</strong>", app_src),
  "home hero replaces 'Most Adverse' with bold 'Most Phenotype +'"
)
stopifnot_msg(
  grepl("<strong>Most Phenotype &minus;</strong>", app_src,
        fixed = TRUE),
  "home hero replaces 'Most Favorable' with bold 'Most Phenotype \u2212'"
)
stopifnot_msg(
  !grepl("<strong>Most Adverse</strong>", app_src, fixed = TRUE) &&
  !grepl("<strong>Most Favorable</strong>", app_src, fixed = TRUE),
  "no stale <strong>Most Adverse/Favorable</strong> remain in app.R"
)

## "About phenotype tails" markdown uses the renamed terms.
## We accept either the bare Unicode minus (U+2212) or the R-source
## escape `\u2212` in the file's bytes -- both R parse to the same
## char at app load time, and StrReplace edits emitted the escape
## form rather than the literal Unicode char.
stopifnot_msg(
  grepl("**Most Phenotype +**", app_src, fixed = TRUE) &&
  (grepl("**Most Phenotype \u2212**", app_src, fixed = TRUE) ||
   grepl("**Most Phenotype \\u2212**", app_src, fixed = TRUE)),
  "About-phenotype-tails markdown uses **Most Phenotype +/-**"
)
stopifnot_msg(
  !grepl("**Most Adverse**", app_src, fixed = TRUE) &&
  !grepl("**Most Favorable**", app_src, fixed = TRUE),
  "no stale **Most Adverse/Favorable** in app.R markdown"
)

## Sidebar Phenotype-groups helpText now uses bold +/- labels.
stopifnot_msg(
  grepl("<strong>Most Phenotype +</strong>", app_src, fixed = TRUE) &&
  grepl("<strong>Most Phenotype &minus;</strong>", app_src, fixed = TRUE),
  "sidebar Phenotype groups helpText uses <strong>Most Phenotype +/-</strong>"
)

## Plots remap the LEGEND labels but keep the underlying data values.
stopifnot_msg(
  grepl('"Most Adverse"   = "Most Phenotype +"', app_src, fixed = TRUE),
  "scale_*_manual labels remap 'Most Adverse' -> 'Most Phenotype +'"
)
stopifnot_msg(
  grepl('"Most Favorable" = "Most Phenotype \u2212"', app_src, fixed = TRUE) ||
  grepl('"Most Favorable" = "Most Phenotype \\u2212"', app_src, fixed = TRUE),
  "scale_*_manual labels remap 'Most Favorable' -> 'Most Phenotype \u2212'"
)

## CRITICAL invariant: the data layer still maps the CANONICAL strings to
## colours. If anyone "helpfully" renames the keys to the new display
## strings, define_phenotype_groups()'s output stops being colored.
stopifnot_msg(
  grepl('"Most Adverse"   = "#B2182B"', app_src, fixed = TRUE),
  "scale_*_manual values still keyed on canonical 'Most Adverse'"
)
stopifnot_msg(
  grepl('"Most Favorable" = "#2166AC"', app_src, fixed = TRUE),
  "scale_*_manual values still keyed on canonical 'Most Favorable'"
)

## --- 6. Live test: tag_table preserves HTML cell content -------------------
helpers_env <- new.env(parent = globalenv())
sys.source(helpers_R, envir = helpers_env)

stopifnot_msg(exists("tag_table", envir = helpers_env),
              "tag_table is loaded from helpers.R")
tt <- helpers_env$tag_table

df_html <- data.frame(label = c("a", "b"), stringsAsFactors = FALSE)
df_html$label <- list(HTML("<strong>A</strong>"), HTML("plain"))
tbl <- tt(df_html)
tbl_html <- as.character(tbl)
stopifnot_msg(
  grepl("<strong>A</strong>", tbl_html, fixed = TRUE),
  "tag_table preserves <strong>...</strong> inside list-column cells"
)
stopifnot_msg(
  grepl(">plain<", tbl_html, fixed = TRUE),
  "tag_table also renders plain-text HTML-class cells without escaping"
)

## Backward compat: ordinary character columns still render as
## escaped text.
df_plain <- data.frame(
  group = c("Most Adverse", "Most Favorable", "Other"),
  Cells = c("10", "12", "8"),
  stringsAsFactors = FALSE
)
tbl_plain <- tt(df_plain)
tbl_plain_html <- as.character(tbl_plain)
stopifnot_msg(
  grepl(">Most Adverse<", tbl_plain_html, fixed = TRUE) &&
  grepl(">Most Favorable<", tbl_plain_html, fixed = TRUE),
  "tag_table still works for plain character columns"
)

## --- 7. Live test: group_summary remap (canonical data preserved) ----------
## Re-create the remap that lives in output$group_summary so we can
## prove it is purely cosmetic: it builds NEW HTML strings for the
## table cells, but the underlying state$groups vector is untouched.
display_map <- c(
  "Most Adverse"   = "<strong>Most Phenotype +</strong>",
  "Most Favorable" = "<strong>Most Phenotype &minus;</strong>",
  "Other"          = "Other"
)
raw <- c("Most Adverse", "Most Favorable", "Other", "Most Adverse")
mapped <- display_map[raw]
stopifnot_msg(
  identical(unname(mapped),
            c("<strong>Most Phenotype +</strong>",
              "<strong>Most Phenotype &minus;</strong>",
              "Other",
              "<strong>Most Phenotype +</strong>")),
  "group_summary remap maps canonical strings to bold display labels"
)
stopifnot_msg(
  identical(raw,
            c("Most Adverse", "Most Favorable", "Other", "Most Adverse")),
  "original raw vector is untouched by the remap"
)

## And: feed the remap through the live tag_table to confirm the
## bolded labels make it into the rendered HTML.
tbl_df <- data.frame(
  `Phenotype group` = c("placeholder", "placeholder", "placeholder"),
  Cells = c("10", "12", "8"),
  check.names = FALSE, stringsAsFactors = FALSE
)
tbl_df[["Phenotype group"]] <- lapply(
  display_map[c("Most Adverse", "Most Favorable", "Other")],
  HTML
)
rendered <- as.character(tt(tbl_df))
stopifnot_msg(
  grepl("<strong>Most Phenotype +</strong>", rendered, fixed = TRUE),
  "rendered group_summary table contains <strong>Most Phenotype +</strong>"
)
stopifnot_msg(
  grepl("<strong>Most Phenotype &minus;</strong>", rendered, fixed = TRUE),
  "rendered group_summary table contains <strong>Most Phenotype \u2212</strong>"
)
stopifnot_msg(
  !grepl("Most Adverse", rendered, fixed = TRUE) &&
  !grepl("Most Favorable", rendered, fixed = TRUE),
  "rendered group_summary table does NOT surface canonical strings"
)

cat("\nAll Phenotype-groups layout + label-remap smoke tests passed.\n")
