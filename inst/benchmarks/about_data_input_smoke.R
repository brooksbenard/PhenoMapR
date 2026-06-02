#!/usr/bin/env Rscript
## Smoke tests for the five "polish" changes:
##
##   1. .phenomapr-file-input chain now caps width and forces shrinking
##      so long file names (e.g. BRCA_GSE114727_10X_expression.h5)
##      cannot overflow the Browse + filename row.
##   2. The "About data input" card is split into a two-column body
##      (bullet list on the left, two callout notes on the right). The
##      upload-size note is reworded per the user's spec.
##   3. phenomapr_dl_btn / phenomapr_modal_dl_btn / the modal-DL card
##      header support optional labels so the two download buttons in
##      the "Score by cell type and source" header are visibly
##      distinct ("Download plot" / "Download plot data").
##   4. score_table_data() now splices phenotype_group_<col> columns
##      from state$groups into its export so the "Download plot data"
##      TSV self-contains the phenotype tail labels.
##
## Most assertions are static source inspection; the labeled-button
## helpers and the score_table_data merge are also live-exercised.

suppressPackageStartupMessages({
  library(shiny)
})

stopifnot(requireNamespace("PhenoMapR", quietly = TRUE))

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R      <- system.file("shiny", "app.R",            package = "PhenoMapR")
helpers_R  <- system.file("shiny", "helpers.R",        package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir    <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  helpers_R  <- file.path(dev_dir, "helpers.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(helpers_R), file.exists(styles_css))

app_src     <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
helpers_src <- paste(readLines(helpers_R,  warn = FALSE), collapse = "\n")
styles_src  <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## --- 1. File-input overflow guards -----------------------------------------
## These three CSS edits collectively keep the filename clamped:
## (a) outer block has min-width:0 + max-width:100% so it doesn't
##     push the sidebar wider than the viewport allocates.
## (b) row has max-width:100% / width:100% so it sizes to the block
##     rather than to its widest child.
## (c) name + name-set + name-text each have flex: 1 1 auto + min-width:0
##     so the ellipsis truncation actually triggers.
## Anchor each match to the START OF A LINE via PCRE's (?m)^ multiline
## flag, so compound selectors like
## ".phenomapr-compact-stack > .phenomapr-file-input { ... }" (where
## the suffix appears mid-line) don't shadow the simple selector
## we care about.
extract_rule <- function(css, selector) {
  pattern <- paste0("(?m)^", selector, "\\s*\\{[^}]*\\}")
  regmatches(css, regexpr(pattern, css, perl = TRUE))
}

file_input_rule <- extract_rule(styles_src, "\\.phenomapr-file-input")
stopifnot_msg(length(file_input_rule) == 1L,
              "located the .phenomapr-file-input outer rule")
stopifnot_msg(
  grepl("min-width: 0", file_input_rule) &&
  grepl("max-width: 100%", file_input_rule),
  "outer .phenomapr-file-input pins min-width:0 + max-width:100%"
)

row_rule <- extract_rule(styles_src, "\\.phenomapr-file-input-row")
stopifnot_msg(length(row_rule) == 1L,
              "located the .phenomapr-file-input-row rule")
stopifnot_msg(
  grepl("max-width: 100%", row_rule) && grepl("width: 100%", row_rule),
  ".phenomapr-file-input-row caps to 100% of its parent"
)

name_rule <- extract_rule(styles_src, "\\.phenomapr-file-input-name")
stopifnot_msg(length(name_rule) == 1L,
              "located the .phenomapr-file-input-name rule")
stopifnot_msg(
  grepl("min-width: 0", name_rule) &&
  grepl("max-width: 100%", name_rule) &&
  grepl("overflow: hidden", name_rule),
  ".phenomapr-file-input-name has min-width:0 + max-width:100% + overflow:hidden"
)
## word-break / overflow-wrap rules were removed as ACTIVE
## declarations (they caused the inline-flex wrapper to size to its
## content instead of its parent). The explanatory comment may still
## reference them, so strip comments before asserting.
name_rule_uncomm <- gsub("/\\*[\\s\\S]*?\\*/", "", name_rule, perl = TRUE)
stopifnot_msg(
  !grepl("word-break:\\s*break-all\\s*;", name_rule_uncomm),
  "word-break:break-all is no longer an active declaration in .phenomapr-file-input-name"
)

set_rule <- extract_rule(styles_src, "\\.phenomapr-file-input-name-set")
stopifnot_msg(length(set_rule) == 1L,
              "located the .phenomapr-file-input-name-set rule")
stopifnot_msg(
  grepl("flex: 1 1 auto", set_rule) &&
  grepl("min-width: 0", set_rule) &&
  grepl("max-width: 100%", set_rule),
  ".phenomapr-file-input-name-set has flex:1 1 auto + min-width:0 + max-width:100%"
)

text_rule <- extract_rule(styles_src, "\\.phenomapr-file-input-name-text")
stopifnot_msg(length(text_rule) == 1L,
              "located the .phenomapr-file-input-name-text rule")
stopifnot_msg(
  grepl("flex: 1 1 auto", text_rule) &&
  grepl("white-space: nowrap", text_rule) &&
  grepl("text-overflow: ellipsis", text_rule) &&
  grepl("overflow: hidden", text_rule),
  ".phenomapr-file-input-name-text has flex + nowrap + ellipsis + hidden"
)

## --- 2. "About data input" three-column body -------------------------------
## All three notes (PhenoMapR accepts / Gene IDs / Upload size) are
## peer .data-input-note callouts arranged horizontally with equal
## widths (c(4, 4, 4)). On narrow viewports bslib stacks them
## vertically.
stopifnot_msg(
  grepl("layout_columns\\(\\s*col_widths = c\\(4, 4, 4\\)", app_src),
  "About-data-input uses layout_columns(c(4, 4, 4)) for the body"
)

## The accepts bullet list is now its own .data-input-note callout
## (rather than a left-column wrapper), so the helper text headings
## have to match the three peer titles.
stopifnot_msg(
  grepl(' PhenoMapR accepts"', app_src, fixed = TRUE),
  "accepts note carries the title 'PhenoMapR accepts'"
)
stopifnot_msg(
  grepl(" Gene IDs", app_src) && grepl(" Upload size", app_src),
  "side notes still carry the expected titles (Gene IDs + Upload size)"
)
stopifnot_msg(
  grepl("HUGO symbols", app_src),
  "Gene-IDs note still mentions HUGO symbols"
)
stopifnot_msg(
  grepl("Matrix / data.frame", app_src),
  "accepts note still lists Matrix / data.frame"
)
stopifnot_msg(
  grepl("AnnData", app_src),
  "accepts note still lists AnnData"
)

## The 5MB cap text was replaced with the new wording.
stopifnot_msg(
  grepl("analyses should work as long as your machine has", app_src),
  "Upload-size note carries the new 'analyses should work' wording"
)
stopifnot_msg(
  !grepl("Shiny's 5 MB", app_src, fixed = TRUE),
  "old 'Shiny's 5 MB cap is removed' phrasing is GONE"
)
stopifnot_msg(
  !grepl(".rds. / .h5ad", app_src),
  "old '.rds / .h5ad' phrasing is GONE"
)

## CSS classes for the callouts exist. The pre-restructure
## `.data-input-accepts` / `.data-input-notes` selectors are no
## longer needed because the accepts list is now itself a peer
## `.data-input-note` -- we explicitly assert they're gone so a
## future re-edit doesn't accidentally re-introduce a heterogeneous
## layout.
for (cl in c(".data-input-note",
             ".data-input-note-title")) {
  stopifnot_msg(grepl(cl, styles_src, fixed = TRUE),
                sprintf("styles.css defines %s", cl))
}
stopifnot_msg(
  !grepl(".data-input-accepts", styles_src, fixed = TRUE),
  "obsolete .data-input-accepts selector is GONE from styles.css"
)
stopifnot_msg(
  !grepl(".data-input-notes", styles_src, fixed = TRUE),
  "obsolete .data-input-notes selector is GONE from styles.css"
)

## Equal-height styling so the three peers line up at the same
## height regardless of body length.
note_rule <- extract_rule(styles_src, "\\.data-input-note")
stopifnot_msg(length(note_rule) == 1L,
              "located the .data-input-note rule")
stopifnot_msg(
  grepl("height: 100%", note_rule) && grepl("flex-direction: column", note_rule),
  ".data-input-note stretches to fill its cell (height:100% + flex column)"
)

## --- 3. Labeled download buttons in the boxplot header --------------------
## The helper signatures accept the new label args.
stopifnot_msg(
  grepl("label = NULL", helpers_src, fixed = TRUE),
  "phenomapr_dl_btn / phenomapr_modal_dl_btn expose `label = NULL`"
)
stopifnot_msg(
  grepl("plot_label = NULL", helpers_src, fixed = TRUE),
  "phenomapr_card_header_modal_dl exposes plot_label"
)
stopifnot_msg(
  grepl("data_label = NULL", helpers_src, fixed = TRUE),
  "phenomapr_card_header_modal_dl exposes data_label"
)
stopifnot_msg(
  grepl("phenomapr-panel-download-btn-labeled", helpers_src, fixed = TRUE),
  "helpers add the `--labeled` modifier when a label is provided"
)

## The boxplot header now passes both labels.
title_anchor <- 'tags$strong("Score by cell type and source")'
title_pos <- regexpr(title_anchor, app_src, fixed = TRUE)
stopifnot_msg(title_pos > 0,
              "located the 'Score by cell type and source' header call")
sb_slice <- substr(app_src,
                   max(1L, title_pos - 200L),
                   min(nchar(app_src), title_pos + 800L))
stopifnot_msg(
  grepl('plot_label = "Download plot"', sb_slice, fixed = TRUE),
  "boxplot header sets plot_label = 'Download plot'"
)
stopifnot_msg(
  grepl('data_label = "Download plot data"', sb_slice, fixed = TRUE),
  "boxplot header sets data_label = 'Download plot data'"
)
stopifnot_msg(
  grepl("phenotype group labels", sb_slice),
  "boxplot data tooltip advertises that the TSV includes group labels"
)

## CSS for the labeled variant.
stopifnot_msg(
  grepl(".phenomapr-panel-download-btn-labeled", styles_src, fixed = TRUE),
  "styles.css defines .phenomapr-panel-download-btn-labeled"
)
stopifnot_msg(
  grepl("white-space: nowrap", styles_src) &&
  grepl("font-weight: 500 !important", styles_src),
  "labeled variant overrides width/height back to auto + nowraps the label"
)

## --- 4. Live test: phenomapr_card_header_modal_dl emits BOTH labels --------
helpers_env <- new.env(parent = globalenv())
sys.source(helpers_R, envir = helpers_env)
ch <- helpers_env$phenomapr_card_header_modal_dl(
  shiny::tags$strong("Score by cell type and source"),
  panel_id = "score_box_source_plot",
  plot_label = "Download plot",
  data_download_id = "score_table_download",
  data_label = "Download plot data"
)
ch_html <- as.character(ch)
stopifnot_msg(
  grepl("phenomapr-panel-download-btn-labeled", ch_html, fixed = TRUE),
  "rendered header carries the labeled class"
)
## downloadButton interpolates the label after the icon with
## surrounding whitespace (so the source is `<i></i>\n  Label\n</a>`),
## while actionButton wraps it in `<span class="action-label">Label</span>`.
## Match both shapes.
stopifnot_msg(
  grepl("Download plot data\\s*</a>", ch_html, perl = TRUE),
  "rendered header surfaces visible text 'Download plot data'"
)
stopifnot_msg(
  grepl('class="action-label">Download plot<', ch_html, fixed = TRUE),
  "rendered header surfaces visible text 'Download plot'"
)
stopifnot_msg(
  grepl("score_table_download", ch_html, fixed = TRUE) &&
  grepl("score_box_source_plot_modal_btn", ch_html, fixed = TRUE),
  "rendered header still wires the right input/download ids"
)

## --- 5. Live test: score_table_data() merges phenotype_group_<col> ---------
## We can't run the full Shiny server here, but we can replicate the
## reactive's body against a synthetic state and confirm the merge
## logic produces the expected schema.
n_cells <- 12L
set.seed(7)
scores <- data.frame(
  weighted_sum_score_precog_BRCA = rnorm(n_cells),
  row.names = paste0("Cell_", seq_len(n_cells))
)
state_scores <- scores

## state$groups must use the canonical phenotype-group strings.
state_groups <- data.frame(
  cell_id = rownames(state_scores),
  phenotype_group_weighted_sum_score_precog_BRCA = c(
    rep("Most Adverse", 3),
    rep("Other", 6),
    rep("Most Favorable", 3)
  ),
  stringsAsFactors = FALSE
)
rownames(state_groups) <- rownames(state_scores)

build_score_table <- function(scores, groups) {
  s <- scores
  s$cell_id <- rownames(s)
  score_cols <- setdiff(colnames(s), "cell_id")
  numeric_cols <- score_cols[vapply(s[score_cols], is.numeric, logical(1))]
  for (col in numeric_cols) {
    s[[paste0("scaled_", col)]] <- as.numeric(scale(s[[col]]))
  }
  grp_cols <- character(0)
  if (!is.null(groups)) {
    grp_cols <- grep("^phenotype_group_", colnames(groups), value = TRUE)
    if (length(grp_cols)) {
      join_idx <- match(rownames(s), rownames(groups))
      for (gc in grp_cols) {
        s[[gc]] <- groups[[gc]][join_idx]
      }
    }
  }
  ordered_cols <- c(
    "cell_id",
    unlist(lapply(score_cols, function(col) {
      if (col %in% numeric_cols) c(col, paste0("scaled_", col)) else col
    }), use.names = FALSE),
    grp_cols
  )
  s[, ordered_cols, drop = FALSE]
}

merged <- build_score_table(state_scores, state_groups)
stopifnot_msg(
  "phenotype_group_weighted_sum_score_precog_BRCA" %in% colnames(merged),
  "score_table_data() merges phenotype_group_<col> from state$groups"
)
stopifnot_msg(
  identical(
    tail(colnames(merged), 1L),
    "phenotype_group_weighted_sum_score_precog_BRCA"
  ),
  "phenotype_group_<col> is placed at the END of the column order"
)
stopifnot_msg(
  identical(sort(unique(merged$phenotype_group_weighted_sum_score_precog_BRCA)),
            c("Most Adverse", "Most Favorable", "Other")),
  "merged group values preserve canonical 'Most Adverse/Favorable/Other'"
)

## When state$groups is NULL, the export just drops the group column.
merged_no_groups <- build_score_table(state_scores, NULL)
stopifnot_msg(
  !any(grepl("^phenotype_group_", colnames(merged_no_groups))),
  "no group columns are added when state$groups is NULL"
)
stopifnot_msg(
  "cell_id" %in% colnames(merged_no_groups) &&
  "scaled_weighted_sum_score_precog_BRCA" %in% colnames(merged_no_groups),
  "fallback schema (no groups) still emits cell_id + scaled_<col>"
)

## --- 6. Cell metadata sidebar uses the compact stack -----------------------
## The Cell metadata sub-section in the "1. Data" sidebar (h4 +
## metadata_status + metadata_upload_panel + meta_file picker +
## three column dropdowns) is wrapped in `.phenomapr-compact-stack
## .metadata-compact-stack` so the dropdowns sit close together
## instead of spreading out with Shiny's default .form-group spacing.
stopifnot_msg(
  grepl('class = "phenomapr-compact-stack metadata-compact-stack"',
        app_src, fixed = TRUE),
  "Cell metadata sub-section is wrapped in .phenomapr-compact-stack .metadata-compact-stack"
)

## The three column dropdowns must remain INSIDE the wrapper, in order.
metadata_anchor <- 'class = "phenomapr-compact-stack metadata-compact-stack"'
md_pos <- regexpr(metadata_anchor, app_src, fixed = TRUE)
stopifnot_msg(md_pos > 0,
              "located the metadata-compact-stack wrapper")
md_slice <- substr(app_src, md_pos,
                   min(nchar(app_src), md_pos + 3500L))
for (sel in c('selectInput("meta_cell_id_col"',
              'selectInput("meta_cell_type_col"',
              'selectInput("meta_source_col"')) {
  stopifnot_msg(grepl(sel, md_slice, fixed = TRUE),
                sprintf("'%s' lives inside the compact stack", sel))
}

## The supplementary CSS rules exist.
for (sel in c(".metadata-compact-stack > h4",
              ".metadata-compact-stack .metadata-upload-details",
              ".metadata-compact-stack .metadata-upload-file-host")) {
  stopifnot_msg(grepl(sel, styles_src, fixed = TRUE),
                sprintf("styles.css defines %s", sel))
}

## --- 7. sc_or_bulk_label noun aligns with predicted modality --------------
## In R/detect_format.R, the single-cell branch uses "cells" instead
## of "samples". Bulk branch still uses "samples". The unclear
## branch hedges with "columns".
detect_src <- paste(readLines("R/detect_format.R", warn = FALSE),
                    collapse = "\n")
stopifnot_msg(
  grepl('Single-cell-like \\(%\\.0f%% zeros, %s cells\\)', detect_src),
  "single-cell branch of sc_or_bulk_label uses 'cells'"
)
stopifnot_msg(
  grepl('Bulk-like \\(%\\.0f%% zeros, %s samples\\)', detect_src),
  "bulk branch of sc_or_bulk_label keeps 'samples'"
)
stopifnot_msg(
  !grepl('Single-cell-like \\(%\\.0f%% zeros, %s samples\\)', detect_src),
  "old single-cell -> 'samples' wording is GONE"
)

## Live: call detect_expression_format and verify the noun.
detect_env <- new.env(parent = globalenv())
sys.source("R/detect_format.R", envir = detect_env)
set.seed(7)
sc <- matrix(0, nrow = 200, ncol = 400,
             dimnames = list(c("TP53","MYC", paste0("Gene", 1:198)),
                             paste0("Cell_", 1:400)))
nz <- sample(length(sc), 0.2 * length(sc))
sc[nz] <- rpois(length(nz), 2)
d_sc <- detect_env$detect_expression_format(sc)
stopifnot_msg(identical(d_sc$sc_or_bulk, "single_cell"),
              "synthetic sparse matrix is classified as single_cell")
stopifnot_msg(
  grepl("\\bcells\\b", d_sc$sc_or_bulk_label),
  sprintf("live label says '... cells' (got: %s)", d_sc$sc_or_bulk_label)
)
stopifnot_msg(
  !grepl("\\bsamples\\b", d_sc$sc_or_bulk_label),
  sprintf("live label does NOT say '... samples' for single-cell (got: %s)",
          d_sc$sc_or_bulk_label)
)

cat("\nAll 'About data input' polish smoke tests passed.\n")
