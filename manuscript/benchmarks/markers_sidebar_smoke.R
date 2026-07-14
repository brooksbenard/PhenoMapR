#!/usr/bin/env Rscript
## Smoke tests for the Markers-tab sidebar polish:
##
##   1. Sidebar gets an h4("Marker parameters") heading at the top,
##      mirroring the Cell metadata / Phenotype groups layout.
##   2. The radioButtons("marker_scope") labels are rebranded:
##         "Cohort-wide (adverse vs favorable)"
##            -> "Cohort-wide (phenotype + vs phenotype -)"
##         "Cell-type specific"
##            -> "Cell-type specific (within phenotype groups)"
##      The first two values stayed the same ("phenotype_groups" /
##      "cell_type_specific") so all server logic that pattern-matches
##      them still works.
##   2b. A NEW third value "cell_type_vs_opposite" is exposed by the
##      same radioButtons -- contrasts each (cell type x tail) block
##      against ALL cells in the opposite phenotype tail (any cell
##      type). This recovers markers for cell types that exist in
##      only one tail, where the "within phenotype groups" contrast
##      returns empty. The run_markers observer maps it to
##      celltype_contrast = "vs_opposite_tail".
##   3. The "Marker-gene heatmap" subsection of the sidebar is wrapped
##      in `.phenomapr-compact-stack marker-heatmap-compact-stack`, and
##      styles.css ships the supplementary `.marker-heatmap-compact-stack`
##      rules that tighten the h4 / helpText / Draw button spacing.
##   4. The Marker-gene heatmap subsection exposes a
##      checkboxInput("hm_block_borders") so users can toggle
##      visible black borders around each cell-type x phenotype
##      block. The renderImage path consults this input live and
##      passes block_outline_color = "black" / "white" to
##      plot_phenotype_markers().

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

app_R      <- system.file("shiny", "app.R",            package = "PhenoMapR")
styles_css <- system.file("shiny", "www", "styles.css", package = "PhenoMapR")
dev_dir <- Sys.getenv("PHENOMAPR_SHINY_DIR", unset = "")
if (nzchar(dev_dir) && dir.exists(dev_dir)) {
  app_R      <- file.path(dev_dir, "app.R")
  styles_css <- file.path(dev_dir, "www", "styles.css")
}
stopifnot(file.exists(app_R), file.exists(styles_css))
app_src    <- paste(readLines(app_R,      warn = FALSE), collapse = "\n")
styles_src <- paste(readLines(styles_css, warn = FALSE), collapse = "\n")

## Locate the "5. Markers" sidebar so we can scope every assertion to
## just that nav_panel. Without scoping, an "adverse vs favorable" hit
## elsewhere in the app (e.g. the Phenotype-groups help text or the
## About-markers blurb) would mask a regression in the actual label.
anchor <- '" 5. Markers")'
ap <- regexpr(anchor, app_src, fixed = TRUE)
stopifnot_msg(ap > 0, "located the '5. Markers' nav_panel anchor")
## The sidebar lives in the next ~3000 chars after the nav_step span.
slice <- substr(app_src, ap, min(nchar(app_src), ap + 5000L))

## --- 1. h4("Marker parameters") heading at the top of the sidebar ---------
stopifnot_msg(
  grepl('h4("Marker parameters")', slice, fixed = TRUE),
  "Markers sidebar opens with h4(\"Marker parameters\")"
)
## ...and it appears BEFORE radioButtons("marker_scope").
h4_pos <- regexpr('h4("Marker parameters")', slice, fixed = TRUE)
rb_pos <- regexpr('radioButtons(\n          "marker_scope"', slice, fixed = TRUE)
stopifnot_msg(h4_pos > 0 && rb_pos > 0 && h4_pos < rb_pos,
              "h4 heading sits ABOVE the marker_scope radioButtons")

## --- 2. Rebranded scope labels ---------------------------------------------
## R parses "\u2212" inside a double-quoted string as the Unicode
## minus character at runtime, but the on-disk source contains the
## literal six-character escape sequence "\u2212". readLines() does
## NOT interpret the escape, so we assert against the literal source
## form rather than the post-parse character. Both shapes are
## acceptable -- we accept either to be tolerant of file rewrites
## that re-encode the escape as the literal U+2212 glyph.
has_unicode_minus_label <- function(s) {
  grepl('"Cohort-wide (phenotype + vs phenotype \\u2212)"',
        s, fixed = TRUE) ||
  grepl(paste0('"Cohort-wide (phenotype + vs phenotype ',
               "\u2212", ')"'),
        s, fixed = TRUE) ||
  grepl('phenomapr_phenotype_plus\\("phenotype \\+"\\)', s) ||
  grepl('phenomapr_phenotype_minus\\("phenotype \\\\u2212"\\)', s)
}
ct_specific_label <- "\"Cell-type specific (within phenotype groups)\""
stopifnot_msg(
  has_unicode_minus_label(slice),
  "cohort-wide scope label rebranded to 'phenotype + vs phenotype -' (Unicode minus)"
)
stopifnot_msg(
  grepl(ct_specific_label, slice, fixed = TRUE),
  "cell-type-specific scope label clarified with '(within phenotype groups)'"
)
## Underlying values stay canonical so server-side comparisons keep working.
stopifnot_msg(
  grepl('"phenotype_groups"', slice, fixed = TRUE) &&
  grepl('"cell_type_specific"', slice, fixed = TRUE),
  "underlying marker_scope values stay 'phenotype_groups' / 'cell_type_specific'"
)
## Old wording is GONE from the Markers sidebar slice.
stopifnot_msg(
  !grepl('"Cohort-wide (adverse vs favorable)"', slice, fixed = TRUE),
  "old 'Cohort-wide (adverse vs favorable)' label is GONE from sidebar"
)
stopifnot_msg(
  !grepl('"Cell-type specific"\\s*=', slice),
  "old standalone 'Cell-type specific' label is GONE from sidebar"
)

## --- 3. Heatmap sub-section is wrapped in the compact stack ---------------
stopifnot_msg(
  grepl('class = "phenomapr-compact-stack marker-heatmap-compact-stack"',
        slice, fixed = TRUE),
  "Marker-gene heatmap subsection is wrapped in the compact stack"
)

## All six expected children live inside the wrapper, in order.
wrap_anchor <- 'class = "phenomapr-compact-stack marker-heatmap-compact-stack"'
wp <- regexpr(wrap_anchor, slice, fixed = TRUE)
stopifnot_msg(wp > 0, "located the marker-heatmap-compact-stack opening")
sub_slice <- substr(slice, wp, min(nchar(slice), wp + 3000L))
for (needle in c(
  'h4("Marker-gene heatmap")',
  'helpText(',
  'numericInput("hm_top_n"',
  'numericInput("hm_n_labels"',
  'checkboxInput("hm_block_borders"',
  'actionButton("draw_marker_heatmap"'
)) {
  stopifnot_msg(
    grepl(needle, sub_slice, fixed = TRUE),
    sprintf("compact stack contains: %s", needle)
  )
}

## --- 2b. Third marker scope option (cell_type_vs_opposite) ---------------
stopifnot_msg(
  grepl('"cell_type_vs_opposite"', slice, fixed = TRUE),
  "marker_scope radioButtons exposes 'cell_type_vs_opposite' value"
)
stopifnot_msg(
  grepl("Cell-type", slice, fixed = TRUE) &&
    grepl("opposite-tail cells", slice, fixed = TRUE),
  "marker_scope third option label mentions 'Cell-type' and 'opposite-tail cells'"
)
## Server-side decode: the run_markers observer maps the new UI value
## to celltype_contrast = "vs_opposite_tail" before calling
## find_phenotype_markers().
stopifnot_msg(
  grepl('"cell_type_vs_opposite"', app_src, fixed = TRUE) &&
    grepl('"vs_opposite_tail"', app_src, fixed = TRUE),
  "run_markers maps cell_type_vs_opposite -> vs_opposite_tail in app.R"
)

## --- 2c. Stringency-based ordering of the three marker_scope options ------
## The radioButtons should list the three options in order of *increasing
## stringency*: cohort-wide -> cell-type x opposite-tail -> cell-type
## within phenotype groups. This is the order users encounter when reading
## the sidebar top-to-bottom and matches the "About markers" card text.
v_phenotype <- regexpr('"phenotype_groups"',      slice, fixed = TRUE)
v_opposite  <- regexpr('"cell_type_vs_opposite"', slice, fixed = TRUE)
v_specific  <- regexpr('"cell_type_specific"',    slice, fixed = TRUE)
stopifnot_msg(
  v_phenotype > 0 && v_opposite > 0 && v_specific > 0,
  "all three marker_scope values are present in the sidebar slice"
)
stopifnot_msg(
  v_phenotype < v_opposite,
  "stringency order: 'phenotype_groups' precedes 'cell_type_vs_opposite'"
)
stopifnot_msg(
  v_opposite < v_specific,
  paste0("stringency order: 'cell_type_vs_opposite' (mid) precedes ",
         "'cell_type_specific' (most stringent)")
)
## And the "About markers" card mirrors that ordering: the
## cell_type_vs_opposite description ("opposite-tail cells") appears
## before the within-phenotype-groups description in the markdown.
about_idx <- regexpr("About markers", app_src, fixed = TRUE)
stopifnot(about_idx > 0L)
about_window <- substr(app_src, about_idx,
                       min(nchar(app_src), about_idx + 4000L))
a_opposite <- regexpr("opposite-tail cells", about_window, fixed = TRUE)
a_specific <- regexpr("within phenotype groups",
                      about_window, fixed = TRUE)
stopifnot_msg(
  a_opposite > 0 && a_specific > 0,
  "About markers card describes both per-cell-type contrasts"
)
stopifnot_msg(
  a_opposite < a_specific,
  paste0("About markers card mirrors the sidebar ordering: ",
         "opposite-tail (mid) before within-phenotype-groups (most stringent)")
)

## --- 4. Live block-border toggle plumbing --------------------------------
## renderImage / download paths consult input$hm_block_borders and
## override outline_marker_blocks + block_outline_color before
## delegating to plot_phenotype_markers(). The checkbox flip therefore
## immediately re-renders the heatmap with thick black outlines around
## each cell-type x phenotype block, without rebuilding the cached
## expression matrix.
stopifnot_msg(
  grepl("input\\$hm_block_borders", app_src),
  "app.R reads input$hm_block_borders to drive the live border toggle"
)
stopifnot_msg(
  grepl("block_outline_color", app_src, fixed = TRUE),
  "app.R passes block_outline_color through to plot_phenotype_markers()"
)
stopifnot_msg(
  grepl("outline_marker_blocks", app_src, fixed = TRUE),
  "app.R passes outline_marker_blocks through to plot_phenotype_markers()"
)

## --- 4. CSS supplementary rules exist --------------------------------------
for (sel in c(
  ".marker-heatmap-compact-stack > h4",
  ".marker-heatmap-compact-stack .help-block",
  ".marker-heatmap-compact-stack > .btn",
  ".marker-heatmap-compact-stack > .action-button"
)) {
  stopifnot_msg(
    grepl(sel, styles_src, fixed = TRUE),
    sprintf("styles.css defines %s", sel)
  )
}

cat("\nAll Markers sidebar polish smoke tests passed.\n")
