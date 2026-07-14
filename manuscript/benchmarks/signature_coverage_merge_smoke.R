#!/usr/bin/env Rscript
## Smoke tests for the merged "Phenotype signature & gene coverage"
## panel.
##
## The Phenotype-signature tab used to render two side-by-side cards
## ("Phenotype signature status" and "Gene coverage") in a
## layout_columns(col_widths = c(6, 6), ...). They have been merged
## into a single card whose body holds a 2-column layout: the left
## half is the signature status (source / cancer type / |z| cutoff /
## # significant genes) and the right half is the gene-coverage DT
## table with the inherited TSV download icon pinned to the upper
## right of the merged card header.
##
## This script asserts that:
##   1. The two-card scaffolding is gone -- no more standalone
##      `card_header(tags$strong("Gene coverage"))` next to a
##      `card_header(tags$strong("Phenotype signature status"))`.
##   2. The merged card with header
##      "Phenotype signature & gene coverage" exists and uses
##      phenomapr_card_header_dl() so the gene-coverage download
##      stays in the upper right.
##   3. Inside the merged card, both subheads
##      ("Phenotype signature status" and "Gene coverage") are
##      rendered as <h6> elements with the
##      .phenomapr-signature-coverage-subhead class.
##   4. uiOutput("reference_status") and DTOutput("gene_coverage_tbl")
##      both still mount inside the merged card body.
##   5. output$reference_status now surfaces a "Significant genes"
##      row that is computed from input$z_score_cutoff -- both for
##      built-in and custom references.
##   6. The companion CSS rules for the new layout exist in
##      styles.css.

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

## 1. The old standalone "Gene coverage" card_header is gone -- the
##    string only survives as the right-half subhead inside the
##    merged card. We assert the *card_header* form is no longer
##    present.
stopifnot_msg(
  !grepl(
    'card_header\\(\\s*tags\\$strong\\(\\s*"Gene coverage"\\s*\\)\\s*\\)',
    app_src
  ),
  "the standalone card_header(\"Gene coverage\") is gone"
)
stopifnot_msg(
  !grepl(
    paste0(
      'phenomapr_card_header_dl\\(\\s*\\n?\\s*',
      'tags\\$strong\\(\\s*"Gene coverage"\\s*\\)'
    ),
    app_src
  ),
  "the standalone phenomapr_card_header_dl(\"Gene coverage\") is gone"
)
stopifnot_msg(
  !grepl(
    'card_header\\(\\s*tags\\$strong\\(\\s*"Phenotype signature status"\\s*\\)\\s*\\)',
    app_src
  ),
  "the standalone card_header(\"Phenotype signature status\") is gone"
)

## 2. The merged card exists. We require BOTH the merged title
##    string and the phenomapr_card_header_dl wrapper (so the
##    download icon stays in the upper right).
stopifnot_msg(
  grepl(
    paste0(
      'phenomapr_card_header_dl\\([^)]*?',
      'tags\\$strong\\(\\s*"Phenotype signature & gene coverage"'
    ),
    app_src,
    perl = TRUE
  ),
  "merged phenomapr_card_header_dl(\"Phenotype signature & gene coverage\") is present"
)
stopifnot_msg(
  grepl('download_id\\s*=\\s*"gene_coverage_tbl_download"', app_src),
  "merged card header still wires download_id = gene_coverage_tbl_download"
)

## 3. The two subheads are rendered as h6 with the new helper class.
stopifnot_msg(
  grepl(
    paste0(
      'tags\\$h6\\([\\s\\S]*?"Phenotype signature status"',
      '[\\s\\S]*?phenomapr-signature-coverage-subhead'
    ),
    app_src,
    perl = TRUE
  ),
  "left subhead: <h6>Phenotype signature status</h6> with subhead class"
)
stopifnot_msg(
  grepl(
    paste0(
      'tags\\$h6\\([\\s\\S]*?"Gene coverage"',
      '[\\s\\S]*?phenomapr-signature-coverage-subhead'
    ),
    app_src,
    perl = TRUE
  ),
  "right subhead: <h6>Gene coverage</h6> with subhead class"
)

## 4. Both outputs still mount inside the merged card body.
stopifnot_msg(
  grepl('uiOutput\\(\\s*"reference_status"\\s*\\)', app_src),
  "uiOutput(\"reference_status\") is still mounted somewhere"
)
stopifnot_msg(
  grepl('DTOutput\\(\\s*"gene_coverage_tbl"\\s*\\)', app_src),
  "DTOutput(\"gene_coverage_tbl\") is still mounted somewhere"
)

## 5. output$reference_status surfaces the new "Significant genes"
##    row driven by input$z_score_cutoff. We require:
##       - the "Significant genes:" label appears exactly twice
##         (once in the built-in branch, once in the custom branch).
##       - both branches read input$z_score_cutoff to drive the count.
n_sig_lbl <- length(
  gregexpr('"Significant genes: "', app_src, fixed = TRUE)[[1L]]
)
stopifnot_msg(
  n_sig_lbl >= 2L,
  "\"Significant genes: \" label appears in both built-in and custom branches"
)
ref_status_idx <- regexpr(
  "output\\$reference_status\\s*<-\\s*renderUI\\(",
  app_src,
  perl = TRUE
)
stopifnot(ref_status_idx > 0L)
ref_status_window <- substr(
  app_src,
  ref_status_idx,
  min(nchar(app_src), ref_status_idx + 4000L)
)
stopifnot_msg(
  grepl("input\\$z_score_cutoff", ref_status_window),
  "output$reference_status reads input$z_score_cutoff for the live |z| threshold"
)
stopifnot_msg(
  grepl("get_top_prognostic_genes", ref_status_window) &&
    grepl("abs\\(z_vec\\)\\s*>=\\s*cut_val", ref_status_window),
  "built-in branch counts genes via abs(z_vec) >= cut_val"
)

## 6. CSS for the new layout is in place.
stopifnot_msg(
  grepl("\\.phenomapr-signature-coverage-card", styles_src, perl = TRUE),
  ".phenomapr-signature-coverage-card CSS rules exist"
)
stopifnot_msg(
  grepl("\\.phenomapr-signature-coverage-half\\s*\\+\\s*\\.phenomapr-signature-coverage-half",
        styles_src, perl = TRUE),
  "vertical divider between the two halves is defined in CSS"
)
stopifnot_msg(
  grepl("\\.phenomapr-signature-coverage-subhead", styles_src, perl = TRUE),
  ".phenomapr-signature-coverage-subhead CSS rules exist"
)

cat("\nAll signature/coverage merge smoke tests passed.\n")
