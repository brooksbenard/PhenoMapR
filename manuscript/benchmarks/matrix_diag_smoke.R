#!/usr/bin/env Rscript
## Smoke tests for the 1. Data tab matrix-diagnostics block.
##
## Re-creates the renderUI() logic from inst/shiny/app.R that builds the
## "Matrix diagnostics" panel (gene-ID style, expression format, sc vs
## bulk, optional cleanup block) so we can validate the rendered HTML
## across the cases the user is most likely to hit:
##   1. Raw bulk counts with HUGO genes -> cleanup block (log2(CPM+1)).
##   2. Single-cell raw counts (sparse) -> Seurat-style cleanup default.
##   3. Already log-normalized matrix -> "matrix looks ready" pill.
##   4. ENSG-prefixed matrix              -> HUGO cleanup highlighted.

suppressPackageStartupMessages({
  library(shiny)
})

stopifnot(requireNamespace("PhenoMapR", quietly = TRUE))
stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

## Pull the live detector from the dev tree so we don't accidentally call
## an outdated installed copy.
devtools::load_all(quiet = TRUE)

## Re-implement the renderUI logic verbatim (kept in lock-step with
## inst/shiny/app.R::output$expr_matrix_diagnostics). When the
## production renderUI changes, update this helper to match.
render_diag <- function(expr) {
  d <- PhenoMapR::detect_expression_format(expr, verbose = FALSE)

  fmt_class <- switch(d$format,
    raw_counts = "diag-warn", cpm_or_tpm = "diag-warn",
    z_scaled = "diag-warn",   log_normalized = "diag-ok",
    "diag-info")
  id_class <- switch(d$gene_id_kind,
    hugo = "diag-ok",     ensembl = "diag-warn",
    mixed = "diag-warn",  "diag-info")
  id_label <- switch(d$gene_id_kind,
    hugo    = sprintf("HUGO symbols (%d/%d)", d$n_hugo_like, d$n_genes),
    ensembl = sprintf("Ensembl IDs (%d/%d ENSG-prefixed)",
                      d$n_ensembl, d$n_genes),
    mixed   = sprintf("Mixed (%d ENSG + %d HUGO-like)",
                      d$n_ensembl, d$n_hugo_like),
    "Unknown gene-ID style")

  diag_row <- function(class_, label, value, detail = NULL) {
    tags$div(class = paste("diag-row", class_),
      tags$div(class = "diag-label", label),
      tags$div(class = "diag-value", value),
      if (!is.null(detail))
        tags$div(class = "diag-detail", detail))
  }

  needs_cleanup <- d$gene_id_kind %in% c("ensembl", "mixed") ||
    d$n_dup > 0L ||
    d$format %in% c("raw_counts", "cpm_or_tpm")
  default_mode <- if (identical(d$sc_or_bulk, "single_cell"))
    "single_cell" else "bulk"

  cleanup_block <- if (needs_cleanup) {
    tags$div(class = "diag-cleanup",
      tags$div(class = "diag-cleanup-title",
               tags$strong(" Cleanup options")),
      checkboxInput("diag_do_hugo",
        "Clean gene IDs to approved HUGO symbols (HGNChelper)",
        value = d$gene_id_kind %in% c("ensembl", "mixed") || d$n_dup > 0L),
      checkboxInput("diag_do_collapse",
        "Collapse duplicate gene rows by mean", value = d$n_dup > 0L),
      checkboxInput("diag_do_normalize",
        "Log-normalize expression",
        value = d$format %in% c("raw_counts", "cpm_or_tpm")),
      radioButtons("diag_norm_mode", "Normalization mode",
        choices = c("Auto-detect (recommended)" = "auto",
                    "Single cell (Seurat log1p of library-size-scaled counts)" = "single_cell",
                    "Bulk (log2(CPM + 1))" = "bulk"),
        selected = if (identical(d$sc_or_bulk, "unclear")) "auto" else default_mode,
        inline = FALSE),
      actionButton("diag_run_clean", "Clean & normalize",
                   class = "btn-primary btn-sm", width = "100%"))
  } else {
    tags$div(class = "diag-cleanup diag-cleanup-clean",
             tags$span(" Matrix looks ready -- no cleanup recommended."))
  }

  tags$div(class = "expr-matrix-diag",
    tags$div(class = "diag-title", tags$strong(" Matrix diagnostics")),
    diag_row(id_class, "Gene IDs", id_label,
             detail = if (d$n_dup > 0L)
               sprintf("%d duplicate gene ID(s)", d$n_dup) else NULL),
    diag_row(fmt_class, "Expression format",
             sprintf("%s (%s confidence)",
                     d$format_label, d$format_confidence)),
    diag_row("diag-info", "Data type", d$sc_or_bulk_label),
    cleanup_block,
    if (length(d$recommendations))
      tags$ul(class = "diag-recs",
              lapply(d$recommendations, tags$li)))
}

html_of <- function(x) as.character(x)

## ---- Case 1: bulk raw counts, HUGO genes ---------------------------------
set.seed(1)
bulk <- matrix(rpois(200 * 30, 5), nrow = 200,
               dimnames = list(c("TP53","MYC","EGFR","BRCA1","CDKN2A",
                                 paste0("GENE", 1:195)),
                               paste0("S", 1:30)))
h <- html_of(render_diag(bulk))
stopifnot_msg(grepl("expr-matrix-diag", h),
              "case 1 renders the diagnostics wrapper")
stopifnot_msg(grepl("Raw counts", h),
              "case 1 reports raw counts")
stopifnot_msg(grepl("HUGO symbols", h),
              "case 1 reports HUGO gene IDs")
stopifnot_msg(grepl("Bulk-like", h),
              "case 1 reports bulk-like data type")
stopifnot_msg(grepl("diag_run_clean", h, fixed = TRUE),
              "case 1 shows the Clean & normalize button")
stopifnot_msg(grepl('value="bulk".*checked|checked.*value="bulk"', h),
              "case 1 selects bulk normalization mode by default")

## ---- Case 2: single-cell raw counts (sparse) -----------------------------
set.seed(7)
sc <- matrix(0, nrow = 200, ncol = 400,
             dimnames = list(c("TP53","MYC", paste0("Gene", 1:198)),
                             paste0("Cell_", 1:400)))
nz <- sample(length(sc), 0.2 * length(sc))
sc[nz] <- rpois(length(nz), 2)
h <- html_of(render_diag(sc))
stopifnot_msg(grepl("Single-cell-like", h),
              "case 2 reports single-cell-like data type")
stopifnot_msg(grepl('value="single_cell".*checked|checked.*value="single_cell"', h),
              "case 2 selects single_cell normalization mode by default")

## ---- Case 3: already log-normalized -> no cleanup ------------------------
set.seed(3)
logn <- matrix(runif(200 * 30, 0, 8), nrow = 200,
               dimnames = list(paste0("G", 1:200), paste0("S", 1:30)))
h <- html_of(render_diag(logn))
stopifnot_msg(grepl("diag-cleanup-clean", h, fixed = TRUE),
              "case 3 hides cleanup and shows the 'ready' chip")
stopifnot_msg(!grepl("diag_run_clean", h, fixed = TRUE),
              "case 3 does not render the Clean & normalize button")

## ---- Case 4: ENSG-prefixed rownames + counts -----------------------------
set.seed(8)
ensg <- matrix(rpois(200 * 30, 5), nrow = 200,
               dimnames = list(sprintf("ENSG%011d", 1:200),
                               paste0("S", 1:30)))
h <- html_of(render_diag(ensg))
stopifnot_msg(grepl("Ensembl IDs", h),
              "case 4 reports Ensembl gene IDs")
stopifnot_msg(grepl('id="diag_do_hugo".*checked|checked.*id="diag_do_hugo"', h),
              "case 4 enables the HUGO cleanup checkbox by default")
stopifnot_msg(grepl("Gene IDs look like Ensembl", h),
              "case 4 surfaces the Ensembl -> HUGO recommendation")

cat("\nAll matrix-diagnostics renderUI smoke tests passed.\n")
