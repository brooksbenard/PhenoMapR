#!/usr/bin/env Rscript
## Quick offline smoke tests for the latest Score-tab UI changes:
##   1. plot_score_distribution() output still accepts +geom_vline() / +annotate()
##      overlays (the new "dashed tail thresholds" feature in app.R).
##   2. The new three-row "Scoring summary" UI block renders the expected
##      sds-row / sds-label / sds-value classes for the (data, signature,
##      cancer-type) triple under a representative set of input states.

suppressPackageStartupMessages({
  library(ggplot2)
  library(shiny)
})

stopifnot(requireNamespace("PhenoMapR", quietly = TRUE))

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

## --- 1. plot_score_distribution + dashed quantile overlay ----------------
set.seed(1)
v   <- as.numeric(scale(rnorm(500)))
df  <- data.frame(score = v)
ttl <- "PhenoMapR score distribution (smoke)"
p   <- PhenoMapR::plot_score_distribution(df, score_column = "score",
                                          main = ttl, base_size = 13)

pct  <- 0.05
q_lo <- stats::quantile(v, pct,       na.rm = TRUE, names = FALSE)
q_hi <- stats::quantile(v, 1 - pct,   na.rm = TRUE, names = FALSE)

p2 <- p +
  ggplot2::geom_vline(xintercept = c(q_lo, q_hi),
                      linetype = "dashed", colour = "#264653",
                      linewidth = 0.7) +
  ggplot2::annotate("text", x = q_lo, y = Inf, label = "Tail = 5%",
                    hjust = 1.05, vjust = 1.3, size = 3.2, colour = "#264653") +
  ggplot2::annotate("text", x = q_hi, y = Inf, label = "Tail = 5%",
                    hjust = -0.05, vjust = 1.3, size = 3.2, colour = "#264653")

stopifnot_msg(inherits(p2, "ggplot"),
              "+geom_vline / +annotate yields a ggplot")

## Build the plot so geoms actually get instantiated and confirm a vline layer
## made it into the rendered output.
gb <- suppressWarnings(ggplot2::ggplot_build(p2))
stopifnot_msg(length(gb$data) >= 1,
              "ggplot_build succeeds on the overlaid histogram")

geom_classes <- vapply(p2$layers,
                       function(L) class(L$geom)[1L], character(1))
stopifnot_msg("GeomVline" %in% geom_classes,
              "histogram has a GeomVline layer for tail cutoffs")
stopifnot_msg("GeomText" %in% geom_classes,
              "histogram has a GeomText layer for tail labels")

## --- 2. Re-create the score_data_status renderUI logic (mirrors app.R) ----

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

make_status <- function(expr_summary, reference_choice, custom_source = NULL,
                        cancer_type = NULL) {
  s  <- expr_summary
  rc <- reference_choice %||% "(none)"

  if (is.null(s)) {
    data_row <- tags$div(
      class = "sds-row sds-row-empty",
      tags$div(class = "sds-icon", icon("circle-info")),
      tags$div(class = "sds-content",
               tags$div(class = "sds-label", "Input data"),
               tags$em("Load an expression dataset in 1. Data first."))
    )
  } else {
    kind_label <- switch(s$kind %||% "loaded",
      seurat  = "Seurat object", sce = "SingleCellExperiment",
      spatial = "Spatial dataset", matrix = "Expression matrix",
      anndata = "AnnData", sprintf("%s object", s$kind %||% "loaded"))
    bits <- character(0)
    if (!is.na(s$default_assay %||% NA_character_) && nzchar(s$default_assay)) {
      bits <- c(bits, sprintf("assay: %s", s$default_assay))
    }
    if (length(s$layers_avail %||% character(0))) {
      bits <- c(bits, sprintf("layers: %s",
                              paste(s$layers_avail, collapse = ", ")))
    }
    data_row <- tags$div(
      class = "sds-row sds-row-ok",
      tags$div(class = "sds-icon", icon("check-circle")),
      tags$div(class = "sds-content",
               tags$div(class = "sds-label", "Input data"),
               tags$div(class = "sds-value", kind_label),
               if (length(bits))
                 tags$div(class = "sds-detail",
                          paste(bits, collapse = "  \u00B7  ")))
    )
  }

  phen_label <- switch(rc,
    "precog"           = "PRECOG (built-in)",
    "tcga"             = "TCGA (built-in)",
    "pediatric_precog" = "Pediatric PRECOG (built-in)",
    "ici_precog"       = "ICI PRECOG (built-in)",
    "_custom"          = {
      cs <- custom_source %||% "upload"
      if (identical(cs, "derive"))
        "Custom signature (derived from your bulk + phenotype)"
      else
        "Custom signature (uploaded)"
    },
    sprintf("%s", rc))
  phen_icon <- if (identical(rc, "_custom")) "flask" else "dna"
  phen_row <- tags$div(class = "sds-row sds-row-ok",
    tags$div(class = "sds-icon", icon(phen_icon)),
    tags$div(class = "sds-content",
             tags$div(class = "sds-label", "Phenotype source"),
             tags$div(class = "sds-value", phen_label)))

  if (identical(rc, "_custom")) {
    ct_row <- tags$div(class = "sds-row sds-row-na",
      tags$div(class = "sds-icon", icon("minus")),
      tags$div(class = "sds-content",
               tags$div(class = "sds-label", "Cancer / tissue type"),
               tags$em("Not applicable for custom signatures")))
  } else {
    ct <- cancer_type
    if (is.null(ct) || !nzchar(ct)) {
      ct_row <- tags$div(class = "sds-row sds-row-empty",
        tags$div(class = "sds-icon", icon("circle-question")),
        tags$div(class = "sds-content",
                 tags$div(class = "sds-label", "Cancer / tissue type"),
                 tags$em("Select one in 2. Phenotype.")))
    } else {
      ct_row <- tags$div(class = "sds-row sds-row-ok",
        tags$div(class = "sds-icon", icon("ribbon")),
        tags$div(class = "sds-content",
                 tags$div(class = "sds-label", "Cancer / tissue type"),
                 tags$div(class = "sds-value", ct)))
    }
  }

  tags$div(class = "score-data-status", data_row, phen_row, ct_row)
}

html_of <- function(x) as.character(x)

## ---- Case A: no data, defaults --------------------------------------------
h <- html_of(make_status(expr_summary = NULL,
                         reference_choice = "precog",
                         cancer_type = "Breast cancer"))
stopifnot_msg(grepl("score-data-status", h),
              "case A renders score-data-status wrapper")
stopifnot_msg(grepl("Load an expression dataset", h),
              "case A surfaces the no-data prompt")
stopifnot_msg(grepl("PRECOG \\(built-in\\)", h),
              "case A renders the PRECOG signature label")
stopifnot_msg(grepl("Breast cancer", h),
              "case A renders the cancer type")

## ---- Case B: Seurat object + TCGA + cancer type ---------------------------
expr_seu <- list(kind = "seurat", default_assay = "RNA",
                 layers_avail = c("data", "counts"))
h <- html_of(make_status(expr_summary = expr_seu,
                         reference_choice = "tcga",
                         cancer_type = "Lung adenocarcinoma"))
stopifnot_msg(grepl("Seurat object", h),
              "case B labels the input as Seurat")
stopifnot_msg(grepl("assay: RNA", h, fixed = TRUE),
              "case B surfaces the detected assay")
stopifnot_msg(grepl("layers: data, counts", h, fixed = TRUE),
              "case B surfaces the available layers")
stopifnot_msg(grepl("TCGA \\(built-in\\)", h),
              "case B renders the TCGA signature label")
stopifnot_msg(grepl("Lung adenocarcinoma", h),
              "case B renders the cancer type")

## ---- Case C: Matrix + custom (derive) signature ---------------------------
expr_mat <- list(kind = "matrix")
h <- html_of(make_status(expr_summary = expr_mat,
                         reference_choice = "_custom",
                         custom_source = "derive"))
stopifnot_msg(grepl("Expression matrix", h),
              "case C labels the input as a matrix")
stopifnot_msg(grepl("Custom signature \\(derived", h),
              "case C labels the signature as custom (derived)")
stopifnot_msg(grepl("Not applicable for custom signatures", h),
              "case C marks cancer-type as N/A for custom signatures")

## ---- Case D: AnnData + Pediatric PRECOG + no cancer-type yet --------------
expr_ad <- list(kind = "anndata", default_assay = NA_character_,
                layers_avail = character(0))
h <- html_of(make_status(expr_summary = expr_ad,
                         reference_choice = "pediatric_precog",
                         cancer_type = ""))
stopifnot_msg(grepl("AnnData", h),
              "case D labels the input as AnnData")
stopifnot_msg(grepl("Pediatric PRECOG \\(built-in\\)", h),
              "case D labels the signature as Pediatric PRECOG")
stopifnot_msg(grepl("Select one in 2\\. Phenotype", h),
              "case D prompts the user to pick a cancer type")

cat("\nAll Score-tab UI smoke tests passed.\n")
