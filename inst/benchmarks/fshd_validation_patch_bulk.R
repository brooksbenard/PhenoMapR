#!/usr/bin/env Rscript
## Rebuild bulk validation cache (GSE115650, GSE56787, GSE26852) with fixed
## expression/metadata alignment. Uses vignettes/.geo_cache when present.
## Run from repo root:
##   Rscript inst/benchmarks/fshd_validation_patch_bulk.R

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
repo_root <- if (length(script_path) == 1L && nzchar(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
} else {
  getwd()
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) repo_root <- getwd()
out_dir <- file.path(repo_root, "vignettes")
geo_cache <- file.path(out_dir, ".geo_cache")
val_base <- if (dir.exists(geo_cache)) geo_cache else file.path(tempdir(), "fshd_val_cache")
dir.create(val_base, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  } else {
    library(PhenoMapR)
  }
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GEOquery)
  library(Biobase)
})

fpkm_to_tpm <- function(fpkm) {
  fpkm <- as.matrix(fpkm)
  cs <- colSums(fpkm, na.rm = TRUE)
  cs[cs <= 0] <- NA_real_
  sweep(fpkm, 2L, cs, "/") * 1e6
}

collapse_ensembl_tpm_to_symbol <- function(tpm_ens_rows_samples_cols) {
  ens <- rownames(tpm_ens_rows_samples_cols)
  ens_base <- sub("\\.[0-9]+$", "", ens)
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = ens_base, column = "SYMBOL",
    keytype = "ENSEMBL", multiVals = "first"
  )
  df <- as.data.frame(tpm_ens_rows_samples_cols)
  df$gene_symbol <- sym
  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), , drop = FALSE]
  sample_cols <- setdiff(names(df), "gene_symbol")
  agg <- stats::aggregate(
    df[, sample_cols, drop = FALSE],
    by = list(gene_symbol = df$gene_symbol), FUN = mean, na.rm = TRUE
  )
  mat <- as.matrix(agg[, -1, drop = FALSE])
  rownames(mat) <- agg$gene_symbol
  mat
}

collapse_ensembl_counts_to_symbol_sum <- function(cnt) {
  ens <- rownames(cnt)
  if (any(grepl("^ENSG", ens))) {
    ens_base <- sub("\\.[0-9]+$", "", ens)
    sym <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = ens_base, column = "SYMBOL",
      keytype = "ENSEMBL", multiVals = "first"
    )
  } else {
    sym <- ens
  }
  df <- as.data.frame(cnt)
  df$gene_symbol <- sym
  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), , drop = FALSE]
  sample_cols <- setdiff(names(df), "gene_symbol")
  agg <- stats::aggregate(
    df[, sample_cols, drop = FALSE],
    by = list(gene_symbol = df$gene_symbol), FUN = sum, na.rm = TRUE
  )
  mat <- as.matrix(agg[, -1, drop = FALSE])
  rownames(mat) <- agg$gene_symbol
  mat
}

pick_first_existing_file <- function(dir, patterns) {
  for (pat in patterns) {
    hits <- list.files(dir, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

infer_sample_key <- function(pd) {
  if ("geo_accession" %in% names(pd)) return(as.character(pd$geo_accession))
  rownames(pd)
}

match_pd_to_expr <- function(expr_colnames, pd) {
  key <- infer_sample_key(pd)
  m <- match(expr_colnames, key)
  if (!anyNA(m)) return(m)
  if ("title" %in% names(pd)) {
    m2 <- match(expr_colnames, as.character(pd$title))
    if (!anyNA(m2)) return(m2)
  }
  rep(NA_integer_, length(expr_colnames))
}

align_pd_to_expr <- function(expr, pd) {
  m <- match_pd_to_expr(colnames(expr), pd)
  if (anyNA(m)) {
    bad <- colnames(expr)[is.na(m)]
    stop(
      "Expression columns not found in GEO metadata: ",
      paste(head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) " ..." else ""
    )
  }
  list(expr = expr, pd = pd[m, , drop = FALSE])
}

collapse_characteristics <- function(pd) {
  ch_cols <- grep("^characteristics", names(pd), value = TRUE)
  if (length(ch_cols) == 0) return(rep("", nrow(pd)))
  apply(pd[, ch_cols, drop = FALSE], 1, function(x) paste(x, collapse = " | "))
}

score_bulk_matrix <- function(expr, ref_fshd) {
  val_score_col <- paste0("weighted_sum_score_", colnames(ref_fshd)[1])
  sc <- PhenoMap(
    expression = expr,
    reference = ref_fshd,
    z_score_cutoff = 2,
    verbose = FALSE
  )
  sc[, val_score_col, drop = TRUE]
}

message("Using GEO cache: ", val_base)

pheno_bulk <- readRDS(file.path(out_dir, "FSHD_GSE140261_phenotype.rds"))
tpm_sym <- readRDS(file.path(out_dir, "FSHD_GSE140261_bulk_expression.rds"))
ref_fshd <- derive_reference_from_bulk(
  bulk_expression = tpm_sym,
  phenotype = pheno_bulk,
  sample_id_column = "sample_id",
  phenotype_column = "group",
  phenotype_type = "binary",
  normalize = TRUE,
  binary_positive_reference = "first",
  verbose = FALSE
)

validation_plot_df <- NULL
validation_wilcox <- list()
append_bulk <- function(cohort, group, score) {
  validation_plot_df <<- rbind(
    validation_plot_df,
    data.frame(cohort = cohort, group = group, score = score, stringsAsFactors = FALSE)
  )
}

message("GSE115650 ...")
tryCatch({
  if (!dir.exists(file.path(val_base, "GSE115650"))) {
    suppressMessages(GEOquery::getGEOSuppFiles("GSE115650", makeDirectory = TRUE, baseDir = val_base))
  }
  fp <- pick_first_existing_file(
    file.path(val_base, "GSE115650"),
    c("GSE115650.*FPKM.*csv\\.gz$", "FPKM.*csv\\.gz$")
  )
  fpkm <- read.csv(gzfile(fp), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  tpm <- collapse_ensembl_tpm_to_symbol(fpkm_to_tpm(fpkm))
  g <- suppressMessages(GEOquery::getGEO("GSE115650", GSEMatrix = TRUE)[[1]])
  pd <- Biobase::pData(g)
  meta <- tolower(collapse_characteristics(pd))
  pd$group <- NA_character_
  pd$group[grepl("control|normal|healthy", meta)] <- "Control"
  pd$group[grepl("fshd|facioscapulohumeral", meta)] <- "FSHD"
  aligned <- align_pd_to_expr(tpm, pd)
  pd <- aligned$pd
  pd$group <- factor(pd$group, levels = c("Control", "FSHD"))
  x <- score_bulk_matrix(aligned$expr, ref_fshd)
  ok <- stats::complete.cases(x, pd$group) & !is.na(pd$group)
  if (sum(ok) > 0 && length(unique(pd$group[ok])) >= 2) {
    validation_wilcox[["GSE115650"]] <- stats::wilcox.test(x[ok] ~ pd$group[ok])
  }
  append_bulk("GSE115650", pd$group, x)
}, error = function(e) message("GSE115650 failed: ", conditionMessage(e)))

message("GSE56787 ...")
tryCatch({
  if (!dir.exists(file.path(val_base, "GSE56787"))) {
    suppressMessages(GEOquery::getGEOSuppFiles("GSE56787", makeDirectory = TRUE, baseDir = val_base))
  }
  fp <- pick_first_existing_file(
    file.path(val_base, "GSE56787"),
    c("GSE56787.*biopsy.*gene\\.counts.*csv\\.gz$", "biopsy.*gene\\.counts.*csv\\.gz$")
  )
  cnt <- read.csv(gzfile(fp), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  cpm <- sweep(collapse_ensembl_counts_to_symbol_sum(cnt), 2, colSums(cnt), "/") * 1e6
  g <- suppressMessages(GEOquery::getGEO("GSE56787", GSEMatrix = TRUE)[[1]])
  pd <- Biobase::pData(g)
  biopsy <- grepl("muscle biopsy", tolower(collapse_characteristics(pd)))
  pd <- pd[biopsy, , drop = FALSE]
  pd$group <- ifelse(grepl("fshd", tolower(as.character(pd$source_name_ch1))), "FSHD", "Control")
  pd$group <- factor(pd$group, levels = c("Control", "FSHD"))
  aligned <- align_pd_to_expr(cpm, pd)
  pd <- aligned$pd
  x <- score_bulk_matrix(aligned$expr, ref_fshd)
  ok <- stats::complete.cases(x, pd$group)
  if (sum(ok) > 0 && length(unique(pd$group[ok])) >= 2) {
    validation_wilcox[["GSE56787"]] <- stats::wilcox.test(x[ok] ~ pd$group[ok])
  }
  append_bulk("GSE56787", pd$group, x)
}, error = function(e) message("GSE56787 failed: ", conditionMessage(e)))

message("GSE26852 ...")
tryCatch({
  g <- suppressMessages(GEOquery::getGEO("GSE26852", GSEMatrix = TRUE)[[1]])
  pd <- Biobase::pData(g)
  ds <- pd$characteristics_ch1.3
  keep <- grepl("Normal Muscle", ds) | grepl("FSHD T2-STIR", ds)
  pd <- pd[keep, , drop = FALSE]
  pd$group <- ifelse(grepl("FSHD", pd$characteristics_ch1.3), "FSHD", "Control")
  pd$group <- factor(pd$group, levels = c("Control", "FSHD"))
  fd <- Biobase::fData(g)
  e <- Biobase::exprs(g)[, rownames(pd), drop = FALSE]
  sym <- fd$Symbol
  ok <- !is.na(sym) & nzchar(sym)
  e <- e[ok, , drop = FALSE]
  sym <- sym[ok]
  dfl <- data.frame(sym = sym, e, check.names = FALSE)
  agg <- stats::aggregate(. ~ sym, data = dfl, mean)
  rownames(agg) <- agg$sym
  mat <- as.matrix(agg[, -1, drop = FALSE])
  mat <- apply(mat, 2, function(v) v / sum(v) * 1e6)
  x <- score_bulk_matrix(mat, ref_fshd)
  ok_w <- stats::complete.cases(x, pd$group)
  if (sum(ok_w) > 0 && length(unique(pd$group[ok_w])) >= 2) {
    validation_wilcox[["GSE26852"]] <- stats::wilcox.test(x[ok_w] ~ pd$group[ok_w])
  }
  append_bulk("GSE26852", pd$group, x)
}, error = function(e) message("GSE26852 failed: ", conditionMessage(e)))

if (is.null(validation_plot_df) || nrow(validation_plot_df) == 0) {
  stop("Bulk validation patch produced no rows.")
}

validation_plot_df$score_scaled <- ave(
  validation_plot_df$score,
  validation_plot_df$cohort,
  FUN = function(x) as.numeric(scale(x))
)
bulk_path <- file.path(out_dir, "FSHD_validation_bulk_plot_df.rds")
saveRDS(validation_plot_df, bulk_path)
message("Saved ", bulk_path, " (", nrow(validation_plot_df), " rows)")

fmt_p <- function(w) {
  if (is.null(w)) return(NA_character_)
  format.pval(w$p.value, digits = 3, eps = 0.001)
}
summary_path <- file.path(out_dir, "FSHD_validation_summary_table.rds")
if (file.exists(summary_path)) {
  sm <- readRDS(summary_path)
  for (geo in names(validation_wilcox)) {
    idx <- which(sm$GEO == geo)
    if (length(idx) == 1L) sm$Wilcoxon_p[idx] <- fmt_p(validation_wilcox[[geo]])
  }
  saveRDS(sm, summary_path)
  message("Updated Wilcoxon p-values in summary table.")
}

for (c in unique(validation_plot_df$cohort)) {
  d <- validation_plot_df[validation_plot_df$cohort == c, ]
  cat(c, ": p =", stats::wilcox.test(d$score ~ d$group)$p.value,
      " FSHD>Control =", mean(d$score[d$group == "FSHD"]) > mean(d$score[d$group == "Control"]), "\n")
}

saveRDS(validation_wilcox, file.path(out_dir, "FSHD_validation_wilcox.rds"))
message("Done.")
