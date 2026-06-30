#!/usr/bin/env Rscript
## Build cached FSHD validation artifacts for the vignette.
## Run from repo root:
##   Rscript inst/benchmarks/fshd_validation_build_cache.R
##
## Writes vignettes/FSHD_validation_*.rds

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_root <- Sys.getenv("PHENOMAPR_REPO", unset = "")
if (!nzchar(repo_root)) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(script_path) == 1L && nzchar(script_path)) {
    repo_root <- normalizePath(
      file.path(dirname(script_path), "..", ".."),
      mustWork = FALSE
    )
  }
}
if (!nzchar(repo_root) || !file.exists(file.path(repo_root, "DESCRIPTION"))) {
  repo_root <- getwd()
}
out_dir <- file.path(repo_root, "vignettes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  } else {
    library(PhenoMapR)
  }
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery: BiocManager::install('GEOquery')")
}

val_base <- file.path(out_dir, ".geo_cache")
if (!dir.exists(val_base)) {
  val_base <- file.path(tempdir(), "fshd_val_cache")
}
dir.create(val_base, showWarnings = FALSE, recursive = TRUE)

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

parse_ensembl_symbol_rownames <- function(rn) {
  sym <- sub("^[^-]+-", "", rn)
  sym[nzchar(sym)] <- sym[nzchar(sym)]
  sym
}

collapse_symbol_matrix <- function(mat, row_symbols, fun = sum) {
  df <- as.data.frame(mat)
  df$gene_symbol <- row_symbols
  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), , drop = FALSE]
  sample_cols <- setdiff(names(df), "gene_symbol")
  agg <- stats::aggregate(
    df[, sample_cols, drop = FALSE],
    by = list(gene_symbol = df$gene_symbol), FUN = fun, na.rm = TRUE
  )
  out <- as.matrix(agg[, -1, drop = FALSE])
  rownames(out) <- agg$gene_symbol
  out
}

pick_first_existing_file <- function(dir, patterns) {
  for (pat in patterns) {
    hits <- list.files(dir, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

load_wellstone_comprehensive_df <- function(cache_dir) {
  local_path <- file.path(cache_dir, "FSHD_wellstone_comprehensive_df.rda")
  url <- paste0(
    "https://github.com/fredhutch/Wellstone_BiLateral_Biopsy/raw/master/data/",
    "comprehensive_df.rda"
  )
  if (!file.exists(local_path)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    utils::download.file(url, local_path, mode = "wb", quiet = TRUE)
  }
  env <- new.env()
  load(local_path, envir = env)
  env$comprehensive_df
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

## Reorder pData rows to match expression column order; keep expression columns as-is.
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

extract_char_field <- function(pd, field) {
  ch <- collapse_characteristics(pd)
  pat <- paste0("(?i)", field, ":\\s*([^|]+)")
  out <- sub(paste0(".*", pat, ".*"), "\\1", ch, perl = TRUE)
  out[!grepl(pat, ch, perl = TRUE)] <- NA_character_
  trimws(out)
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

derive_fshd_induced_genes <- function(tpm_gse143453) {
  fshd_cols <- grep("FSHD", colnames(tpm_gse143453), value = TRUE)
  ctrl_cols <- grep("Ctrl|Control", colnames(tpm_gse143453), value = TRUE)
  fshd_late <- fshd_cols[grepl("Day_[2-9]|Day[2-9]", fshd_cols)]
  ctrl_late <- ctrl_cols[grepl("Day_[2-9]|Day[2-9]", ctrl_cols)]
  if (length(fshd_late) == 0L || length(ctrl_late) == 0L) {
    stop("Could not identify late differentiation columns in GSE143453 TPM")
  }
  mu_fshd <- rowMeans(tpm_gse143453[, fshd_late, drop = FALSE], na.rm = TRUE)
  mu_ctrl <- rowMeans(tpm_gse143453[, ctrl_late, drop = FALSE], na.rm = TRUE)
  lfc <- log2(pmax(mu_fshd, 0) + 1) - log2(pmax(mu_ctrl, 0) + 1)
  top54 <- names(sort(lfc, decreasing = TRUE))[seq_len(54L)]
  top54 <- unique(c("DUX4", top54))
  top54[seq_len(min(55L, length(top54)))]
}

parse_snrna_cell_meta <- function(cells) {
  day <- ifelse(grepl("Day5|Day_5", cells), "Day5",
                ifelse(grepl("Day3|Day_3", cells), "Day3", "Other"))
  status <- ifelse(grepl("FSHD", cells), "FSHD2", "Control")
  lib <- sub("\\.cell.*$", "", cells)
  replicate <- sub(".*_(R[0-9]+)$", "\\1", lib)
  replicate[!grepl("^R[0-9]+$", replicate)] <- "R1"
  data.frame(
    cell = cells,
    day = day,
    status = status,
    library = lib,
    replicate = replicate,
    stringsAsFactors = FALSE
  )
}

classify_fshd_hi_lo <- function(expr_sym, fshd_genes) {
  genes <- intersect(fshd_genes, rownames(expr_sym))
  if (length(genes) == 0L || ncol(expr_sym) < 4L) {
    return(rep(NA_character_, ncol(expr_sym)))
  }
  enrichment <- colMeans(log1p(expr_sym[genes, , drop = FALSE]), na.rm = TRUE)
  km <- stats::kmeans(enrichment, centers = 2, nstart = 25L)
  hi_cluster <- which.max(vapply(split(enrichment, km$cluster), mean, numeric(1)))
  ifelse(km$cluster == hi_cluster, "FSHD-Hi", "FSHD-Lo")
}

save_checkpoint <- function() {
  if (!is.null(validation_plot_df) && nrow(validation_plot_df) > 0) {
    validation_plot_df$score_scaled <- ave(
      validation_plot_df$score,
      validation_plot_df$cohort,
      FUN = function(x) as.numeric(scale(x))
    )
  }
  saveRDS(validation_plot_df, file.path(out_dir, "FSHD_validation_bulk_plot_df.rds"))
  saveRDS(validation_wilcox, file.path(out_dir, "FSHD_validation_wilcox.rds"))
  saveRDS(anatomical_df, file.path(out_dir, "FSHD_validation_anatomical_df.rds"))
  saveRDS(longitudinal_df, file.path(out_dir, "FSHD_validation_longitudinal_df.rds"))
  saveRDS(gse242912_df, file.path(out_dir, "FSHD_validation_gse242912_df.rds"))
  saveRDS(snrna_df, file.path(out_dir, "FSHD_validation_snRNA_scores.rds"))
  if (exists("fshd_induced_genes")) {
    saveRDS(fshd_induced_genes, file.path(out_dir, "FSHD_validation_fshd_induced_genes.rds"))
  }
  if (exists("gse242912_lr_df")) {
    saveRDS(gse242912_lr_df, file.path(out_dir, "FSHD_validation_gse242912_lr_df.rds"))
  }
  if (exists("validation_summary_df")) {
    saveRDS(validation_summary_df, file.path(out_dir, "FSHD_validation_summary_table.rds"))
  }
  gc()
}

## Low-memory GSE143452 snRNA scoring (no Seurat object).
score_gse143452_snrna_lowmem <- function(
    fp,
    ref_fshd,
    fshd_induced_genes,
    days = c("Day3", "Day5"),
    max_cells_per_group = 2000L,
    seed = 1L
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table required for low-memory GSE143452 loading")
  }
  hdr <- names(data.table::fread(fp, nrows = 0L))
  day_pat <- paste(days, collapse = "|")
  target_cols <- grep(paste0("(", day_pat, ")", ".*"), hdr, value = TRUE)
  target_cols <- target_cols[grepl("FSHD|Ctrl|Control", target_cols)]
  fshd_cols <- target_cols[grepl("FSHD", target_cols)]
  ctrl_cols <- target_cols[grepl("Ctrl|Control", target_cols)]
  set.seed(seed)
  subsample_cols <- function(cols, cap) {
    if (length(cols) <= cap) return(cols)
    sample(cols, cap)
  }
  use_cols <- c(
    subsample_cols(fshd_cols, max_cells_per_group),
    subsample_cols(ctrl_cols, max_cells_per_group)
  )
  use_cols <- unique(use_cols)
  use_cols <- c(hdr[1L], use_cols)
  message(
    "GSE143452: reading ", length(use_cols) - 1L,
    " nuclei (", paste(days, collapse = ", "), "; column-selective)"
  )

  dt <- data.table::fread(fp, select = use_cols, data.table = FALSE)
  rn <- dt[[1L]]
  sym <- parse_ensembl_symbol_rownames(rn)
  mat <- as.matrix(dt[, -1L, drop = FALSE])
  rownames(mat) <- sym
  rm(dt)
  gc()

  mat <- collapse_symbol_matrix(mat, sym, fun = sum)
  cells <- colnames(mat)
  meta <- parse_snrna_cell_meta(cells)

  ref_syms <- rownames(ref_fshd)[abs(ref_fshd[, 1L]) >= 2]
  score_genes <- intersect(ref_syms, rownames(mat))
  score_mat <- mat[score_genes, , drop = FALSE]

  val_score_col <- paste0("weighted_sum_score_", colnames(ref_fshd)[1L])
  scores <- PhenoMap(
    expression = score_mat,
    reference = ref_fshd,
    z_score_cutoff = 2,
    verbose = FALSE
  )

  fshd_state <- rep(NA_character_, length(cells))
  names(fshd_state) <- cells
  for (dy in days) {
    fshd_cells <- cells[meta$day == dy & meta$status == "FSHD2"]
    if (length(fshd_cells) > 0L) {
      fshd_state[fshd_cells] <- classify_fshd_hi_lo(
        mat[, fshd_cells, drop = FALSE],
        fshd_induced_genes
      )
    }
  }
  induced_present <- intersect(fshd_induced_genes, rownames(mat))
  n_detected <- colSums(mat[induced_present, , drop = FALSE] > 0, na.rm = TRUE)

  out <- data.frame(
    meta,
    score = scores[cells, val_score_col, drop = TRUE],
    n_fshd_genes_detected = as.integer(n_detected[cells]),
    fshd_state = unname(fshd_state[cells]),
    stringsAsFactors = FALSE
  )
  out$score_scaled <- ave(out$score, out$day, FUN = function(x) as.numeric(scale(x)))
  rm(mat, score_mat)
  gc()
  out
}

message("Loading GSE140261 reference ...")
pheno_path <- file.path(out_dir, "FSHD_GSE140261_phenotype.rds")
tpm_path <- file.path(out_dir, "FSHD_GSE140261_bulk_expression.rds")
if (!file.exists(pheno_path) || !file.exists(tpm_path)) {
  if (requireNamespace("googledrive", quietly = TRUE)) {
    options(googledrive_quiet = TRUE)
    googledrive::drive_deauth()
    googledrive::drive_download(
      googledrive::as_id("13jTizdcmWPYMlEnjEIZ-V8U1wQzch2o1"),
      pheno_path, overwrite = TRUE
    )
    googledrive::drive_download(
      googledrive::as_id("14g_uHDVbjDA2Tl6-ttY7yF76s8Sy7CuH"),
      tpm_path, overwrite = TRUE
    )
  } else {
    stop("Missing GSE140261 RDS in vignettes/ and googledrive unavailable.")
  }
}
pheno_bulk <- readRDS(pheno_path)
tpm_sym <- readRDS(tpm_path)
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
val_score_col <- paste0("weighted_sum_score_", colnames(ref_fshd)[1])

validation_plot_df <- NULL
validation_wilcox <- list()
anatomical_df <- NULL
longitudinal_df <- NULL
gse242912_df <- NULL
snrna_df <- NULL

append_bulk <- function(cohort, group, score) {
  validation_plot_df <<- rbind(
    validation_plot_df,
    data.frame(cohort = cohort, group = group, score = score, stringsAsFactors = FALSE)
  )
}

## --- GSE115650 ---------------------------------------------------------------
message("GSE115650 ...")
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE115650", makeDirectory = TRUE, baseDir = val_base))
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
  tpm <- aligned$expr
  pd <- aligned$pd
  pd$group <- factor(pd$group, levels = c("Control", "FSHD"))
  x <- score_bulk_matrix(tpm, ref_fshd)
  names(x) <- infer_sample_key(pd)
  ok <- stats::complete.cases(x, pd$group) & !is.na(pd$group)
  if (sum(ok) > 0 && length(unique(pd$group[ok])) >= 2) {
    validation_wilcox[["GSE115650"]] <- stats::wilcox.test(x[ok] ~ pd$group[ok])
  }
  append_bulk("GSE115650", pd$group, x)
}, error = function(e) message("GSE115650 failed: ", conditionMessage(e)))

## --- GSE56787 ----------------------------------------------------------------
message("GSE56787 ...")
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE56787", makeDirectory = TRUE, baseDir = val_base))
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
  cpm <- aligned$expr
  pd <- aligned$pd
  x <- score_bulk_matrix(cpm, ref_fshd)
  ok <- stats::complete.cases(x, pd$group)
  if (sum(ok) > 0 && length(unique(pd$group[ok])) >= 2) {
    validation_wilcox[["GSE56787"]] <- stats::wilcox.test(x[ok] ~ pd$group[ok])
  }
  append_bulk("GSE56787", pd$group, x)
}, error = function(e) message("GSE56787 failed: ", conditionMessage(e)))

## --- GSE26852 ----------------------------------------------------------------
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

## --- GSE143453 (FSHD-induced gene list for snRNA validation) -----------------
message("GSE143453 ...")
fshd_induced_genes <- NULL
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE143453", makeDirectory = TRUE, baseDir = val_base))
  fp <- pick_first_existing_file(
    file.path(val_base, "GSE143453"),
    c("Nomalized_TPM_table\\.csv\\.gz$", "Normalized_TPM.*csv\\.gz$")
  )
  tpm_full <- read.csv(gzfile(fp), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  fshd_induced_genes <- derive_fshd_induced_genes(tpm_full)
  rm(tpm_full)
}, error = function(e) message("GSE143453 failed: ", conditionMessage(e)))

## --- GSE174301 (anatomical validation source data) -------------------------
message("GSE174301 ...")
gse174301_pd <- NULL
gse174301_tpm <- NULL
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE174301", makeDirectory = TRUE, baseDir = val_base))
  tar_path <- pick_first_existing_file(
    file.path(val_base, "GSE174301"),
    c("RAW\\.tar$")
  )
  td <- tempfile("gse174301")
  dir.create(td)
  untar(tar_path, exdir = td)
  rsem_files <- list.files(td, pattern = "rsem\\.genes\\.results\\.gz$", recursive = TRUE, full.names = TRUE)
  g <- suppressMessages(GEOquery::getGEO("GSE174301", GSEMatrix = TRUE)[[1]])
  pd <- Biobase::pData(g)
  pd$fshd_status <- extract_char_field(pd, "fshd status")
  pd$muscle_source <- extract_char_field(pd, "muscle source")
  pd$day <- extract_char_field(pd, "day of differentiation")
  pd$title <- as.character(pd$title)
  tpm_list <- list()
  for (i in seq_len(nrow(pd))) {
    gsm <- rownames(pd)[i]
    title <- pd$title[i]
    hit <- rsem_files[grepl(title, rsem_files, fixed = TRUE)]
    if (length(hit) == 0) hit <- rsem_files[grepl(gsm, rsem_files, fixed = TRUE)]
    if (length(hit) == 0) next
    res <- read.delim(gzfile(hit[1]), stringsAsFactors = FALSE)
    sym <- res$gene_id
    tpm_list[[gsm]] <- setNames(res$TPM, sym)
  }
  if (length(tpm_list) > 0) {
    all_genes <- unique(unlist(lapply(tpm_list, names)))
    mat <- matrix(0, nrow = length(all_genes), ncol = length(tpm_list),
                  dimnames = list(all_genes, names(tpm_list)))
    for (nm in names(tpm_list)) mat[names(tpm_list[[nm]]), nm] <- tpm_list[[nm]]
    mat <- collapse_ensembl_tpm_to_symbol(mat)
    keep <- intersect(rownames(pd), colnames(mat))
    pd <- pd[keep, , drop = FALSE]
    mat <- mat[, keep, drop = FALSE]
    gse174301_pd <- pd
    gse174301_tpm <- mat
  }
}, error = function(e) message("GSE174301 failed: ", conditionMessage(e)))

save_checkpoint()
message("Checkpoint: bulk validation cohorts saved.")

## --- Anatomical GSE174301 (FSHD, TA vs quad, all timepoints) ----------------
anatomical_df <- NULL
if (!is.null(gse174301_pd) && !is.null(gse174301_tpm)) {
  message("GSE174301 anatomical ...")
  fshd <- gse174301_pd[grepl("FSHD", gse174301_pd$fshd_status, ignore.case = TRUE), , drop = FALSE]
  fshd$muscle <- ifelse(grepl("tibialis|ta", fshd$muscle_source, ignore.case = TRUE), "TA", "Quad")
  fshd <- fshd[fshd$muscle %in% c("TA", "Quad"), , drop = FALSE]
  fshd$day_num <- suppressWarnings(as.numeric(sub(".*Day\\s*([0-9]+).*", "\\1", fshd$day)))
  fshd$day_lab <- ifelse(is.na(fshd$day_num), fshd$day, paste0("Day ", fshd$day_num))
  if (nrow(fshd) > 0) {
    x <- score_bulk_matrix(gse174301_tpm[, rownames(fshd), drop = FALSE], ref_fshd)
    anatomical_df <- data.frame(
      sample = rownames(fshd),
      muscle = fshd$muscle,
      day = fshd$day,
      day_lab = fshd$day_lab,
      day_num = fshd$day_num,
      fshd_status = fshd$fshd_status,
      score = x,
      stringsAsFactors = FALSE
    )
    anatomical_df$score_scaled <- as.numeric(scale(anatomical_df$score))
  }
}

## --- Longitudinal GSE115650 + GSE140261 --------------------------------------
message("Longitudinal ...")
longitudinal_df <- NULL
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE115650", makeDirectory = TRUE, baseDir = val_base))
  fp <- pick_first_existing_file(
    file.path(val_base, "GSE115650"),
    c("GSE115650.*FPKM.*csv\\.gz$", "FPKM.*csv\\.gz$")
  )
  fpkm <- read.csv(gzfile(fp), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  tpm115 <- collapse_ensembl_tpm_to_symbol(fpkm_to_tpm(fpkm))
  g115 <- suppressMessages(GEOquery::getGEO("GSE115650", GSEMatrix = TRUE)[[1]])
  pd115 <- Biobase::pData(g115)
  pd115$patient_id <- extract_char_field(pd115, "sample_id")
  m115 <- match_pd_to_expr(colnames(tpm115), pd115)
  ok115 <- !is.na(m115)
  pd115 <- pd115[m115[ok115], , drop = FALSE]
  tpm115 <- tpm115[, ok115, drop = FALSE]
  meta115 <- tolower(collapse_characteristics(pd115))
  pd115 <- pd115[grepl("fshd|facioscapulohumeral", meta115), , drop = FALSE]
  tpm115_fshd <- tpm115[, grepl("fshd|facioscapulohumeral", meta115), drop = FALSE]
  patients115 <- unique(pd115$patient_id)
  tpm_pat <- do.call(cbind, lapply(patients115, function(pid) {
    cols <- which(pd115$patient_id == pid)
    rowMeans(tpm115_fshd[, cols, drop = FALSE], na.rm = TRUE)
  }))
  rownames(tpm_pat) <- rownames(tpm115_fshd)
  colnames(tpm_pat) <- patients115
  score_v1 <- score_bulk_matrix(tpm_pat, ref_fshd)
  score_v1 <- setNames(as.numeric(score_v1), patients115)

  fshd_v2 <- pheno_bulk[pheno_bulk$group == "FSHD", , drop = FALSE]
  fshd_v2$patient_id <- sub("b$", "", fshd_v2$sample_id)
  score_v2 <- score_bulk_matrix(
    tpm_sym[, fshd_v2$sample_id, drop = FALSE],
    ref_fshd
  )
  names(score_v2) <- fshd_v2$patient_id

  g140 <- suppressMessages(GEOquery::getGEO("GSE140261", GSEMatrix = TRUE)[[1]])
  pd140 <- Biobase::pData(g140)
  pd140$title <- as.character(pd140$title)
  pd140$patient_id <- extract_char_field(pd140, "patient_id")
  pd140$dux4_rlogsum <- suppressWarnings(as.numeric(extract_char_field(pd140, "dux4.rlogsum")))
  pd140$fat_fraction <- suppressWarnings(as.numeric(extract_char_field(pd140, "fat_fraction")))
  pd140$stir_rating <- suppressWarnings(as.numeric(extract_char_field(pd140, "stir_rating")))

  common <- intersect(patients115, names(score_v2))
  if (length(common) > 0) {
    match_v2_num <- function(field) {
      suppressWarnings(as.numeric(extract_char_field(pd140, field)))[
        match(paste0(common, "b"), pd140$title)
      ]
    }
    mean_v1_num <- function(field) {
      vapply(common, function(pid) {
        idx <- which(pd115$patient_id == pid)
        if (length(idx) == 0L) return(NA_real_)
        mean(
          suppressWarnings(as.numeric(extract_char_field(pd115[idx, , drop = FALSE], field))),
          na.rm = TRUE
        )
      }, numeric(1))
    }
    longitudinal_df <- data.frame(
      patient_id = common,
      score_visit1 = score_v1[common],
      score_visit2 = score_v2[common],
      path_score_v1 = mean_v1_num("path.score"),
      inflam_v1 = mean_v1_num("inflam"),
      active_v1 = mean_v1_num("active"),
      stir_v1 = mean_v1_num("stir"),
      t1_v1 = mean_v1_num("t1"),
      fat_fraction_v1 = mean_v1_num("fat.fraction"),
      rnaseq_score_v1 = mean_v1_num("rnaseq.score"),
      path_score_v2 = match_v2_num("pathology.score"),
      inflam_v2 = match_v2_num("inflammation"),
      t1_rating_v2 = match_v2_num("t1_rating"),
      fat_fraction_v2 = match_v2_num("fat_fraction"),
      stir_rating_v2 = match_v2_num("stir_rating"),
      dux4_rlogsum_v2 = match_v2_num("dux4.rlogsum"),
      stringsAsFactors = FALSE
    )
    longitudinal_df$score_visit1_scaled <- as.numeric(scale(longitudinal_df$score_visit1))
    longitudinal_df$score_visit2_scaled <- as.numeric(scale(longitudinal_df$score_visit2))
    longitudinal_df$score_delta <- longitudinal_df$score_visit2 - longitudinal_df$score_visit1
    longitudinal_df$score_delta_scaled <- with(
      longitudinal_df,
      score_visit2_scaled - score_visit1_scaled
    )
    longitudinal_df$delta_fat_fraction <- with(
      longitudinal_df,
      fat_fraction_v2 - fat_fraction_v1
    )
    longitudinal_df$delta_stir <- with(longitudinal_df, stir_rating_v2 - stir_v1)
    longitudinal_df$delta_path_score <- with(
      longitudinal_df,
      path_score_v2 - path_score_v1
    )
    longitudinal_df$delta_inflam <- with(longitudinal_df, inflam_v2 - inflam_v1)
    longitudinal_df$delta_t1 <- with(longitudinal_df, t1_rating_v2 - t1_v1)
    longitudinal_df$score_increased <- longitudinal_df$score_delta > 0
  }
}, error = function(e) message("Longitudinal failed: ", conditionMessage(e)))

save_checkpoint()
message("Checkpoint: anatomical + longitudinal saved.")

## --- GSE242912 bilateral -----------------------------------------------------
message("GSE242912 ...")
gse242912_df <- NULL
gse242912_lr_df <- NULL
tryCatch({
  suppressMessages(GEOquery::getGEOSuppFiles("GSE242912", makeDirectory = TRUE, baseDir = val_base))
  fp_cnt <- pick_first_existing_file(
    file.path(val_base, "GSE242912"),
    c("bilat_gene_counts\\.csv\\.gz$")
  )
  fp_dds <- pick_first_existing_file(
    file.path(val_base, "GSE242912"),
    c("bilat_dds\\.rda\\.gz$")
  )
  cnt <- read.csv(gzfile(fp_cnt), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(cnt) <- cnt[[1]]
  cnt <- cnt[, -1, drop = FALSE]
  cpm <- sweep(collapse_ensembl_counts_to_symbol_sum(as.matrix(cnt)), 2, colSums(cnt), "/") * 1e6
  tmp_rda <- tempfile(fileext = ".rda")
  writeBin(readBin(gzfile(fp_dds, "rb"), "raw", n = 1e9), tmp_rda)
  e242912 <- new.env()
  load(tmp_rda, envir = e242912)
  bilat_dds <- e242912$bilat_dds
  cd <- as.data.frame(SummarizedExperiment::colData(bilat_dds))
  aligned <- align_pd_to_expr(cpm, cd)
  cpm <- aligned$expr
  cd <- aligned$pd
  x <- score_bulk_matrix(cpm, ref_fshd)
  gse242912_df <- data.frame(
    sample_id = cd$sample_id,
    subject = cd$Subject,
    location = as.character(cd$location),
    b_flag = as.logical(cd$b_flag),
    score = x,
    stringsAsFactors = FALSE
  )
  gse242912_df$score_scaled <- as.numeric(scale(gse242912_df$score))

  comprehensive_df <- tryCatch(
    load_wellstone_comprehensive_df(out_dir),
    error = function(e) {
      message("Wellstone comprehensive_df unavailable: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(comprehensive_df)) {
    clin_meta <- comprehensive_df[, c(
      "sample_id", "STIR_RATING", "STIR_status", "FAT_FRACTION", "Cumulative Score"
    )]
    names(clin_meta) <- c(
      "sample_id", "STIR_RATING", "STIR_status", "fat_fraction", "path_score"
    )
    gse242912_df <- merge(
      gse242912_df[, setdiff(
        names(gse242912_df),
        c("STIR_RATING", "STIR_status", "fat_fraction", "path_score")
      ), drop = FALSE],
      clin_meta,
      by = "sample_id",
      all.x = TRUE
    )
  }

  g140 <- suppressMessages(GEOquery::getGEO("GSE140261", GSEMatrix = TRUE)[[1]])
  pd140 <- Biobase::pData(g140)
  pd140$title <- as.character(pd140$title)
  pd140$patient_id <- extract_char_field(pd140, "patient_id")
  clinical_v2 <- data.frame(
    subject = pd140$patient_id,
    dux4_rlogsum_v2 = suppressWarnings(as.numeric(extract_char_field(pd140, "dux4.rlogsum"))),
    fat_fraction_v2 = suppressWarnings(as.numeric(extract_char_field(pd140, "fat_fraction"))),
    stir_rating_v2 = suppressWarnings(as.numeric(extract_char_field(pd140, "stir_rating"))),
    pathology_v2 = extract_char_field(pd140, "pathology"),
    stringsAsFactors = FALSE
  )
  clinical_v2 <- clinical_v2[!duplicated(clinical_v2$subject), , drop = FALSE]

  lr_cols <- c("subject", "score", "score_scaled")
  if ("STIR_status" %in% names(gse242912_df)) {
    lr_cols <- c(
      lr_cols, "STIR_RATING", "STIR_status", "fat_fraction", "path_score"
    )
  }
  suffix_lr_cols <- function(df, sfx) {
    nm <- names(df)
    nm[nm != "subject"] <- paste0(nm[nm != "subject"], "_", sfx)
    nm <- sub("^STIR_RATING_", "stir_rating_", nm)
    nm <- sub("^STIR_status_", "stir_status_", nm)
    names(df) <- nm
    df
  }
  left <- suffix_lr_cols(gse242912_df[gse242912_df$location == "L", lr_cols, drop = FALSE], "L")
  right <- suffix_lr_cols(gse242912_df[gse242912_df$location == "R", lr_cols, drop = FALSE], "R")
  gse242912_lr_df <- merge(left, right, by = "subject")
  gse242912_lr_df$score_diff <- abs(gse242912_lr_df$score_L - gse242912_lr_df$score_R)
  gse242912_lr_df$score_scaled_diff <- abs(
    gse242912_lr_df$score_scaled_L - gse242912_lr_df$score_scaled_R
  )
  med_diff <- stats::median(gse242912_lr_df$score_diff, na.rm = TRUE)
  gse242912_lr_df$score_discordant <- gse242912_lr_df$score_diff > med_diff
  if (all(c("stir_status_L", "stir_status_R") %in% names(gse242912_lr_df))) {
    gse242912_lr_df$stir_discordant <- with(
      gse242912_lr_df,
      !is.na(stir_status_L) & !is.na(stir_status_R) & stir_status_L != stir_status_R
    )
    gse242912_lr_df$score_higher_on_stir_pos <- with(
      gse242912_lr_df,
      ifelse(
        stir_status_L == "STIR+" & stir_status_R == "STIR-",
        score_scaled_L > score_scaled_R,
        ifelse(
          stir_status_R == "STIR+" & stir_status_L == "STIR-",
          score_scaled_R > score_scaled_L,
          NA
        )
      )
    )
    gse242912_lr_df$signed_score_stir_pos_minus_neg <- with(
      gse242912_lr_df,
      ifelse(
        stir_status_L == "STIR+" & stir_status_R == "STIR-",
        score_scaled_L - score_scaled_R,
        ifelse(
          stir_status_R == "STIR+" & stir_status_L == "STIR-",
          score_scaled_R - score_scaled_L,
          NA_real_
        )
      )
    )
    gse242912_lr_df$score_stir_alignment <- with(
      gse242912_lr_df,
      ifelse(
        score_scaled_L > score_scaled_R &
          stir_status_L == "STIR+" & stir_status_R == "STIR-",
        "Higher score on STIR+ side",
        ifelse(
          score_scaled_R > score_scaled_L &
            stir_status_R == "STIR+" & stir_status_L == "STIR-",
          "Higher score on STIR+ side",
          ifelse(
            stir_status_L == stir_status_R,
            paste0("Both ", stir_status_L),
            "STIR discordant; higher score on STIR- side"
          )
        )
      )
    )
    if (all(c("fat_fraction_L", "fat_fraction_R") %in% names(gse242912_lr_df))) {
      gse242912_lr_df$signed_fat_high_minus_low <- with(
        gse242912_lr_df,
        ifelse(
          score_scaled_L > score_scaled_R,
          fat_fraction_L - fat_fraction_R,
          fat_fraction_R - fat_fraction_L
        )
      )
    }
    if (all(c("path_score_L", "path_score_R") %in% names(gse242912_lr_df))) {
      gse242912_lr_df$signed_path_high_minus_low <- with(
        gse242912_lr_df,
        ifelse(
          score_scaled_L > score_scaled_R,
          path_score_L - path_score_R,
          path_score_R - path_score_L
        )
      )
    }
  }
  gse242912_lr_df <- merge(gse242912_lr_df, clinical_v2, by = "subject", all.x = TRUE)
  gse242912_lr_df <- gse242912_lr_df[, !duplicated(names(gse242912_lr_df)), drop = FALSE]
}, error = function(e) message("GSE242912 failed: ", conditionMessage(e)))

save_checkpoint()
message("Checkpoint: GSE242912 bilateral saved.")

## --- GSE143452 snRNA (low-memory; no Seurat) ---------------------------------
message("GSE143452 snRNA ...")
snrna_df <- NULL
tryCatch({
  if (is.null(fshd_induced_genes)) {
    stop("FSHD-induced gene list unavailable (GSE143453 must load first)")
  }
  suppressMessages(GEOquery::getGEOSuppFiles("GSE143452", makeDirectory = TRUE, baseDir = val_base))
  fp <- pick_first_existing_file(
    file.path(val_base, "GSE143452"),
    c("Seurat\\.Normalize\\.Counts\\.csv\\.gz$")
  )
  snrna_df <- score_gse143452_snrna_lowmem(
    fp,
    ref_fshd,
    fshd_induced_genes = fshd_induced_genes,
    days = c("Day3", "Day5"),
    max_cells_per_group = 2000L
  )
}, error = function(e) message("GSE143452 failed: ", conditionMessage(e)))

## --- Summary table -----------------------------------------------------------
fmt_p <- function(w) {
  if (is.null(w)) return(NA_character_)
  format(w$p.value, digits = 3, scientific = TRUE)
}

n_by <- function(cohort, grp = NULL) {
  sub <- validation_plot_df[validation_plot_df$cohort == cohort, , drop = FALSE]
  if (is.null(grp)) return(nrow(sub))
  sum(sub$group == grp, na.rm = TRUE)
}

validation_summary_df <- data.frame(
  GEO = c(
    "GSE140261", "GSE122873", "GSE115650", "GSE56787", "GSE26852",
    "GSE174301", "GSE242912", "GSE115650+GSE140261", "GSE143452"
  ),
  Role = c(
    "Reference (bulk)", "Discovery (scRNA)", "Bulk validation", "Bulk validation",
    "Bulk validation", "Anatomical validation", "FSHD-only bilateral",
    "Longitudinal tracking", "snRNA validation"
  ),
  Tissue_model = c(
    "Muscle biopsy", "Cultured myocyte", "Muscle biopsy", "Muscle biopsy",
    "Muscle biopsy", "Primary myoblast (FSHD only, all days)",
    "TA muscle biopsy (bilateral)", "Muscle biopsy (paired visits)",
    "Myotube nuclei (snRNA, Day 3 + Day 5)"
  ),
  N = c(
    paste0(nrow(pheno_bulk), " samples"),
    "6898 cells",
    paste0(n_by("GSE115650"), " samples"),
    paste0(n_by("GSE56787"), " samples"),
    paste0(n_by("GSE26852"), " samples"),
    if (is.null(anatomical_df)) "0" else paste0(nrow(anatomical_df), " samples"),
    if (is.null(gse242912_df)) "0" else paste0(nrow(gse242912_df), " samples"),
    if (is.null(longitudinal_df)) "0" else paste0(nrow(longitudinal_df), " paired subjects"),
    if (is.null(snrna_df)) "0" else paste0(nrow(snrna_df), " nuclei")
  ),
  Contrast = c(
    "FSHD vs Control (training)", "Adverse vs favorable cell tails",
    "FSHD vs Control", "FSHD vs Control (biopsy only)", "FSHD vs normal muscle",
    "TA vs quadriceps (FSHD1/2, by day)", "L-R concordance; STIR+ vs STIR- (Wellstone MRI)",
    "Visit I vs Visit II (within subject)", "FSHD2 vs Control; FSHD-Hi vs Lo (DUX4+54 genes)"
  ),
  Normalization = c(
    "TPM", "log-normalized counts", "FPKM→TPM", "Counts→CPM", "Illumina→symbol, scaled",
    "RSEM TPM", "Counts→CPM", "FPKM→TPM / TPM", "Seurat normalized counts"
  ),
  Statistic = c(
    "Reference z-scores", "Wilcoxon (cell scores by disease)",
    fmt_p(validation_wilcox[["GSE115650"]]),
    fmt_p(validation_wilcox[["GSE56787"]]),
    fmt_p(validation_wilcox[["GSE26852"]]),
    if (is.null(anatomical_df)) NA_character_ else {
      f2d12 <- anatomical_df[
        anatomical_df$day_lab %in% c("Day 0", "Day 12") &
          anatomical_df$muscle %in% c("TA", "Quad"),
        ,
        drop = FALSE
      ]
      day_p <- c()
      for (dl in c("Day 0", "Day 12")) {
        sub <- f2d12[f2d12$day_lab == dl, , drop = FALSE]
        if (length(unique(sub$muscle)) >= 2L) {
          day_p <- c(
            day_p,
            tryCatch(
              stats::wilcox.test(
                sub$score_scaled[sub$muscle == "TA"],
                sub$score_scaled[sub$muscle == "Quad"],
                alternative = "greater"
              )$p.value,
              error = function(e) NA_real_
            )
          )
        }
      }
      if (length(day_p) > 0L) {
        paste0(
          "Day0/12 TA>Quad p=",
          paste(vapply(day_p, function(p) format(p, digits = 2), character(1)), collapse = ", ")
        )
      } else {
        NA_character_
      }
    },
    if (is.null(gse242912_df)) NA_character_ else if (!"STIR_status" %in% names(gse242912_df)) {
      lr <- gse242912_df
      left <- lr[lr$location == "L", c("subject", "score")]
      right <- lr[lr$location == "R", c("subject", "score")]
      names(left)[2] <- names(right)[2] <- "score"
      merged <- merge(left, right, by = "subject", suffixes = c("_L", "_R"))
      if (nrow(merged) >= 3L) {
        paste0("L-R r=", round(stats::cor(merged$score_L, merged$score_R), 2))
      } else {
        NA_character_
      }
    } else {
      stir_ok <- gse242912_df$STIR_status %in% c("STIR+", "STIR-")
      w <- tryCatch(
        stats::wilcox.test(
          gse242912_df$score_scaled[stir_ok & gse242912_df$STIR_status == "STIR+"],
          gse242912_df$score_scaled[stir_ok & gse242912_df$STIR_status == "STIR-"]
        ),
        error = function(e) NULL
      )
      stir_cor <- tryCatch(
        stats::cor(
          gse242912_df$score_scaled[stir_ok],
          gse242912_df$STIR_RATING[stir_ok],
          use = "complete.obs"
        ),
        error = function(e) NA_real_
      )
      paste0(
        "STIR+ vs STIR- p=", fmt_p(w),
        "; r(STIR rating)=", round(stir_cor, 2)
      )
    },
    if (is.null(longitudinal_df)) NA_character_ else {
      ok <- is.finite(longitudinal_df$score_delta) & is.finite(longitudinal_df$dux4_rlogsum_v2)
      if (sum(ok) >= 5L) {
        ct <- suppressWarnings(stats::cor.test(
          longitudinal_df$score_delta[ok],
          longitudinal_df$dux4_rlogsum_v2[ok],
          method = "spearman"
        ))
        paste0("paired n=", nrow(longitudinal_df), "; Δscore~DUX4 ρ=", round(ct$estimate, 2))
      } else {
        paste0("paired n=", nrow(longitudinal_df))
      }
    },
    if (is.null(snrna_df)) NA_character_ else {
      d5 <- snrna_df[snrna_df$day == "Day5", , drop = FALSE]
      hi_lo <- d5[d5$status == "FSHD2" & d5$fshd_state %in% c("FSHD-Hi", "FSHD-Lo"), , drop = FALSE]
      paste0(
        "D5 FSHD2 vs Ctrl p=",
        format(tryCatch(stats::wilcox.test(score ~ status, data = d5)$p.value, error = function(e) NA), digits = 2),
        "; Hi vs Lo n=", nrow(hi_lo)
      )
    }
  ),
  Notes = c(
    "Training reference for ref_fshd",
    "Primary PhenoMapR cell discovery (GSE122873)",
    "Independent biopsy cohort",
    "Myotube/myoblast samples excluded",
    "Other myopathies excluded",
    "Global scaling across all FSHD samples; Quad–TA pairs at Day 0 and Day 12",
    "STIR+/STIR- from Wellstone comprehensive_df; paired discordance n=5",
    "Visit II samples overlap reference training set",
    "Held-out snRNA; FSHD-Hi/Lo by DUX4 + 54 gene enrichment clustering"
  ),
  stringsAsFactors = FALSE
)

if (!is.null(validation_plot_df) && nrow(validation_plot_df) > 0) {
  bulk_val_cohorts <- c("GSE115650", "GSE56787", "GSE26852")
  validation_plot_df <- validation_plot_df[
    validation_plot_df$cohort %in% bulk_val_cohorts,
    ,
    drop = FALSE
  ]
  validation_plot_df$score_scaled <- ave(
    validation_plot_df$score,
    validation_plot_df$cohort,
    FUN = function(x) as.numeric(scale(x))
  )
}

save_checkpoint()
message("Wrote validation cache to ", out_dir)
print(validation_summary_df[, c("GEO", "Role", "N", "Statistic")])
