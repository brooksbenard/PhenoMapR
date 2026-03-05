# Download GSE205154 (bulk PDAC), clean, run PhenoMapR scoring + custom reference + KM
# Run from package root: Rscript data-raw/download_and_run_gse205154.R
# Requires: GEOquery, PhenoMapR, survival (optional: survminer for pretty KM)

# Load PhenoMapR from source if running from package root
if (file.exists("DESCRIPTION")) {
  d <- read.dcf("DESCRIPTION")
  if (d[1, "Package"] == "PhenoMapR") {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools", repos = "https://cloud.r-project.org")
    devtools::load_all(".")
  }
}
library(GEOquery)
library(PhenoMapR)

options(timeout = 300)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
dir.create("vignettes/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Download GEO series and TMM matrix ----
message("Downloading GSE205154 series and supplementary...")
getGEOSuppFiles("GSE205154", makeDirectory = TRUE, baseDir = "data-raw", filter_regex = "TMM")
gse <- getGEO("GSE205154", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
pheno <- Biobase::pData(eset)

# Parse tumor type (Primary vs Metastatic)
pheno$tumor_type <- ifelse(grepl("Primary", pheno$characteristics_ch1), "Primary", "Metastatic")
pheno$sample_id <- pheno$geo_accession

# Read TMM matrix (genes x samples). File has ID, ensembl_gene_id, hgnc_symbol, ..., then sample columns (ST-*)
tmm_path <- list.files("data-raw/GSE205154", pattern = "TMM", full.names = TRUE)
if (length(tmm_path) == 0) stop("TMM file not found. Check data-raw/GSE205154/")
message("Reading ", tmm_path[1])
tmm <- read.delim(
  if (grepl("\\.gz$", tmm_path[1])) gzfile(tmm_path[1]) else tmm_path[1],
  check.names = FALSE, stringsAsFactors = FALSE
)
# Use HGNC symbols as rownames so they match PRECOG reference
if (!"hgnc_symbol" %in% colnames(tmm)) stop("TMM file must have hgnc_symbol column")
sample_cols <- grep("^ST-", colnames(tmm), value = TRUE)
if (length(sample_cols) == 0) sample_cols <- setdiff(colnames(tmm), c("ID", "ensembl_gene_id", "hgnc_symbol", "gene_biotype", "gene_name", "TMM_pass_filter"))
bulk_mat <- as.matrix(tmm[, sample_cols])
rownames(bulk_mat) <- tmm$hgnc_symbol
# Remove duplicate symbols (keep first) and empty/NA symbols
dup <- duplicated(rownames(bulk_mat)) | rownames(bulk_mat) == "" | is.na(rownames(bulk_mat))
bulk_mat <- bulk_mat[!dup, ]
mode(bulk_mat) <- "numeric"
# Map TMM sample IDs (ST-*) to GEO accessions if present in pheno
if (!any(pheno$sample_id %in% colnames(bulk_mat))) {
  if ("title" %in% colnames(pheno)) {
    match_title <- match(colnames(bulk_mat), pheno$title)
    if (any(!is.na(match_title))) {
      new_names <- colnames(bulk_mat)
      new_names[!is.na(match_title)] <- pheno$sample_id[match_title[!is.na(match_title)]]
      colnames(bulk_mat) <- new_names
    }
  }
}
keep <- intersect(colnames(bulk_mat), pheno$sample_id)
if (length(keep) == 0) keep <- colnames(bulk_mat)
bulk_mat <- bulk_mat[, keep, drop = FALSE]
pheno <- pheno[pheno$sample_id %in% keep, ]
if (nrow(pheno) == 0) pheno <- data.frame(sample_id = keep, tumor_type = NA, stringsAsFactors = FALSE)
rownames(pheno) <- pheno$sample_id

# ---- 2. Synthetic survival for demonstration (GEO has no survival; real data in dbGaP phs003597) ----
set.seed(42)
# Simulate: higher PhenoMapR score -> worse survival (for demo only)
pheno$survival_time <- runif(nrow(pheno), 6, 60)
pheno$survival_event <- rbinom(nrow(pheno), 1, 0.6)
pheno_file <- "inst/extdata/GSE205154_phenotype.rds"
saveRDS(pheno[, c("sample_id", "tumor_type", "survival_time", "survival_event")], pheno_file)
message("Saved phenotype to ", pheno_file)

# ---- 3. Score with PRECOG Pancreatic and Pancreatic_Metastasis ----
message("Scoring with PRECOG Pancreatic...")
scores_pancreatic <- score_expression(
  expression = bulk_mat,
  reference = "precog",
  cancer_type = "Pancreatic",
  verbose = TRUE
)
message("Scoring with PRECOG Pancreatic_Metastasis...")
scores_met <- score_expression(
  expression = bulk_mat,
  reference = "precog",
  cancer_type = "Pancreatic_Metastasis",
  verbose = TRUE
)

# Get score columns by name (e.g. weighted_sum_score_precog_Pancreatic)
col_pan <- grep("Pancreatic$", colnames(scores_pancreatic), value = TRUE)[1]
col_met <- grep("Pancreatic_Metastasis", colnames(scores_met), value = TRUE)[1]
scores_df <- data.frame(
  sample_id = rownames(scores_pancreatic),
  score_Pancreatic = if (!is.na(col_pan)) scores_pancreatic[[col_pan]] else NA_real_,
  score_Pancreatic_Metastasis = if (!is.na(col_met)) scores_met[[col_met]] else NA_real_,
  stringsAsFactors = FALSE
)
dat <- merge(pheno, scores_df, by = "sample_id")
saveRDS(dat, "inst/extdata/GSE205154_scores.rds")
message("Saved scores to inst/extdata/GSE205154_scores.rds")

# ---- 4. Kaplan-Meier by score (primary tumors only) ----
dat_primary <- dat[dat$tumor_type == "Primary", ]
if (nrow(dat_primary) >= 10) {
  dat_primary$score_grp <- ifelse(
    dat_primary$score_Pancreatic >= median(dat_primary$score_Pancreatic),
    "High", "Low"
  )
  fit <- survival::survfit(
    survival::Surv(survival_time, survival_event) ~ score_grp,
    data = dat_primary
  )
  pdf("vignettes/figures/GSE205154_km_primary.pdf", width = 5, height = 4)
  plot(fit, col = c("blue", "red"), lwd = 2, xlab = "Time", ylab = "Survival probability",
       main = "GSE205154 Primary PDAC (score stratified)")
  legend("bottomleft", legend = c("Low score", "High score"), col = c("blue", "red"), lwd = 2, bty = "n")
  dev.off()
  message("Saved Kaplan-Meier to vignettes/figures/GSE205154_km_primary.pdf")
  lr <- survival::survdiff(survival::Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
  message("Log-rank p-value: ", round(1 - pchisq(lr$chisq, 1), 4))
}

# ---- 5. Custom reference from bulk + survival ----
# Expression: samples x genes for derive_reference_from_bulk
bulk_samples_rows <- t(bulk_mat)
rownames(bulk_samples_rows) <- as.character(rownames(bulk_samples_rows))
pheno_surv <- pheno[pheno$sample_id %in% rownames(bulk_samples_rows), c("sample_id", "survival_time", "survival_event")]
pheno_surv$sample_id <- as.character(pheno_surv$sample_id)
pheno_surv <- pheno_surv[match(rownames(bulk_samples_rows), pheno_surv$sample_id), , drop = FALSE]
pheno_surv <- pheno_surv[!is.na(pheno_surv$sample_id), , drop = FALSE]
bulk_samples_rows <- bulk_samples_rows[rownames(bulk_samples_rows) %in% pheno_surv$sample_id, , drop = FALSE]
if (nrow(pheno_surv) < 10) {
  message("Skipping custom reference (need >= 10 samples with phenotype); have ", nrow(pheno_surv))
} else {
  message("Deriving custom reference from expression + survival (n=", nrow(pheno_surv), ")...")
  ref_custom <- derive_reference_from_bulk(
    bulk_expression = bulk_samples_rows,
    phenotype = pheno_surv,
    sample_id_column = "sample_id",
    phenotype_type = "survival",
    survival_time = "survival_time",
    survival_event = "survival_event",
    gene_axis = "cols",
    verbose = TRUE
  )
  saveRDS(ref_custom, "inst/extdata/GSE205154_custom_reference.rds")
  message("Saved custom reference to inst/extdata/GSE205154_custom_reference.rds")
  scores_custom <- score_expression(
    expression = bulk_mat,
    reference = ref_custom,
    z_score_cutoff = 2,
    verbose = TRUE
  )
  saveRDS(scores_custom, "inst/extdata/GSE205154_scores_custom.rds")
  message("Saved custom scores to inst/extdata/GSE205154_scores_custom.rds")
}
message("Done. Outputs in inst/extdata/ and vignettes/figures/")
