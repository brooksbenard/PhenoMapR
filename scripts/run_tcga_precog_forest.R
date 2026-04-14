#!/usr/bin/env Rscript
#
# Bulk TCGA survival stratification using PhenoMapR PRECOG meta-z signatures.
# - For each TCGA cancer type, compute PhenoMapR weighted-sum score using the
#   concordant Adult PRECOG cancer type (from the package reference coverage map)
# - Stratify patients by the median score (high vs low)
# - Test survival with Cox PH + log-rank, then compile results into a forest plot
#
# This script expects TCGA inputs to already be present locally:
# - PANCAN clinical follow-up: data/tcga/clinical_PANCAN_patient_with_followup.tsv
#   OS: vital_status (Dead/Alive), days_to_death, days_to_last_followup.
#   Progression columns are empty for DLBC, TGCT, READ in this snapshot — OS only (see
#   results/tcga_outcomes_methods_note.txt).
# - Expression TPM files: data/tcga/TCGA_<TCGA_CODE>_tpm.fullIDs.remapped.tsv.gz
#
# It does NOT download from Google Drive (downloads can be environment-dependent).
# For TCGA-CDR outcomes (OS, DSS, DFI, PFI from TCGA-CDR-SupplementalTableS1.xlsx) use
# scripts/run_tcga_precog_forest_incremental.R instead.

suppressPackageStartupMessages({
  library(data.table)
  library(PhenoMapR)
  library(survival)
  library(ggplot2)
})

z_score_cutoff <- 2
reference <- "precog"

tpm_dir <- "data/tcga"
clinical_file <- file.path(tpm_dir, "clinical_PANCAN_patient_with_followup.tsv")
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clinical_file)) {
  stop("Missing PANCAN clinical file: ", clinical_file)
}

cat("Loading PANCAN clinical:", clinical_file, "\n")
clin <- fread(clinical_file, showProgress = FALSE)
req <- c("bcr_patient_barcode", "acronym", "vital_status", "days_to_death", "days_to_last_followup")
if (!all(req %in% names(clin))) {
  stop("PANCAN clinical TSV must include: ", paste(req, collapse = ", "))
}

na_pancan_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "[Not Available]", "[Not Applicable]", "NA", "Unknown")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

build_survival_pancan <- function(clin_sub) {
  if (nrow(clin_sub) == 0L) {
    return(data.frame(
      patient = character(),
      os_time = numeric(),
      os_event = integer(),
      stringsAsFactors = FALSE
    ))
  }
  vs <- tolower(trimws(as.character(clin_sub$vital_status)))
  dtd <- na_pancan_num(clin_sub$days_to_death)
  dtlfu <- na_pancan_num(clin_sub$days_to_last_followup)
  tt <- ifelse(vs == "dead", dtd, ifelse(vs == "alive", dtlfu, NA_real_))
  ev <- ifelse(vs == "dead", 1L, ifelse(vs == "alive", 0L, NA_integer_))
  out <- data.frame(
    patient = as.character(clin_sub$bcr_patient_barcode),
    os_time = tt,
    os_event = ev,
    stringsAsFactors = FALSE
  )
  ok <- is.finite(out$os_time) & !is.na(out$os_event)
  out <- out[ok, , drop = FALSE]
  out <- out[!duplicated(out$patient), , drop = FALSE]
  out
}

# Build concordance map (Adult PRECOG full name -> TCGA canonical code)
# Taken from inst/figures/make_reference_coverage.R in this repo.
precog_to_tcga <- c(
  Adrenocortical = "ACC", Appendix = "Appendix", Bladder = "BLCA", Brain_Astrocytoma = "LGG",
  Brain_Glioblastoma = "GBM", Brain_Glioma = "LGG", Brain_Medulloblastoma = "MB",
  Brain_Meningioma = "Meningioma", Brain_Neuroblastoma = "NB", Breast = "BRCA",
  Breast_Metastasis = "BRCA_Met", Cervical = "CESC", Colon = "COAD", Colon_Metastasis = "COAD_Met",
  Colon_Rectal = "READ", Desmoid = "Desmoid", Gastric = "STAD", Germ_cell_tumors = "TGCT",
  Head_and_neck = "HNSC", Head_and_neck_Larynx = "HNSC", Head_and_neck_Oral_SCC = "HNSC",
  Hematopoietic_AML = "LAML", Hematopoietic_B.ALL = "ALL",
  Hematopoietic_Burkitt_lymphoma = "Burkitt", Hematopoietic_CLL = "CLL",
  Hematopoietic_DLBCL = "DLBC", Hematopoietic_DLBCL_CHOP.treated = "DLBC",
  Hematopoietic_DLBCL_RCHOP.treated = "DLBC", Hematopoietic_FL = "FL",
  Hematopoietic_Mantle_cell_lymphoma = "MCL", Hematopoietic_Multiple_myeloma = "MM",
  Hematopoietic_Peripheral_T.cell_Lymphoma = "PTCL", Kidney = "KIRC", Liver = "LIHC",
  Liver_Primary = "LIHC", Lung_ADENO = "LUAD", Lung_LCC = "Lung_LCC", Lung_SCC = "LUSC",
  Lung_SCLC = "SCLC", Melanoma = "SKCM", Melanoma_Metastasis = "SKCM_Met", Mesothelioma = "MESO",
  Ovarian = "OV", Ovarian_Epithelial = "OV", Pancreatic = "PAAD",
  Pancreatic_Metastasis = "PAAD_Met", Prostate = "PRAD", Sarcoma_Ewing_sarcoma = "ESFT",
  Sarcoma_Liposarcoma = "SARC", Sarcoma_Osteosarcoma = "OSTEO", Sarcoma_Uterine = "UCS"
)

# For scoring we need a PRECOG label for each TCGA code present locally.
# If multiple PRECOG labels map to the same TCGA code, prefer primary/untreated labels.
concordant_precog_label <- function(tcga_code) {
  labels <- names(precog_to_tcga)[precog_to_tcga == tcga_code]
  if (length(labels) == 0) return(NA_character_)

  # Prefer non-metastatic, then non-primary-suffix labels.
  non_met <- labels[!grepl("Metastasis", labels, fixed = TRUE)]
  labels <- if (length(non_met) > 0) non_met else labels

  non_primary <- labels[!grepl("_Primary", labels, fixed = TRUE)]
  if (length(non_primary) > 0) labels <- non_primary

  labels[[1]]
}

sample_to_patient <- function(sample_ids) {
  parts <- strsplit(sample_ids, ".", fixed = TRUE)
  vapply(parts, function(x) paste(x[1], x[2], x[3], sep = "-"), character(1))
}

score_patient_level <- function(tpm_file, precog_label, z_score_cutoff) {
  # Correctness-first implementation:
  # load the full TPM matrix and score via PhenoMapR::PhenoMap().
  # This avoids subtle parsing/alignment issues when subsetting via awk.
  dt <- fread(tpm_file, showProgress = FALSE, fill = TRUE)
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]

  scores <- PhenoMap(
    expression = mat,
    reference = reference,
    cancer_type = precog_label,
    z_score_cutoff = z_score_cutoff,
    verbose = FALSE
  )

  score_col <- colnames(scores)[1]
  sample_ids <- rownames(scores)
  score_df <- data.frame(
    patient = sample_to_patient(sample_ids),
    score = as.numeric(scores[[score_col]]),
    stringsAsFactors = FALSE
  )
  aggregate(score ~ patient, data = score_df, FUN = mean)
}

cat("Scanning local TCGA TPM files in:", tpm_dir, "\n")
tpm_files <- list.files(
  tpm_dir,
  pattern = "^TCGA_[A-Z0-9]+_tpm\\.fullIDs\\.remapped\\.tsv\\.gz$",
  full.names = TRUE
)
if (length(tpm_files) == 0) {
  stop("No matching TPM files found in ", tpm_dir)
}

extract_tcga_code <- function(path) {
  sub("^TCGA_([^_]+)_tpm.*$", "\\1", basename(path))
}

tcga_codes <- unique(vapply(tpm_files, extract_tcga_code, character(1)))

cat("Found TCGA codes:", paste(tcga_codes, collapse = ", "), "\n")

results <- list()
for (tcga_code in sort(tcga_codes)) {
  precog_label <- concordant_precog_label(tcga_code)
  if (is.na(precog_label)) {
    cat("Skipping ", tcga_code, ": no concordant PRECOG label.\n", sep = "")
    next
  }

  tpm_file <- tpm_files[extract_tcga_code(tpm_files) == tcga_code][1]
  cat("\n---", tcga_code, "->", precog_label, "---\n")

  clin_sub <- clin[acronym == tcga_code]
  os_df <- build_survival_pancan(clin_sub)
  endpoint_lab <- "OS"
  if (nrow(os_df) < 20) {
    cat("Skipping ", tcga_code, ": too few ", endpoint_lab, " patients (n=", nrow(os_df), ").\n", sep = "")
    next
  }

  score_patient <- score_patient_level(tpm_file, precog_label, z_score_cutoff = z_score_cutoff)
  names(score_patient) <- c("patient", "score")

  dat <- merge(os_df, score_patient, by = "patient")
  if (nrow(dat) < 20) {
    cat("Skipping ", tcga_code, ": too few matched patients after merge (n=", nrow(dat), ").\n", sep = "")
    next
  }

  med <- median(dat$score, na.rm = TRUE)
  dat$score_grp <- factor(ifelse(dat$score >= med, "High", "Low"), levels = c("Low", "High"))

  # Cox + log-rank
  fit <- coxph(Surv(os_time, os_event) ~ score_grp, data = dat)
  coef_hi <- coef(fit)[["score_grpHigh"]]
  hr <- exp(coef_hi)
  ci <- exp(confint(fit)[1, ])

  lr <- survdiff(Surv(os_time, os_event) ~ score_grp, data = dat)
  p_lr <- 1 - pchisq(lr$chisq, 1)

  results[[tcga_code]] <- data.frame(
    tcga_code = tcga_code,
    precog_label = precog_label,
    survival_endpoint = endpoint_lab,
    n_patients = nrow(dat),
    median_score = med,
    hr = hr,
    ci_low = ci[1],
    ci_high = ci[2],
    p_logrank = p_lr,
    stringsAsFactors = FALSE
  )
}

if (length(results) == 0) {
  stop("No results to plot. Check file availability and concordance mapping.")
}

res_df <- do.call(rbind, results)
res_df <- res_df[order(res_df$hr, decreasing = TRUE), , drop = FALSE]
res_df$z_score_cutoff <- z_score_cutoff
res_df$reference <- reference

res_tsv <- file.path(out_dir, "tcga_precog_forest_results.tsv")
plot_data_tsv <- file.path(out_dir, "tcga_precog_forest_plot_data.tsv")
fwrite(res_df, res_tsv, sep = "\t")
fwrite(res_df, plot_data_tsv, sep = "\t")
cat("\nWrote results:", res_tsv, "\n")
cat("Wrote forest plot data (same rows as simple figure):", plot_data_tsv, "\n")
cat("For per-sample PANCAN + scores + full forest, run scripts/run_tcga_precog_forest_incremental.R\n")

# Forest plot (log scale on HR)
p <- ggplot(res_df, aes(y = reorder(tcga_code, hr), x = hr)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  labs(
    title = "TCGA survival by median PhenoMapR PRECOG score",
    x = "Hazard ratio (High vs Low; log scale)",
    y = "TCGA cancer type"
  ) +
  theme_minimal(base_size = 12)

forest_png <- file.path(out_dir, "tcga_precog_forest.png")
ggsave(forest_png, p, width = 7.5, height = 4.5, dpi = 200)
cat("Wrote forest plot:", forest_png, "\n")

cat("\nTop results (sorted by HR):\n")
print(res_df[, c("tcga_code", "precog_label", "survival_endpoint", "n_patients", "hr", "ci_low", "ci_high", "p_logrank")])

