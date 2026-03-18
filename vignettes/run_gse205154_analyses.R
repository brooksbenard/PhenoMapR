# Run GSE205154 vignette analyses and save plots (no pandoc required)
# Run from package root: Rscript Vignettes/run_gse205154_analyses.R

# Run from package root so Vignettes/ and data files are found
if (basename(getwd()) == "Vignettes") setwd("..")
if (file.exists("DESCRIPTION") && read.dcf("DESCRIPTION")[1, "Package"] == "PhenoMapR") {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
  devtools::load_all(".")
}
suppressPackageStartupMessages(library(PhenoMapR))
suppressPackageStartupMessages(library(survival))

vignette_dir <- "Vignettes"
info_path <- file.path(vignette_dir, "GSE205154.info.txt")
matrix_path <- file.path(vignette_dir, "GSE205154.GPL20301.matrix.txt")
if (!file.exists(info_path)) {
  vignette_dir <- "."
  info_path <- "GSE205154.info.txt"
  matrix_path <- "GSE205154.GPL20301.matrix.txt"
}

cat("\n========== Loading data ==========\n")
info <- read.delim(info_path, stringsAsFactors = FALSE, check.names = FALSE)
colnames(info)[colnames(info) == "Array"] <- "sample_id"
info$tumor_type <- info$Tumor_Type
info$survival_time <- as.numeric(info$OS_Time)
info$survival_event <- as.integer(info$OS_Status)
info <- info[!is.na(info$survival_time) & !is.na(info$survival_event), ]

mat <- read.delim(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)
# Use gene symbol (Name); keep first of duplicates and drop empty/NA
mat$Name[mat$Name == "" | is.na(mat$Name)] <- NA
mat <- mat[!is.na(mat$Name), ]
dup <- duplicated(mat$Name)
mat <- mat[!dup, ]
rownames(mat) <- mat$Name
sample_cols <- grep("^GSM", colnames(mat), value = TRUE)
bulk_mat <- as.matrix(mat[, sample_cols])
mode(bulk_mat) <- "numeric"

keep <- intersect(colnames(bulk_mat), info$sample_id)
bulk_mat <- bulk_mat[, keep, drop = FALSE]
pheno <- info[info$sample_id %in% keep, ]
pheno <- pheno[match(keep, pheno$sample_id), ]
rownames(pheno) <- pheno$sample_id

cat("Expression:", nrow(bulk_mat), "genes x", ncol(bulk_mat), "samples\n")
cat("Phenotype:", nrow(pheno), "samples (Primary:", sum(pheno$tumor_type == "Primary"), ", Met:", sum(pheno$tumor_type == "Met"), ")\n")

cat("\n========== PRECOG scoring ==========\n")
scores_primary <- score_expression(expression = bulk_mat, reference = "precog", cancer_type = "Pancreatic", verbose = TRUE)
scores_met <- score_expression(expression = bulk_mat, reference = "precog", cancer_type = "Pancreatic_Metastasis", verbose = TRUE)
col_pan <- grep("Pancreatic$", colnames(scores_primary), value = TRUE)[1]
col_met <- grep("Pancreatic_Metastasis", colnames(scores_met), value = TRUE)[1]
scores_df <- data.frame(
  sample_id = rownames(scores_primary),
  score_Pancreatic = if (!is.na(col_pan)) scores_primary[[col_pan]] else NA_real_,
  score_Pancreatic_Metastasis = if (!is.na(col_met)) scores_met[[col_met]] else NA_real_,
  stringsAsFactors = FALSE
)
dat <- merge(pheno, scores_df, by = "sample_id")

cat("\n========== Primary: KM by Pancreatic score ==========\n")
dat_primary <- dat[dat$tumor_type == "Primary", ]
dat_primary$score_grp <- ifelse(dat_primary$score_Pancreatic >= median(dat_primary$score_Pancreatic, na.rm = TRUE), "High", "Low")
fit_primary <- survfit(Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
plot(fit_primary, col = c("#2166AC", "#B2182B"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Primary: KM by PRECOG Pancreatic score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("#2166AC", "#B2182B"), lwd = 2, bty = "n")
lr_primary <- survdiff(Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
cat("Primary tumors — Log-rank p-value:", round(1 - pchisq(lr_primary$chisq, 1), 4), "\n")

cat("\n========== Metastatic: KM by Pancreatic_Metastasis score ==========\n")
dat_met <- dat[dat$tumor_type == "Met", ]
dat_met$score_grp <- ifelse(dat_met$score_Pancreatic_Metastasis >= median(dat_met$score_Pancreatic_Metastasis, na.rm = TRUE), "High", "Low")
fit_met <- survfit(Surv(survival_time, survival_event) ~ score_grp, data = dat_met)
plot(fit_met, col = c("#2166AC", "#B2182B"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Metastatic: KM by PRECOG Pancreatic_Metastasis score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("#2166AC", "#B2182B"), lwd = 2, bty = "n")
lr_met <- survdiff(Surv(survival_time, survival_event) ~ score_grp, data = dat_met)
cat("Metastatic — Log-rank p-value:", round(1 - pchisq(lr_met$chisq, 1), 4), "\n")

cat("\n========== Custom reference from survival ==========\n")
bulk_samples_rows <- t(bulk_mat)
rownames(bulk_samples_rows) <- as.character(rownames(bulk_samples_rows))
pheno_surv <- pheno[, c("sample_id", "survival_time", "survival_event")]
pheno_surv$sample_id <- as.character(pheno_surv$sample_id)
pheno_surv <- pheno_surv[pheno_surv$sample_id %in% rownames(bulk_samples_rows), , drop = FALSE]
pheno_surv <- pheno_surv[match(rownames(bulk_samples_rows), pheno_surv$sample_id), , drop = FALSE]
pheno_surv <- pheno_surv[!is.na(pheno_surv$sample_id), , drop = FALSE]
bulk_samples_rows <- bulk_samples_rows[rownames(bulk_samples_rows) %in% pheno_surv$sample_id, , drop = FALSE]
ref_custom <- derive_reference_from_bulk(
  bulk_expression = bulk_samples_rows,
  phenotype = pheno_surv,
  sample_id_column = "sample_id",
  phenotype_type = "survival",
  survival_time = "survival_time",
  survival_event = "survival_event",
  verbose = TRUE
)
cat("Custom reference: first 5 genes\n")
print(head(ref_custom, 5))

cat("\n========== Score with custom reference ==========\n")
scores_custom <- score_expression(expression = bulk_mat, reference = ref_custom, z_score_cutoff = 2, verbose = TRUE)
col_custom <- grep("survival_z|weighted_sum", colnames(scores_custom), value = TRUE)[1]
dat$score_custom <- scores_custom[match(dat$sample_id, rownames(scores_custom)), col_custom]

cat("\n========== Primary: KM by custom score ==========\n")
dat_primary <- dat[dat$tumor_type == "Primary", ]
dat_primary$custom_grp <- ifelse(dat_primary$score_custom >= median(dat_primary$score_custom, na.rm = TRUE), "High", "Low")
fit_primary_c <- survfit(Surv(survival_time, survival_event) ~ custom_grp, data = dat_primary)
pdf(file.path(fig_dir, "GSE205154_km_primary_custom.pdf"), width = 5, height = 4)
plot(fit_primary_c, col = c("blue", "red"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Primary: KM by custom survival-based score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("blue", "red"), lwd = 2, bty = "n")
dev.off()
lr_primary_c <- survdiff(Surv(survival_time, survival_event) ~ custom_grp, data = dat_primary)
cat("Primary (custom) — Log-rank p-value:", round(1 - pchisq(lr_primary_c$chisq, 1), 4), "\n")

cat("\n========== Metastatic: KM by custom score ==========\n")
dat_met <- dat[dat$tumor_type == "Met", ]
dat_met$custom_grp <- ifelse(dat_met$score_custom >= median(dat_met$score_custom, na.rm = TRUE), "High", "Low")
fit_met_c <- survfit(Surv(survival_time, survival_event) ~ custom_grp, data = dat_met)
plot(fit_met_c, col = c("#2166AC", "#B2182B"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Metastatic: KM by custom survival-based score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("#2166AC", "#B2182B"), lwd = 2, bty = "n")
lr_met_c <- survdiff(Surv(survival_time, survival_event) ~ custom_grp, data = dat_met)
cat("Metastatic (custom) — Log-rank p-value:", round(1 - pchisq(lr_met_c$chisq, 1), 4), "\n")

cat("\nDone. KM plots are shown above (not saved). For HTML vignette with plots, render the Rmd files.\n")
