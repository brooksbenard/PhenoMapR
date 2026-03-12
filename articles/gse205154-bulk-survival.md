# Scoring bulk gene expression samples with PhenoMapR

## Overview

This vignette is an **example of using PhenoMapR to score bulk RNA-seq
samples** and to **show that the resulting score stratifies outcome**.
PhenoMapR has built-in pan-cancer outcomes meta-z score signatures from
[PRECOG](https://precog.stanford.edu/). To demonstrate that PhenoMapR
can assign prognostic risk to patient samples profiled by bulk RNAseq,
we use the
[GSE205154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154)
dataset: **289 primary and metastatic pancreatic ductal adenocarcinoma
(PDAC)** samples with gene expression and overall survival (Sears et
al.). Pre-processed gene expression and cleaned clinical annotations
were obtained from [PRECOG](https://precog.stanford.edu/). After scoring
samples with PhenoMapR’s built-in PRECOG references, we stratify by
score (high vs low) and demonstrate that the PhenoMapR score separates
survival in both primary and metastatic groups.

**Data files** (downloadable from
[PRECOG](https://precog.stanford.edu/)):

- `GSE205154.GPL20301.matrix.txt` — normalized expression (genes ×
  samples)
- `GSE205154.info.txt` — sample annotations (Tumor_Type: Primary / Met,
  OS_Time, OS_Status, Specimen_Site when available)

## 1. Load bulk expression and clinical annotations for GSE205154

Load the expression matrix and clinical info so we can score samples
with PhenoMapR and link scores to outcome.

``` r
knitr::opts_chunk$set(fig.width = 6, fig.height = 6, fig.align = "center")
suppressPackageStartupMessages(library(PhenoMapR))

# Vignette data: local paths, then Google Drive (https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n)
vignette_dir <- if (dir.exists("vignettes")) "vignettes" else if (dir.exists("Vignettes")) "Vignettes" else "."
info_path <- file.path(vignette_dir, "GSE205154.info.txt")
matrix_path <- file.path(vignette_dir, "GSE205154.GPL20301.matrix.txt")
if (!file.exists(info_path)) {
  info_path <- "GSE205154.info.txt"
  matrix_path <- "GSE205154.GPL20301.matrix.txt"
}
if (!file.exists(matrix_path) && !nzchar(Sys.getenv("CI", "")) && requireNamespace("googledrive", quietly = TRUE)) {
  googledrive::drive_deauth()
  googledrive::drive_download(googledrive::as_id("1Vk4KCQWF9ikpAuMsjFzVDCoy1TzDl2rN"), matrix_path, overwrite = TRUE)
}
if (!file.exists(info_path) && !nzchar(Sys.getenv("CI", "")) && requireNamespace("googledrive", quietly = TRUE)) {
  googledrive::drive_deauth()
  googledrive::drive_download(googledrive::as_id("1omAA2kfVn-nyyZfcc4vBhRFogC6cuoNQ"), info_path, overwrite = TRUE)
}
if (!file.exists(matrix_path)) {
  u <- Sys.getenv("PHENOMAPR_GSE205154_MATRIX_URL", "")
  if (nzchar(u)) tryCatch({ download.file(u, matrix_path, mode = "wb", quiet = TRUE) }, error = function(e) NULL)
}
if (!file.exists(info_path)) {
  u <- Sys.getenv("PHENOMAPR_GSE205154_INFO_URL", "")
  if (nzchar(u)) tryCatch({ download.file(u, info_path, mode = "wb", quiet = TRUE) }, error = function(e) NULL)
}
has_data <- file.exists(info_path) && file.exists(matrix_path)
knitr::opts_chunk$set(eval = has_data)
if (!has_data) {
  message("GSE205154 data files not found. See Vignettes/README.md for download instructions.")
} else {
  # Clinical annotations: Array = sample ID, Tumor_Type = Primary/Met, OS_Time, OS_Status, Specimen_Site (when available)
  info <- read.delim(info_path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info)[colnames(info) == "Array"] <- "sample_id"
  info$tumor_type <- info$Tumor_Type
  site_col <- intersect(c("Specimen_Site", "Specimen.Site"), colnames(info))[1]
  if (!is.na(site_col)) info$specimen_site <- info[[site_col]]
  info$survival_time <- as.numeric(info$OS_Time) / 365.25  # convert to years
  info$survival_event <- as.integer(info$OS_Status)
  info <- info[!is.na(info$survival_time) & !is.na(info$survival_event), ]

  # Expression matrix: Gene (ENSG), Name (symbol), then GSM columns
  mat <- read.delim(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)
  # Use gene symbol (Name); keep first of duplicates, drop empty/NA
  mat$Name[mat$Name == "" | is.na(mat$Name)] <- NA
  mat <- mat[!is.na(mat$Name), ]
  mat <- mat[!duplicated(mat$Name), ]
  rownames(mat) <- mat$Name
  sample_cols <- grep("^GSM", colnames(mat), value = TRUE)
  bulk_mat <- as.matrix(mat[, sample_cols])
  mode(bulk_mat) <- "numeric"

  # Align to samples in clinical info
  keep <- intersect(colnames(bulk_mat), info$sample_id)
  bulk_mat <- bulk_mat[, keep, drop = FALSE]
  pheno <- info[info$sample_id %in% keep, ]
  pheno <- pheno[match(keep, pheno$sample_id), ]
  rownames(pheno) <- pheno$sample_id

  message("Expression: ", nrow(bulk_mat), " genes × ", ncol(bulk_mat), " samples")
  message("Phenotype: ", nrow(pheno), " samples (Primary: ", sum(pheno$tumor_type == "Primary"), ", Met: ", sum(pheno$tumor_type == "Met"), ")")
}
```

    ## Expression: 38562 genes × 289 samples

    ## Phenotype: 289 samples (Primary: 218, Met: 71)

## 2. Score bulk samples with PhenoMapR

Score all samples with PhenoMapR using the built-in meta-z score
signatures for Adult PRECOG **Pancreatic** (primary) and
**Pancreatic_Metastasis**. By assigning a prognostic risk score to the
samples, we should be able to identify patients expressing more
adversely prognostic signatures as well as more favorably prognostic
signatures and this should stratify outcomes.

``` r
# To assign a PhenoMap score to each sample, pass the bulk gene expression matrix
# to the expression argument, select the reference database to score against (can
# be one of: precog, tcga, pediatric_precog, or ici_precog, and specify the cancer
# type to score against (use list_cancer_types(reference = "precog") to see available
# cancer types)

# Score all samples using the Primary PAAD signature
scores_primary <- PhenoMap(
  expression = bulk_mat,
  reference = "precog",
  cancer_type = "Pancreatic",
  verbose = TRUE
)
```

    ## Detected input type: matrix

    ## 7376 genes used for scoring against PancreaticCalculating scores...
    ## Completed scoring for Pancreatic

``` r
# Score all samples using the Metastatic PAAD signature
scores_met <- PhenoMap(
  expression = bulk_mat,
  reference = "precog",
  cancer_type = "Pancreatic_Metastasis",
  verbose = TRUE
)
```

    ## Detected input type: matrix

    ## 2997 genes used for scoring against Pancreatic_MetastasisCalculating scores...
    ## Completed scoring for Pancreatic_Metastasis

``` r
# make a combined dataframe of sample clinical annotations and PhenoMap scores 
col_pan <- grep("Pancreatic$", colnames(scores_primary), value = TRUE)[1]
col_met <- grep("Pancreatic_Metastasis", colnames(scores_met), value = TRUE)[1]
scores_df <- data.frame(
  sample_id = rownames(scores_primary),
  score_Pancreatic = if (!is.na(col_pan)) scores_primary[[col_pan]] else NA_real_,
  score_Pancreatic_Metastasis = if (!is.na(col_met)) scores_met[[col_met]] else NA_real_,
  stringsAsFactors = FALSE
)
dat <- merge(pheno, scores_df, by = "sample_id")
```

``` r
# Summary histogram: Primary and Metastatic PhenoMapR scores (z-scaled)
suppressPackageStartupMessages(library(ggplot2))
score_long <- data.frame(
  score = c(dat$score_Pancreatic, dat$score_Pancreatic_Metastasis),
  type = rep(c("Pancreatic (primary)", "Pancreatic_Metastasis"), each = nrow(dat))
)
score_long$score_z <- as.numeric(scale(score_long$score))
ggplot(score_long, aes(x = score_z, fill = type)) +
  geom_histogram(alpha = 0.8, position = "identity", bins = 30) +
  scale_fill_manual(values = c("Pancreatic (primary)" = "tan", "Pancreatic_Metastasis" = "darkred")) +
  labs(x = "PhenoMapR score (z-scaled)", y = "Count", fill = "Signature", title = "PhenoMapR score distribution") +
  theme_minimal()
```

![](gse205154-bulk-survival_files/figure-html/score-histogram-1.png)

## 3. Primary vs Metastatic outcomes

Before seeing if PhenoMapR scores stratify outcomes, let’s see if the
primary and metastatic patients in the study show differences in
outcomes.

``` r
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
# Unnamed vector: order matches strata (Met first, then Primary). Named palettes trigger scale warnings.
pal_tumor <- c("darkred", "tan")  # Metastatic, Primary

fit_primary <- survfit(Surv(survival_time, survival_event) ~ tumor_type, data = dat)
lr_primary <- survdiff(Surv(survival_time, survival_event) ~ tumor_type, data = dat)
pval_primary <- 1 - pchisq(lr_primary$chisq, 1)
cox_primary <- coxph(Surv(survival_time, survival_event) ~ tumor_type, data = dat)
hr_primary <- exp(coef(cox_primary))[1]
ci_primary <- as.vector(exp(confint(cox_primary)))
label_primary <- sprintf("p = %s\nHR = %.2f (95%% CI: %.2f-%.2f)",
  format.pval(pval_primary, digits = 2, eps = 0.001), hr_primary, ci_primary[1], ci_primary[2])
max_time_primary <- max(dat$survival_time, na.rm = TRUE)
ggsurvplot(fit_primary, data = dat, palette = pal_tumor, risk.table = FALSE,
           title = "Survival by Primary vs Metastatic status",
           xlab = "Time (years)", ylab = "Survival probability", legend.title = "Status",
           legend.labs = c("Metastatic", "Primary"), legend = "right",
           pval = label_primary, pval.coord = c(max_time_primary * 0.5, 0.95), pval.size = 3.5)
```

![](gse205154-bulk-survival_files/figure-html/primary-vs-met-1.png)

It seems like patients with metastatic disease trend towards worse
outcomes but it’s not statistically significant.

## 4. Primary tumors: PhenoMapR stratifies outcome

Stratify primary samples by PhenoMapR score (high vs low median split)
and show that the **Pancreatic** score stratifies survival.

``` r
# Unnamed vector: order matches strata (High first, then Low). Named palettes can trigger scale warnings.
pal_km <- c("#B2182B", "#2166AC")  # High (adverse), Low (favorable)
dat_primary <- dat[dat$tumor_type == "Primary", ]
dat_primary$score_grp <- ifelse(
  dat_primary$score_Pancreatic >= median(dat_primary$score_Pancreatic, na.rm = TRUE),
  "High", "Low"
)
fit_primary <- survfit(Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
lr_primary <- survdiff(Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
pval_primary <- 1 - pchisq(lr_primary$chisq, 1)
cox_primary <- coxph(Surv(survival_time, survival_event) ~ score_grp, data = dat_primary)
hr_primary <- exp(coef(cox_primary))[1]
ci_primary <- as.vector(exp(confint(cox_primary)))
label_primary <- sprintf("p = %s\nHR = %.2f (95%% CI: %.2f-%.2f)",
  format.pval(pval_primary, digits = 2, eps = 0.001), hr_primary, ci_primary[1], ci_primary[2])
max_time_primary <- max(dat_primary$survival_time, na.rm = TRUE)
ggsurvplot(fit_primary, data = dat_primary, palette = pal_km, risk.table = FALSE,
           title = "Survival by PhenoMapR score in Primary tumors",
           xlab = "Time (years)", ylab = "Survival probability", legend.title = "PhenoMapR Score",
           legend.labs = c("High", "Low"), legend = "right",
           pval = label_primary, pval.coord = c(max_time_primary * 0.5, 0.95), pval.size = 3.5)
```

![](gse205154-bulk-survival_files/figure-html/km-primary-1.png)

## 5. Metastatic samples: PhenoMapR stratifies outcome

Stratify metastatic samples by PhenoMapR score (high vs low median
split) and show that the **Pancreatic_Metastasis** score stratifies
survival in this subset.

``` r
dat_met <- dat[dat$tumor_type == "Met", ]
dat_met$score_grp <- ifelse(
  dat_met$score_Pancreatic_Metastasis >= median(dat_met$score_Pancreatic_Metastasis, na.rm = TRUE),
  "High", "Low"
)
fit_met <- survfit(Surv(survival_time, survival_event) ~ score_grp, data = dat_met)
lr_met <- survdiff(Surv(survival_time, survival_event) ~ score_grp, data = dat_met)
pval_met <- 1 - pchisq(lr_met$chisq, 1)
cox_met <- coxph(Surv(survival_time, survival_event) ~ score_grp, data = dat_met)
hr_met <- exp(coef(cox_met))[1]
ci_met <- as.vector(exp(confint(cox_met)))
label_met <- sprintf("p = %s\nHR = %.2f (95%% CI: %.2f-%.2f)",
  format.pval(pval_met, digits = 2, eps = 0.001), hr_met, ci_met[1], ci_met[2])
max_time_met <- max(dat_met$survival_time, na.rm = TRUE)
ggsurvplot(fit_met, data = dat_met, palette = pal_km, risk.table = FALSE,
           title = "Survival by PhenoMapR score in Metastatic Patients",
           xlab = "Time (years)", ylab = "Survival probability", legend.title = "PhenoMapR Score",
           legend.labs = c("High", "Low"), legend = "right",
           pval = label_met, pval.coord = c(max_time_met * 0.5, 0.95), pval.size = 3.5)
```

![](gse205154-bulk-survival_files/figure-html/km-metastatic-1.png)

## 6. Additional example: GSE253260 bulk expression

In many bulk expression studies, only the expression matrix is available
without matched outcomes. Here we show how to start from GEO, preprocess
expression and clinical data, and then score and stratify patients using
PhenoMapR.

``` r
study <- "GSE253260"
expr_tpm <- NULL
gset <- NULL

# ---- Primary: Google Drive RDS (ExpressionSet) ----
vignette_dir <- if (dir.exists("vignettes")) "vignettes" else if (dir.exists("Vignettes")) "Vignettes" else "."
if (!dir.exists(vignette_dir)) {
  dir.create(vignette_dir, recursive = TRUE, showWarnings = FALSE)
}
gse_rds_path <- file.path(vignette_dir, "GSE253260.rds")

primary_ok <- tryCatch({
  if (!requireNamespace("googledrive", quietly = TRUE)) stop("googledrive not available")
  googledrive::drive_deauth()
  googledrive::drive_download(
    googledrive::as_id("1ZBiFtZDTBKSpbLCgw6z_-qDF7WF2Sam5"),
    path = gse_rds_path,
    overwrite = TRUE
  )
  gset <- readRDS(gse_rds_path)
  if (is.list(gset) && length(gset) > 0 && !inherits(gset, "ExpressionSet")) {
    gset <- gset[[1]]
  }
  if (!inherits(gset, "ExpressionSet")) stop("RDS is not an ExpressionSet")
  expr <- Biobase::exprs(gset)
  if (!is.matrix(expr) || length(dim(expr)) != 2) stop("exprs(gset) not a matrix")
  if (nrow(expr) < ncol(expr) && nrow(expr) > 0) {
    expr <- t(expr)
    message("Transposed GSE253260 expression to genes × samples.")
  }
  fd <- Biobase::fData(gset)
  if (!is.null(fd) && nrow(fd) == nrow(expr)) {
    sym_col <- grep("symbol|SYMBOL|Gene.symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(sym_col)) {
      symbols <- as.character(fd[[sym_col]])
      if (requireNamespace("HGNChelper", quietly = TRUE)) {
        checked <- HGNChelper::checkGeneSymbols(symbols, unmapped.as.na = FALSE)
        if ("Suggested.Symbol" %in% names(checked)) symbols <- checked$Suggested.Symbol
      }
      symbols[symbols == "" | is.na(symbols)] <- NA
      if (sum(!is.na(symbols)) > 0) {
        keep <- !is.na(symbols)
        expr <- expr[keep, , drop = FALSE]
        symbols <- symbols[keep]
        rownames(expr) <- make.unique(symbols)
      }
    }
  }
  expr <- as.matrix(expr)
  expr[expr < 0] <- 0
  expr_tpm <<- as.matrix(apply(expr, 2, function(x) {
    denom <- sum(x, na.rm = TRUE)
    if (denom == 0) x else 1e6 * x / denom
  }))
  nrow(expr_tpm) > 0 && ncol(expr_tpm) > 0
}, error = function(e) {
  message("Primary (Drive RDS) failed: ", conditionMessage(e))
  FALSE
})
```

    ## File downloaded:

    ## • GSE253260.rds <id: 1ZBiFtZDTBKSpbLCgw6z_-qDF7WF2Sam5>

    ## Saved locally as:

    ## • ./GSE253260.rds

    ## Loading required namespace: Biobase

``` r
# ---- Fallback: GEO metadata + raw counts URL ----
if (!primary_ok && (is.null(expr_tpm) || nrow(expr_tpm) == 0 || ncol(expr_tpm) == 0)) {
  message("Downloading GEO metadata for: ", study)
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    geo_ok <- tryCatch({
      geo_data <- GEOquery::getGEO(GEO = study, GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
      if (length(geo_data) > 0) {
        gset <<- geo_data[[1]]
      }
      length(geo_data) > 0
    }, error = function(e) {
      message("getGEO failed: ", conditionMessage(e))
      FALSE
    })
  } else {
    geo_ok <- FALSE
  }

  if (geo_ok && requireNamespace("data.table", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    counts_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
      "&acc=", study,
      "&file=", study, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    )
    message("Attempting to download expression matrix from: ", counts_url)
    count_ok <- tryCatch({
      expression <- data.table::fread(counts_url, header = TRUE, colClasses = "integer")
      expression <- as.matrix(expression, rownames = 1)
      suppressPackageStartupMessages(library("org.Hs.eg.db", character.only = TRUE))
      gene_symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = rownames(expression),
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
      rownames(expression) <- gene_symbols
      expression <- expression[!is.na(rownames(expression)), , drop = FALSE]
      expression <- as.matrix(expression)
      expression[expression < 0] <- 0
      expr_tpm <<- as.matrix(apply(expression, 2, function(x) {
        denom <- sum(x, na.rm = TRUE)
        if (denom == 0) x else 1e6 * x / denom
      }))
      nrow(expr_tpm) > 0 && ncol(expr_tpm) > 0
    }, error = function(e) {
      message("Raw counts download/processing failed: ", conditionMessage(e))
      FALSE
    })
    if (count_ok) message("GSE253260 expression loaded from GEO raw counts (fallback).")
  }
}
```

    ## Downloading GEO metadata for: GSE253260

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

    ## Found 1 file(s)

    ## GSE253260_series_matrix.txt.gz

    ## Annotation GPL not available, so will use submitter GPL instead

    ## 

    ## Attempting to download expression matrix from: https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&acc=GSE253260&file=GSE253260_raw_counts_GRCh38.p13_NCBI.tsv.gz

    ## 'select()' returned 1:1 mapping between keys and columns

    ## GSE253260 expression loaded from GEO raw counts (fallback).

``` r
if (is.null(expr_tpm) || nrow(expr_tpm) == 0 || ncol(expr_tpm) == 0) {
  stop("Could not load GSE253260 expression (Drive RDS and GEO raw counts both failed or unavailable).")
}
if (is.null(gset)) {
  stop("GSE253260 clinical metadata (gset) not available; need Drive RDS or successful getGEO.")
}

message("GSE253260 TPM matrix: ", nrow(expr_tpm), " genes × ", ncol(expr_tpm), " samples")
```

    ## GSE253260 TPM matrix: 37674 genes × 317 samples

``` r
pheno <- Biobase::pData(gset)
pheno$geo_accession <- rownames(pheno)

pick_col <- function(df, patterns) {
  cols <- colnames(df)
  hit <- cols[apply(outer(patterns, cols, Vectorize(grepl, SIMPLIFY = FALSE), ignore.case = TRUE),
                    2, any)]
  if (length(hit)) hit[1] else NA_character_
}

id_col   <- pick_col(pheno, c("sample", "patient", "id", "accession"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
age_col  <- pick_col(pheno, c("age"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
sex_col  <- pick_col(pheno, c("sex", "gender"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
os_t_col <- pick_col(pheno, c("os_days", "os.time", "os", "overall.survival"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
os_e_col <- pick_col(pheno, c("os_censor", "os_event", "os.status", "os_event"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
pfs_t_col <- pick_col(pheno, c("pfs_days", "pfs.time", "pfs"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
pfs_e_col <- pick_col(pheno, c("pfs_censor", "pfs_event", "pfs.status"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
stage_col <- pick_col(pheno, c("stage", "disease", "resectable", "borderline", "locally", "meta"))
```

    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical
    ## Warning in FUN(newX[, i], ...): coercing argument of type 'list' to logical

``` r
sample_id <- if (!is.na(id_col)) as.character(pheno[[id_col]]) else pheno$geo_accession

clin <- data.frame(
  geo_accession = pheno$geo_accession,
  sample_id     = sample_id,
  age           = if (!is.na(age_col)) pheno[[age_col]] else NA,
  sex           = if (!is.na(sex_col)) pheno[[sex_col]] else NA,
  OS_Days       = if (!is.na(os_t_col)) as.numeric(pheno[[os_t_col]]) else NA_real_,
  OS_Censor     = if (!is.na(os_e_col)) as.integer(pheno[[os_e_col]]) else NA_integer_,
  PFS_Days      = if (!is.na(pfs_t_col)) as.numeric(pheno[[pfs_t_col]]) else NA_real_,
  PFS_Censor    = if (!is.na(pfs_e_col)) as.integer(pheno[[pfs_e_col]]) else NA_integer_,
  stage_raw     = if (!is.na(stage_col)) as.character(pheno[[stage_col]]) else NA_character_,
  stringsAsFactors = FALSE
)

stage_map <- function(x) {
  x_low <- tolower(x)
  dplyr::case_when(
    grepl("meta", x_low) ~ "Metastatic",
    grepl("loc", x_low)  ~ "Locally Advanced",
    grepl("border", x_low) ~ "Borderline",
    grepl("res", x_low)  ~ "Resectable",
    TRUE ~ NA_character_
  )
}
clin$stage_group <- factor(stage_map(clin$stage_raw),
                           levels = c("Resectable", "Borderline", "Locally Advanced", "Metastatic"))

common_samples <- intersect(colnames(expr_tpm), clin$sample_id)
clin <- clin[clin$sample_id %in% common_samples, , drop = FALSE]
clin <- clin[match(common_samples, clin$sample_id), , drop = FALSE]
rownames(clin) <- clin$sample_id
expr_tpm <- expr_tpm[, common_samples, drop = FALSE]
```

``` r
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

dat_os <- subset(clin, !is.na(OS_Days) & !is.na(OS_Censor) & !is.na(stage_group))
if (nrow(dat_os) >= 10) {
  fit_os <- survfit(Surv(OS_Days, OS_Censor) ~ stage_group, data = dat_os)
  ggsurvplot(
    fit_os,
    data = dat_os,
    title = "GSE253260: Overall survival by disease stage",
    xlab = "Time (days)",
    ylab = "Survival probability",
    legend.title = "Stage",
    risk.table = FALSE
  )
}

dat_pfs <- subset(clin, !is.na(PFS_Days) & !is.na(PFS_Censor) & !is.na(stage_group))
if (nrow(dat_pfs) >= 10) {
  fit_pfs <- survfit(Surv(PFS_Days, PFS_Censor) ~ stage_group, data = dat_pfs)
  ggsurvplot(
    fit_pfs,
    data = dat_pfs,
    title = "GSE253260: Progression-free survival by disease stage",
    xlab = "Time (days)",
    ylab = "PFS probability",
    legend.title = "Stage",
    risk.table = FALSE
  )
}
```

``` r
suppressPackageStartupMessages(library(PhenoMapR))

non_met_idx <- which(clin$stage_group != "Metastatic" & !is.na(clin$stage_group))
met_idx     <- which(clin$stage_group == "Metastatic")

scores_nonmet <- NULL
scores_met    <- NULL

if (length(non_met_idx) > 0) {
  expr_nonmet <- expr_tpm[, non_met_idx, drop = FALSE]
  scores_nonmet <- PhenoMap(
    expression  = expr_nonmet,
    reference   = "precog",
    cancer_type = "Pancreatic",
    verbose     = TRUE
  )
}
```

    ## Detected input type: matrix

    ## 7426 genes used for scoring against PancreaticCalculating scores...
    ## Completed scoring for Pancreatic

``` r
if (length(met_idx) > 0) {
  expr_met <- expr_tpm[, met_idx, drop = FALSE]
  scores_met <- PhenoMap(
    expression  = expr_met,
    reference   = "precog",
    cancer_type = "Pancreatic_Metastasis",
    verbose     = TRUE
  )
}
```

    ## Detected input type: matrix

    ## 3007 genes used for scoring against Pancreatic_MetastasisCalculating scores...
    ## Completed scoring for Pancreatic_Metastasis

``` r
clin$score_Pancreatic <- NA_real_
clin$score_Pancreatic_Met <- NA_real_

if (!is.null(scores_nonmet)) {
  col_pan <- grep("Pancreatic$", colnames(scores_nonmet), value = TRUE)[1]
  if (!is.na(col_pan)) {
    clin$score_Pancreatic[non_met_idx] <- scores_nonmet[, col_pan]
  }
}

if (!is.null(scores_met)) {
  col_met <- grep("Pancreatic_Metastasis", colnames(scores_met), value = TRUE)[1]
  if (!is.na(col_met)) {
    clin$score_Pancreatic_Met[met_idx] <- scores_met[, col_met]
  }
}
```

``` r
dat_nonmet <- subset(clin, stage_group != "Metastatic" & !is.na(stage_group) & !is.na(score_Pancreatic))

if (nrow(dat_nonmet) > 0) {
  ggplot(dat_nonmet, aes(x = stage_group, y = score_Pancreatic, fill = stage_group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.size = 0.5) +
    theme_minimal() +
    labs(
      title = "GSE253260: PRECOG Pancreatic score by disease stage (non-metastatic)",
      x = "Stage", y = "PhenoMapR score"
    ) +
    theme(legend.position = "none")
}
```

![](gse205154-bulk-survival_files/figure-html/gse253260-nonmet-plots-1.png)

``` r
dat_nonmet_os <- subset(dat_nonmet, !is.na(OS_Days) & !is.na(OS_Censor))
if (nrow(dat_nonmet_os) >= 10) {
  median_cut <- median(dat_nonmet_os$score_Pancreatic, na.rm = TRUE)
  dat_nonmet_os$score_grp <- ifelse(dat_nonmet_os$score_Pancreatic >= median_cut, "High", "Low")

  fit_nonmet <- survfit(Surv(OS_Days, OS_Censor) ~ score_grp, data = dat_nonmet_os)
  ggsurvplot(
    fit_nonmet,
    data = dat_nonmet_os,
    title = "GSE253260 (non-metastatic): OS by PRECOG Pancreatic score (median split)",
    xlab = "Time (days)", ylab = "Survival probability",
    legend.title = "Score group", legend.labs = c("High", "Low"),
    risk.table = FALSE
  )
}
```

``` r
dat_met <- subset(clin, stage_group == "Metastatic" & !is.na(score_Pancreatic_Met))

dat_met_os <- subset(dat_met, !is.na(OS_Days) & !is.na(OS_Censor))
if (nrow(dat_met_os) >= 10) {
  median_cut_met <- median(dat_met_os$score_Pancreatic_Met, na.rm = TRUE)
  dat_met_os$score_grp <- ifelse(dat_met_os$score_Pancreatic_Met >= median_cut_met, "High", "Low")

  fit_met <- survfit(Surv(OS_Days, OS_Censor) ~ score_grp, data = dat_met_os)
  ggsurvplot(
    fit_met,
    data = dat_met_os,
    title = "GSE253260 (metastatic): OS by PRECOG Pancreatic_Metastasis score (median split)",
    xlab = "Time (days)", ylab = "Survival probability",
    legend.title = "Score group", legend.labs = c("High", "Low"),
    risk.table = FALSE
  )
}
```

## 7. Summary

PhenoMapR can take a bulk expression dataset and assign prognostic risk
scores and significantly stratify outcomes in both primary and
metastatic disease in PAAD. In this case, we had access to outcomes
annotations for the cohort but it is reasonable to assume many datasets
profiled with bulk gene expression do not have outcomes data available
or reported. In these cases, PhenoMapR can nominate samples that are on
the extremes of the phenotype space (e.g. survival) and help define
sample groupings.

## 8. References

- **GSE205154**: Sears et al. Bulk RNA-Seq gene expression for 289
  primary and metastatic PDAC samples. *Nature Cancer* (2025). [GEO:
  GSE205154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154).
- **PRECOG**: Benard B et al. PRECOG update: an augmented resource of
  clinical outcome associations with gene expression for adult,
  pediatric, and immunotherapy cohorts. *Nucleic Acids Research* (2025).
  [precog.stanford.edu](https://precog.stanford.edu/); [PubMed:
  40909678](https://pubmed.ncbi.nlm.nih.gov/40909678/).

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] org.Hs.eg.db_3.22.0  AnnotationDbi_1.72.0 IRanges_2.44.0      
    ##  [4] S4Vectors_0.48.0     Biobase_2.70.0       BiocGenerics_0.56.0 
    ##  [7] generics_0.1.4       survminer_0.5.2      ggpubr_0.6.3        
    ## [10] survival_3.8-6       ggplot2_4.0.2        PhenoMapR_0.1.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.3.0                   gridExtra_2.3              
    ##  [3] httr2_1.2.2                 rlang_1.1.7                
    ##  [5] magrittr_2.0.4              otel_0.2.0                 
    ##  [7] matrixStats_1.5.0           compiler_4.5.3             
    ##  [9] RSQLite_2.4.6               png_0.1-8                  
    ## [11] systemfonts_1.3.2           vctrs_0.7.1                
    ## [13] crayon_1.5.3                pkgconfig_2.0.3            
    ## [15] fastmap_1.2.0               backports_1.5.0            
    ## [17] XVector_0.50.0              labeling_0.4.3             
    ## [19] rmarkdown_2.30              tzdb_0.5.0                 
    ## [21] ragg_1.5.1                  bit_4.6.0                  
    ## [23] purrr_1.2.1                 xfun_0.56                  
    ## [25] cachem_1.1.0                jsonlite_2.0.0             
    ## [27] blob_1.3.0                  DelayedArray_0.36.0        
    ## [29] broom_1.0.12                R6_2.6.1                   
    ## [31] bslib_0.10.0                RColorBrewer_1.1-3         
    ## [33] limma_3.66.0                car_3.1-5                  
    ## [35] GenomicRanges_1.62.1        jquerylib_0.1.4            
    ## [37] Seqinfo_1.0.0               SummarizedExperiment_1.40.0
    ## [39] knitr_1.51                  R.utils_2.13.0             
    ## [41] readr_2.2.0                 rentrez_1.2.4              
    ## [43] Matrix_1.7-4                splines_4.5.3              
    ## [45] tidyselect_1.2.1            abind_1.4-8                
    ## [47] yaml_2.3.12                 curl_7.0.0                 
    ## [49] lattice_0.22-9              tibble_3.3.1               
    ## [51] KEGGREST_1.50.0             withr_3.0.2                
    ## [53] S7_0.2.1                    evaluate_1.0.5             
    ## [55] desc_1.4.3                  xml2_1.5.2                 
    ## [57] Biostrings_2.78.0           pillar_1.11.1              
    ## [59] MatrixGenerics_1.22.0       carData_3.0-6              
    ## [61] hms_1.1.4                   scales_1.4.0               
    ## [63] glue_1.8.0                  tools_4.5.3                
    ## [65] data.table_1.18.2.1         ggsignif_0.6.4             
    ## [67] GEOquery_2.78.0             fs_1.6.7                   
    ## [69] XML_3.99-0.22               grid_4.5.3                 
    ## [71] tidyr_1.3.2                 googledrive_2.1.2          
    ## [73] Formula_1.2-5               cli_3.6.5                  
    ## [75] rappdirs_0.3.4              textshaping_1.0.5          
    ## [77] gargle_1.6.1                S4Arrays_1.10.1            
    ## [79] dplyr_1.2.0                 gtable_0.3.6               
    ## [81] R.methodsS3_1.8.2           rstatix_0.7.3              
    ## [83] sass_0.4.10                 digest_0.6.39              
    ## [85] SparseArray_1.10.9          htmlwidgets_1.6.4          
    ## [87] farver_2.1.2                memoise_2.0.1              
    ## [89] htmltools_0.5.9             pkgdown_2.2.0              
    ## [91] R.oo_1.27.1                 lifecycle_1.0.5            
    ## [93] httr_1.4.8                  statmod_1.5.1              
    ## [95] bit64_4.6.0-1
