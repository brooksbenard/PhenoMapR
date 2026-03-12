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
without matched outcomes. Here we show how to score an independent bulk
expression dataset (GSE253260) against the PRECOG primary pancreatic
signature.

``` r
if (!requireNamespace("googledrive", quietly = TRUE)) {
  stop("The 'googledrive' package is required for the GSE253260 example. Install with install.packages('googledrive').")
}

vignette_dir <- if (dir.exists("vignettes")) "vignettes" else if (dir.exists("Vignettes")) "Vignettes" else "."
if (!dir.exists(vignette_dir)) {
  dir.create(vignette_dir, recursive = TRUE, showWarnings = FALSE)
}

expr_path_253260 <- file.path(vignette_dir, "GSE253260_expression.rds")
info_path_253260 <- file.path(vignette_dir, "GSE253260_info.rds")

googledrive::drive_deauth() # Disable Google sign-in requirement

# Download expression matrix (GSE253260 bulk expression)
googledrive::drive_download(
  googledrive::as_id("1YuZQjGY6CTt-uicxRqYzp9t_tnIuQN4R"),
  path = expr_path_253260,
  overwrite = TRUE
)
```

    ## File downloaded:

    ## • GSE253260_expression.rds <id: 1YuZQjGY6CTt-uicxRqYzp9t_tnIuQN4R>

    ## Saved locally as:

    ## • ./GSE253260_expression.rds

``` r
# Download optional phenotype/info RDS (not required for scoring in this section)
googledrive::drive_download(
  googledrive::as_id("1Tpb8JC-2wO0Qppi5kM1vRuYTwBD1vY8f"),
  path = info_path_253260,
  overwrite = TRUE
)
```

    ## File downloaded:

    ## • GSE253260_info.rds <id: 1Tpb8JC-2wO0Qppi5kM1vRuYTwBD1vY8f>

    ## Saved locally as:

    ## • ./GSE253260_info.rds

``` r
bulk_253260 <- readRDS(expr_path_253260)
if (!is.matrix(bulk_253260) && !is.data.frame(bulk_253260)) {
  stop("GSE253260 object must be a matrix or data.frame with genes in rows")
}
if (is.data.frame(bulk_253260)) bulk_253260 <- as.matrix(bulk_253260)

message("GSE253260 expression: ", nrow(bulk_253260), " genes × ", ncol(bulk_253260), " samples")
```

    ## GSE253260 expression: 317 genes × 28858 samples

``` r
# Score samples with PRECOG primary pancreatic signature
scores_253260 <- PhenoMap(
  expression = bulk_253260,
  reference = "precog",
  cancer_type = "Pancreatic",
  verbose = TRUE
)
```

    ## Detected input type: matrix

    ## Warning in calculate_weighted_scores(expression_matrix = expr_info$matrix, : No
    ## common genes found between expression and reference data for Pancreatic

``` r
score_col_253260 <- grep("Pancreatic$", colnames(scores_253260), value = TRUE)[1]
if (!is.na(score_col_253260)) {
  print(plot_score_distribution(
    scores_253260,
    score_column = score_col_253260,
    main = "GSE253260: PRECOG Pancreatic score"
  ))
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] survminer_0.5.2 ggpubr_0.6.3    survival_3.8-6  ggplot2_4.0.2  
    ## [5] PhenoMapR_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10        generics_0.1.4     tidyr_1.3.2        rstatix_0.7.3     
    ##  [5] lattice_0.22-9     digest_0.6.39      magrittr_2.0.4     evaluate_1.0.5    
    ##  [9] grid_4.5.3         RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_2.0.0    
    ## [13] Matrix_1.7-4       backports_1.5.0    Formula_1.2-5      googledrive_2.1.2 
    ## [17] gridExtra_2.3      httr_1.4.8         purrr_1.2.1        scales_1.4.0      
    ## [21] textshaping_1.0.5  jquerylib_0.1.4    abind_1.4-8        cli_3.6.5         
    ## [25] rlang_1.1.7        splines_4.5.3      withr_3.0.2        cachem_1.1.0      
    ## [29] yaml_2.3.12        otel_0.2.0         tools_4.5.3        gargle_1.6.1      
    ## [33] ggsignif_0.6.4     dplyr_1.2.0        curl_7.0.0         broom_1.0.12      
    ## [37] vctrs_0.7.1        R6_2.6.1           lifecycle_1.0.5    fs_1.6.7          
    ## [41] car_3.1-5          htmlwidgets_1.6.4  ragg_1.5.1         pkgconfig_2.0.3   
    ## [45] desc_1.4.3         pkgdown_2.2.0      pillar_1.11.1      bslib_0.10.0      
    ## [49] gtable_0.3.6       glue_1.8.0         systemfonts_1.3.2  xfun_0.56         
    ## [53] tibble_3.3.1       tidyselect_1.2.1   knitr_1.51         farver_2.1.2      
    ## [57] htmltools_0.5.9    carData_3.0-6      rmarkdown_2.30     labeling_0.4.3    
    ## [61] compiler_4.5.3     S7_0.2.1
