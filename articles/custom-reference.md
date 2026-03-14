# GSE205154: Custom survival-based reference with derive_reference_from_bulk

## Overview

[GSE205154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154)
provides bulk RNA-seq and clinical data for 289 primary and metastatic
PDAC samples. This vignette uses the **expression and info files in the
vignette directory** to:

1.  Load expression and survival phenotype.
2.  **Derive a custom gene z-score reference** from expression and
    overall survival with
    [`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md).
3.  Score samples with that reference.
4.  Keep **primary and metastatic analyses separate**: Kaplan–Meier by
    custom score for primary tumors only, then for metastatic samples
    only.

## 1. Load expression and phenotype from vignette directory

``` r
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
  info <- read.delim(info_path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(info)[colnames(info) == "Array"] <- "sample_id"
  info$tumor_type <- info$Tumor_Type
  info$survival_time <- as.numeric(info$OS_Time)
  info$survival_event <- as.integer(info$OS_Status)
  info <- info[!is.na(info$survival_time) & !is.na(info$survival_event), ]

  mat <- read.delim(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)
  mat$Name[mat$Name == "" | is.na(mat$Name)] <- NA
  mat <- mat[!is.na(mat$Name), ]
  mat <- mat[!duplicated(mat$Name), ]
  rownames(mat) <- mat$Name
  sample_cols <- grep("^GSM", colnames(mat), value = TRUE)
  bulk_mat <- as.matrix(mat[, sample_cols])
  mode(bulk_mat) <- "numeric"

  keep <- intersect(colnames(bulk_mat), info$sample_id)
  bulk_mat <- bulk_mat[, keep, drop = FALSE]
  pheno <- info[info$sample_id %in% keep, ]
  pheno <- pheno[match(keep, pheno$sample_id), ]
  rownames(pheno) <- pheno$sample_id
}
```

## 2. Derive custom reference from survival

Expression for
[`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
must be **samples × genes**. We pass `gene_axis = "cols"` so the
function does not transpose.

``` r
bulk_samples_rows <- t(bulk_mat)
rownames(bulk_samples_rows) <- as.character(rownames(bulk_samples_rows))
pheno_surv <- pheno[, c("sample_id", "survival_time", "survival_event")]
pheno_surv$sample_id <- as.character(pheno_surv$sample_id)
pheno_surv <- pheno_surv[pheno_surv$sample_id %in% rownames(bulk_samples_rows), , drop = FALSE]
pheno_surv <- pheno_surv[match(rownames(bulk_samples_rows), pheno_surv$sample_id), , drop = FALSE]
pheno_surv <- pheno_surv[!is.na(pheno_surv$sample_id), , drop = FALSE]
bulk_samples_rows <- bulk_samples_rows[rownames(bulk_samples_rows) %in% pheno_surv$sample_id, , drop = FALSE]

ref_custom <- suppressWarnings(
  derive_reference_from_bulk(
    bulk_expression = bulk_samples_rows,
    phenotype = pheno_surv,
    sample_id_column = "sample_id",
    phenotype_type = "survival",
    survival_time = "survival_time",
    survival_event = "survival_event",
    gene_axis = "cols",
    verbose = TRUE
  )
)
```

    ## Using 289 samples common between expression and phenotype

    ## Cleaning gene symbols to approved HUGO IDs...

    ## Maps last updated on: Sat Nov 16 10:35:32 2024

    ## Collapsed to 38475 unique genes

    ## Expression looks like counts; applying log2(CPM+1)...

    ## Cox PH: 4546 genes had NA z-scores (convergence or low variation)

    ## Derived reference with 38475 genes

``` r
head(ref_custom)
```

    ##          survival_z
    ## A1BG      1.5793134
    ## A1BG-AS1 -0.3185781
    ## A1CF     -1.0188730
    ## A2M      -1.5102313
    ## A2M-AS1  -2.0662227
    ## A2ML1     3.2964255

## 3. Score samples with the custom reference

``` r
scores_custom <- PhenoMap(
  expression = bulk_mat,
  reference = ref_custom,
  z_score_cutoff = 2,
  verbose = TRUE
)
```

    ## Detected input type: matrix

    ## 5317 genes used for scoring against survival_z
    ## Calculating scores...
    ## Completed scoring for survival_z

``` r
col_custom <- grep("survival_z|weighted_sum", colnames(scores_custom), value = TRUE)[1]
dat <- pheno
dat$score_custom <- scores_custom[match(dat$sample_id, rownames(scores_custom)), col_custom]
```

## 4. Primary tumors only: Kaplan–Meier by custom score

``` r
suppressPackageStartupMessages(library(survival))
dat_primary <- dat[dat$tumor_type == "Primary", ]
dat_primary$custom_grp <- ifelse(
  dat_primary$score_custom >= median(dat_primary$score_custom, na.rm = TRUE),
  "High", "Low"
)
fit_primary <- survfit(Surv(survival_time, survival_event) ~ custom_grp, data = dat_primary)
plot(fit_primary, col = c("blue", "red"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Primary: Kaplan–Meier by custom survival-based score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("blue", "red"), lwd = 2, bty = "n")
```

![](custom-reference_files/figure-html/km-primary-1.png)

``` r
lr_primary <- survdiff(Surv(survival_time, survival_event) ~ custom_grp, data = dat_primary)
message("Primary — Log-rank p-value: ", round(1 - pchisq(lr_primary$chisq, 1), 4))
```

    ## Primary — Log-rank p-value: 0

## 5. Metastatic samples only: Kaplan–Meier by custom score

``` r
dat_met <- dat[dat$tumor_type == "Met", ]
dat_met$custom_grp <- ifelse(
  dat_met$score_custom >= median(dat_met$score_custom, na.rm = TRUE),
  "High", "Low"
)
fit_met <- survfit(Surv(survival_time, survival_event) ~ custom_grp, data = dat_met)
plot(fit_met, col = c("blue", "red"), lwd = 2, xlab = "Time", ylab = "Survival probability",
     main = "GSE205154 Metastatic: Kaplan–Meier by custom survival-based score")
legend("bottomleft", legend = c("Low score", "High score"), col = c("blue", "red"), lwd = 2, bty = "n")
```

![](custom-reference_files/figure-html/km-metastatic-1.png)

``` r
lr_met <- survdiff(Surv(survival_time, survival_event) ~ custom_grp, data = dat_met)
message("Metastatic — Log-rank p-value: ", round(1 - pchisq(lr_met$chisq, 1), 4))
```

    ## Metastatic — Log-rank p-value: 1e-04

## 6. References

- **GSE205154**: [Bulk RNA-Seq for 289 primary and metastatic PDAC
  samples](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205154)
  (GEO). Place `GSE205154.GPL20301.matrix.txt` and `GSE205154.info.txt`
  in the vignette directory.
- **PhenoMapR**:
  [`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
  for cohort-specific survival z-scores; see
  [`?derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md).

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
    ## [1] survival_3.8-6  PhenoMapR_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-4          jsonlite_2.0.0        dplyr_1.2.0          
    ##  [4] compiler_4.5.3        tidyselect_1.2.1      jquerylib_0.1.4      
    ##  [7] splines_4.5.3         systemfonts_1.3.2     textshaping_1.0.5    
    ## [10] yaml_2.3.12           fastmap_1.2.0         lattice_0.22-9       
    ## [13] R6_2.6.1              generics_0.1.4        knitr_1.51           
    ## [16] htmlwidgets_1.6.4     tibble_3.3.1          desc_1.4.3           
    ## [19] bslib_0.10.0          pillar_1.11.1         rlang_1.1.7          
    ## [22] cachem_1.1.0          splitstackshape_1.4.8 xfun_0.56            
    ## [25] fs_1.6.7              sass_0.4.10           otel_0.2.0           
    ## [28] cli_3.6.5             pkgdown_2.2.0         magrittr_2.0.4       
    ## [31] digest_0.6.39         grid_4.5.3            lifecycle_1.0.5      
    ## [34] vctrs_0.7.1           evaluate_1.0.5        glue_1.8.0           
    ## [37] data.table_1.18.2.1   HGNChelper_0.8.15     ragg_1.5.1           
    ## [40] rmarkdown_2.30        tools_4.5.3           pkgconfig_2.0.3      
    ## [43] htmltools_0.5.9
