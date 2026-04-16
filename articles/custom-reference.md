# Derive a custom reference signature with PhenoMapR

## Overview

Some rare cancer types are not represented in TCGA or PRECOG and are
therefore not available for the built-in reference scoring in
**`PhenoMapR`**. Additionally, you might have a proprietary dataset
(cancer or not) with bulk gene expression and sample annotations that
you would like to use as a phenotype signature reference. For these
cases, **`PhenoMapR`** provides functionality to take a gene expression
file and a sample annotation file and will generate gene-level z-scores
for your phenotype of interest. This is performed with the
**[`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)**
function.

Here, we use the
[GSE253260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253260)
dataset to demonstrate how **`PhenoMapR`** generates phenotype
signatures for use in downstream
**[`PhenoMap()`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md)**
scoring. GSE253260 contains bulk RNA-seq data for 317 PDAC samples with
survival and other clinical annotations. [\[1\]](#ref1)

This vignette uses the following pre-processed data:

- `GSE253260_expression.rds` (expression matrix, samples × genes)
- `GSE253260_info.rds` (clinical metadata)

1.  Load expression and sample annotations and visualize metadata.
2.  **Derive a custom gene z-score reference** from expression and
    overall survival with
    **[`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)**.
3.  Score samples with the resulting custom reference.
4.  Stratify GSE253260 overall survival by its **own** custom signature
    score (median split).
5.  Inspect cohort-wide score distribution by disease stage and overall
    survival by stage.
6.  **KM within each disease stage** using the cohort-wide score and a
    stage-specific median split.
7.  **Derive stage-specific** custom references and re-stratify survival
    within each stage by stage-specific median.

## 1. Load expression and phenotype data

Load the pre-processed bulk gene expression and sample annotations. This
dataset contains a range of patients across different disease stages,
and we will get back to this later in the vignette.

``` r
suppressPackageStartupMessages({
  library(PhenoMapR)
  library(ggplot2)
  library(dplyr)
})

# Vignette data: Google Drive
# requireNamespace("googledrive", quietly = TRUE)
googledrive::drive_deauth()
options(googledrive_quiet = TRUE)

# Load expression data
googledrive::drive_download(
  googledrive::as_id("1NXLzKx0O6q-9Bnl_qPwA0hhijIB9jk-6"),
  path = "GSE253260_expression.rds",
  overwrite = TRUE
)

bulk_mat <- readRDS("GSE253260_expression.rds")

# Load sample annotation data
googledrive::drive_download(
  googledrive::as_id("1j0syd8soOQNkudhGnc407TVCJuuvm4jk"),
  path = "GSE253260_info.rds",
  overwrite = TRUE
)

pheno <- readRDS("GSE253260_info.rds") %>%
  dplyr::filter(geo_accession %in% rownames(bulk_mat))

# Plot cohort demographics
dat_demo <- pheno %>%
  transmute(
    Age = suppressWarnings(as.numeric(Age)),
    Sex = as.character(Sex),
    ECOG = suppressWarnings(as.numeric(ECOG)),
    Stage = suppressWarnings(as.factor(Disease_Stage))
  )

if (nrow(dat_demo) > 0) {
 # Age histogram with mean line and text
mean_age <- mean(subset(dat_demo, is.finite(Age))$Age, na.rm = TRUE)

p_age <- ggplot(subset(dat_demo, is.finite(Age)), aes(x = Age)) +
  geom_histogram(bins = 30, fill = "grey75", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = mean_age, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = mean_age, y = Inf, 
           label = paste0("Mean = ", round(mean_age, 1)), 
           vjust = 1.5, hjust = 1, color = "black", size = 3.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "# of Patients", title = "Age")

# Sex bar plot with counts on top
p_sex <- ggplot(dat_demo, aes(x = Sex)) +
  geom_bar(fill = "grey55", color = "white", linewidth = 0.2) +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "# of Patients", title = "Sex") +
  ylim(0, max(table(dat_demo$Sex)) * 1.1)  # Add space for labels

# ECOG bar plot with counts on top
dat_ecog <- subset(dat_demo, is.finite(ECOG))
dat_ecog$ECOG <- factor(dat_ecog$ECOG)

p_ecog <- ggplot(dat_ecog, aes(x = ECOG)) +
  geom_bar(fill = "grey35", color = "white", linewidth = 0.2) +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3.5, color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "# of Patients", title = "ECOG") +
  ylim(0, max(table(dat_ecog$ECOG)) * 1.1)  # Add space for labels

# Stage bar plot with counts on top
dat_stage <- subset(dat_demo, is.finite(Stage))
dat_stage$Stage <- factor(dat_stage$Stage, levels = c("Resectable", "Borderline", "Locally Advanced", "Metastatic"))

p_stage <- ggplot(dat_stage, aes(x = Stage)) +
  geom_bar(fill = "grey15", color = "white", linewidth = 0.2) +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3.5, color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "# of Patients", title = "Disease Stage") +
  ylim(0, max(table(dat_stage$Stage)) * 1.1)  # Add space for labels
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    print(p_age + p_sex + p_ecog + p_stage + patchwork::plot_layout(widths = c(2, 1, 2, 2)))
  } else {
    print(p_age)
    print(p_sex)
    print(p_ecog)
    print(p_stage)
  }
}
```

![](custom-reference_files/figure-html/load-data%20and%20cohort-demographics-1.png)

## 2. Derive a cohort-wide custom reference signature

Before we perform any subset-specific analyses, we will first use all
patients in the dataset to derive a survival phenotype z-score signature
from the expression data.

Expression for
**[`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)**
must be **genes as rows and samples as columns**. Sample IDs in the
phenotype must match the **column names** of the expression matrix. For
the phenotype dataframe, specify the sample ID column, phenotype type
(here, overall survival), and the column names for time and event.

``` r
ref_custom <- suppressWarnings(derive_reference_from_bulk(
    # derive_reference_from_bulk expects genes × samples; downloaded matrix is samples × genes
    bulk_expression = t(bulk_mat),
    phenotype = pheno,
    sample_id_column = "geo_accession",
    phenotype_type = "survival",
    survival_time = "OS_Time",
    survival_event = "OS_Censor",
    verbose = TRUE
  )
)

# Vizualize the gene z-score signature distribution
PhenoMapR::plot_score_distribution(ref_custom)
```

![](custom-reference_files/figure-html/derive-ref-1.png)

Each dataset and phenotype will generate unique signatures, with
different numbers of genes significantly associated with the phenotype
of interest. For GSE253260, the resulting signature heatmap with top
genes annotated is as follows:

``` r
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Using the output of derive_reference_from_bulk() from above, plot the signature heatmap
plot_reference_signature(ref_custom, row_title = "Survival z-score")
```

![](custom-reference_files/figure-html/unnamed-chunk-1-1.png)

## 3. Score samples with the cohort-wide custom reference

Now that we have our custom phenotype signature for GSE253260, we use
this signature to score all samples in the dataset.

``` r
bulk_genes_rows <- t(bulk_mat)

scores_custom <- PhenoMap(
  expression = bulk_genes_rows,
  reference = ref_custom,
  z_score_cutoff = 2,
  verbose = TRUE
)
```

    ## 1730 genes used for scoring against survival_z
    ## Calculating scores...
    ## Completed scoring for survival_z

``` r
col_custom <- grep("survival_z|weighted_sum", colnames(scores_custom), value = TRUE)[1]

pheno$score_custom <- scores_custom[match(pheno$geo_accession, rownames(scores_custom)), col_custom]
```

Next, we plot the score distributions across the cohort, colored by
median split

``` r
suppressPackageStartupMessages(library(ggplot2))

dat_sc <- pheno[is.finite(pheno$score_custom) & !is.na(pheno$score_custom), , drop = FALSE]
if (nrow(dat_sc) > 0) {
  med_sc <- median(dat_sc$score_custom, na.rm = TRUE)
  dat_sc$score_grp <- factor(ifelse(dat_sc$score_custom >= med_sc, "High", "Low"), levels = c("Low", "High"))
  pal_km <- c("Low" = "#2166AC", "High" = "#B2182B")  # Low (favorable), High (adverse)

  p_hist <- ggplot(dat_sc, aes(x = score_custom)) +
    geom_histogram(aes(fill = score_grp), bins = 20, color = "white", linewidth = 0.2) +
    geom_density(aes(y = after_stat(count), color = score_grp), linewidth = 0.7) +
    scale_fill_manual(values = pal_km, drop = FALSE) +
    scale_color_manual(values = pal_km, drop = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    theme_minimal() +
    labs(
      x = "Cohort-wide custom signature score",
      y = "Number of samples",
      title = "Score distribution",
      fill = "Median split",
      color = "Median split"
    )

  dat_sc <- dat_sc[order(dat_sc$score_custom), , drop = FALSE]
  dat_sc$rank <- seq_len(nrow(dat_sc))
  p_rank <- ggplot(dat_sc, aes(x = rank, y = score_custom)) +
    geom_point(aes(color = score_grp), size = 1.2, alpha = 0.9) +
    scale_color_manual(values = pal_km, drop = FALSE) +
    theme_minimal() +
    labs(
      x = "Samples (ordered by score)",
      y = "Cohort-wide custom signature score",
      title = "Scores across samples",
      color = "Median split"
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    print(p_hist + p_rank + patchwork::plot_layout(widths = c(1.2, 1)))
  } else {
    print(p_hist)
    print(p_rank)
  }
}
```

![](custom-reference_files/figure-html/score-distribution-1.png)

## 4. Survival anlaysis by PhenoMapR score assignment

**`PhenoMapR`** scores samples based on the weighted-sum of a sample’s
gene expression levels and the weights from the phenotype reference
signature. Because the signature used to score samples in the dataset
was derived from the expression associations with survival across the
datasets, the **`PhenoMapR`** score should stratify outcomes. We perform
Kaplan–Meier analysis using a **`PhenoMapR`** median split.

``` r
suppressPackageStartupMessages(library(survival))
pheno$survival_time <- pheno$OS_Time
pheno$survival_event <- pheno$OS_Censor
dat_os <- subset(pheno, !is.na(survival_time) & !is.na(survival_event) & !is.na(score_custom))

  # Convert to months for plotting
  time_months <- dat_os$survival_time
  if (is.finite(max(time_months, na.rm = TRUE)) && max(time_months, na.rm = TRUE) > 400) {
    time_months <- time_months / 30.4375
  }
  dat_os$time_months <- time_months

  dat_os$custom_grp <- factor(
    ifelse(dat_os$score_custom >= median(dat_os$score_custom, na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High")
  )
  fit <- survfit(Surv(time_months, survival_event) ~ custom_grp, data = dat_os)
  lr <- survdiff(Surv(time_months, survival_event) ~ custom_grp, data = dat_os)
  p_lr <- 1 - pchisq(lr$chisq, df = 1)
  cox <- coxph(Surv(time_months, survival_event) ~ custom_grp, data = dat_os)
  hr <- exp(coef(cox))[1]
  ci <- as.vector(exp(confint(cox)))[1:2]
  label <- sprintf("p = %s\nHR = %.1f (%.1f-%.1f)",
                   format.pval(p_lr, digits = 2, eps = 0.001), hr, ci[1], ci[2])
  max_t <- max(dat_os$time_months, na.rm = TRUE)

  if (requireNamespace("survminer", quietly = TRUE)) {
    suppressPackageStartupMessages(library(survminer))
    pal_km <- c("Low" = "#2166AC", "High" = "#B2182B")  # Low (favorable), High (adverse)
    p_surv <- ggsurvplot(
      fit,
      data = dat_os,
      title = "Survival by PhenoMapR Score",
      xlab = "Time (months)",
      ylab = "Survival probability",
      legend.title = "Custom score",
      legend.labs = c("Low", "High"),
      legend = "right",
      palette = unname(pal_km[c("Low", "High")]),
      risk.table = FALSE,
      pval = label,
      pval.coord = c(max_t * 0.4, 0.95),
      pval.size = 3
    ) 
    
    p_surv$plot + theme(plot.title = element_text(hjust = 0.5))
  }
```

![](custom-reference_files/figure-html/km-custom-1.png) The sample
**`PhenoMapR`** scores successfully stratify outcomes across the cohort.

## 5. Overall survival by disease stage

As mentioned earlier (Section 1), this dataset contains patients with a
range of disease stages, reflecting significant differences in
prognostic risk that could influence the **`PhenoMapR`** score
assignment. Before we test whether the **`PhenoMapR`** score is
correlated with disease stage, we first confirm that the four different
disease stages have significantly different overall survival results.

``` r
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(dplyr))

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

pal_stage <- c(
  "Resectable" = "#74B9FF",
  "Borderline" = "#26A69A",
  "Locally Advanced" = "#F39C12",
  "Metastatic" = "#FF6B6B"
)

dat_stage <- data.frame(
  geo_accession = as.character(pheno$geo_accession),
  stage_group = factor(stage_map(pheno$Disease_Stage), levels = names(pal_stage)),
  OS_time_raw = suppressWarnings(as.numeric(pheno$OS_Time)),
  OS_event_raw = suppressWarnings(as.integer(pheno$OS_Censor)),
  stringsAsFactors = FALSE
)

# Heuristic: if time looks like days, convert to months for plotting (as in example figure)
time_months <- dat_stage$OS_time_raw
if (is.finite(max(time_months, na.rm = TRUE)) && max(time_months, na.rm = TRUE) > 400) {
  time_months <- time_months / 30.4375
}
dat_stage$time_months <- time_months

dat_stage <- subset(dat_stage, !is.na(stage_group) & !is.na(time_months) & !is.na(OS_event_raw))
dat_stage$stage_group <- droplevels(dat_stage$stage_group)

if (nrow(dat_stage) >= 10 && nlevels(dat_stage$stage_group) >= 2) {
  fit <- survfit(Surv(time_months, OS_event_raw) ~ stage_group, data = dat_stage)
  lr <- survdiff(Surv(time_months, OS_event_raw) ~ stage_group, data = dat_stage)
  p_lr <- 1 - pchisq(lr$chisq, df = (nlevels(dat_stage$stage_group) - 1))

  cox <- coxph(Surv(time_months, OS_event_raw) ~ stage_group, data = dat_stage)
  hr <- exp(coef(cox))
  ci <- exp(confint(cox))
  hr_lines <- paste0(
    gsub("^stage_group", "", names(hr)),
    ": HR ", sprintf("%.1f", hr),
    " (", sprintf("%.1f", ci[, 1]), "-", sprintf("%.1f", ci[, 2]), ")"
  )
  label <- paste0(
    "p ", format.pval(p_lr, digits = 2, eps = 0.001),
    "\n", paste(hr_lines, collapse = "\n")
  )
  max_t <- max(dat_stage$time_months, na.rm = TRUE)

  if (requireNamespace("survminer", quietly = TRUE)) {
    suppressPackageStartupMessages(library(survminer))
    pal_vec <- unname(pal_stage[levels(dat_stage$stage_group)])
    p_surv_4 <- ggsurvplot(
      fit,
      data = dat_stage,
      title = "Survival by Disease Stage",
      palette = pal_vec,
      risk.table = FALSE,
      legend.title = "Stage",
      legend = "right",
      legend.labs = levels(dat_stage$stage_group),
      xlab = "Time (months)",
      ylab = "Overall survival probability",
      pval = label,
      pval.coord = c(max_t * 0.25, 0.90),
      pval.size = 3
    )
    p_surv_4$plot + theme(plot.title = element_text(hjust = 0.5))
    # print(p_surv_4)
  } else {
    plot(fit, col = pal_stage[levels(dat_stage$stage_group)], lwd = 2,
         xlab = "Time (months)", ylab = "Overall survival probability")
    legend("topright", legend = levels(dat_stage$stage_group),
           col = pal_stage[levels(dat_stage$stage_group)], lwd = 2, bty = "n")
    mtext(paste0("p = ", format.pval(p_lr, digits = 2, eps = 0.001)), side = 4, adj = 1, line = -1, cex = 0.8)
  }
}
```

![](custom-reference_files/figure-html/km-by-stage-1.png) Indeed, the
disease stages provide robust survival stratification.

## 6. Cohort-wide custom score distribution by disease stage

To check if patients with a more adverse **`PhenoMapR`** score are
enriched in patients with more aggressive disease, we plot the
**`PhenoMapR`** score distribution across the four disease stage groups.

``` r
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggsignif))

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

pal_stage <- c(
  "Resectable" = "#74B9FF",
  "Borderline" = "#26A69A",
  "Locally Advanced" = "#F39C12",
  "Metastatic" = "#FF6B6B"
)

dat_hist <- dat_os %>%
  mutate(
    geo_accession = as.character(geo_accession),
    stage_group = factor(stage_map(pheno$Disease_Stage[match(geo_accession, pheno$geo_accession)]),
                         levels = names(pal_stage))
  ) %>%
  filter(!is.na(stage_group))

if (nrow(dat_hist) > 0) {
  library(ggsignif)

p_dist <- ggplot(dat_hist, aes(x = stage_group, y = scale(score_custom), fill = stage_group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = pal_stage, drop = FALSE) +
  geom_signif(
    comparisons = combn(levels(dat_hist$stage_group), 2, simplify = FALSE),
    test = "wilcox.test",
    map_signif_level = TRUE,  # TRUE converts p to *, **, ***
    textsize = 3,
    step_increase = 0.08,  # Vertical spacing between brackets
    tip_length = 0.01
  ) +
  theme_pubr(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.title      = element_text(hjust = 0.5)
  ) +
  labs(
    title = "PhenoMapR Score Across Stage",
    x     = NULL,
    y     = "Cohort-wide signature score",
    fill  = "Stage"
  )
p_dist
}
```

![](custom-reference_files/figure-html/custom-score-hist-1.png)
Interestingly, it does not appear that any of the disease risk groups
are enriched for favorable or adverse **`PhenoMapR`** scores. This
suggests that **`PhenoMapR`** identifies potentially stage agnostic
prognostic signals and might be able to refine the risk within each
stage.

## 7. Survival within each disease stage (stage-specific median split)

We next ask whether the cohort-wide custom score stratifies survival
**within** each disease stage. For each stage we use the **median of the
cohort-wide score within that stage**: samples above that stage’s median
are “High”, below are “Low”. KM is then run separately in each stage.

``` r
if (!requireNamespace("survminer", quietly = TRUE)) {
  message("Package 'survminer' not installed; skipping within-stage KM.")
} else if (!requireNamespace("patchwork", quietly = TRUE)) {
  message("Package 'patchwork' not installed; skipping within-stage KM.")
} else {
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(patchwork))
  suppressPackageStartupMessages(library(dplyr))

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
  pal_stage <- c(
    "Resectable" = "#74B9FF",
    "Borderline" = "#26A69A",
    "Locally Advanced" = "#F39C12",
    "Metastatic" = "#FF6B6B"
  )

  # Build per-sample data: stage, OS (months), cohort-wide score
  os_time_raw <- suppressWarnings(as.numeric(pheno$OS_Time))
  os_event_raw <- suppressWarnings(as.integer(pheno$OS_Censor))
  time_mo <- os_time_raw
  if (is.finite(max(time_mo, na.rm = TRUE)) && max(time_mo, na.rm = TRUE) > 400) {
    time_mo <- time_mo / 30.4375
  }
  dat_ws <- data.frame(
    geo_accession = as.character(pheno$geo_accession),
    stage_group = factor(stage_map(pheno$Disease_Stage), levels = names(pal_stage)),
    time_months = time_mo,
    OS_event = os_event_raw,
    score_custom = pheno$score_custom,
    stringsAsFactors = FALSE
  )
  dat_ws <- subset(dat_ws, !is.na(stage_group) & !is.na(time_months) & !is.na(OS_event) & is.finite(score_custom))

  plots_cohort <- list()
  for (st in names(pal_stage)) {
    dat_st <- subset(dat_ws, stage_group == st)
    if (nrow(dat_st) < 8) next
    # Median within this stage; split High/Low by that median
    med_stage <- median(dat_st$score_custom, na.rm = TRUE)
    dat_st$grp <- factor(
      ifelse(dat_st$score_custom >= med_stage, "High", "Low"),
      levels = c("Low", "High")
    )
    if (length(unique(dat_st$grp)) < 2) next

    fit <- survfit(Surv(time_months, OS_event) ~ grp, data = dat_st)
    lr <- survdiff(Surv(time_months, OS_event) ~ grp, data = dat_st)
    p_lr <- 1 - pchisq(lr$chisq, df = 1)
    cox <- coxph(Surv(time_months, OS_event) ~ grp, data = dat_st)
    hr <- exp(coef(cox))[1]
    ci <- as.vector(exp(confint(cox)))[1:2]
    label <- sprintf("p = %s\nHR = %.1f (%.1f-%.1f)",
                     format.pval(p_lr, digits = 2, eps = 0.001), hr, ci[1], ci[2])
    max_t <- max(dat_st$time_months, na.rm = TRUE)

    col_high <- pal_stage[[st]]
    col_low <- grDevices::adjustcolor(col_high, alpha.f = 0.55)

    p <- ggsurvplot(
      fit,
      data = dat_st,
      palette = unname(c(col_low, col_high)),
      risk.table = FALSE,
      legend.title = NULL,
      legend.labs = NULL,
      title = st,
      xlab = "Time (months)",
      ylab = "Overall survival probability",
      pval = label,
      pval.coord = c(max_t * 0.45, 0.95),
      pval.size = 4
    )
    p$plot <- p$plot + ggplot2::aes(linetype = grp) +
      ggplot2::scale_linetype_manual(
        values = c("Low" = "solid", "High" = "dashed"),
        name = "Score strata",
        labels = c("Low" = "Below median", "High" = "Above median"),
        guide = ggplot2::guide_legend(
          override.aes = list(color = c("grey70", "grey20"))
        )
      ) +
      ggplot2::guides(color = "none")
    plots_cohort[[st]] <- p$plot
  }

  if (length(plots_cohort) > 0) {
    ordered <- plots_cohort[intersect(names(pal_stage), names(plots_cohort))]
    print(wrap_plots(ordered, ncol = 2, guides = "collect") &
            ggplot2::theme(legend.position = "bottom"))
  }
}
```

![](custom-reference_files/figure-html/km-within-stage-cohort-median-1.png)
It does not appear that the cohort-wide signature provides robust
stratification of samples within each disease stage.

## 8. Stage-specific custom references and within-stage re-stratification

We now derive **stage-specific** survival signatures: for each disease
stage we fit a Cox model on that stage only, obtain a stage-specific
z-score reference, score only that stage’s samples with that reference,
and stratify by the **stage-specific** median (High/Low within that
stage). This tests whether a signature tailored to each stage can
further stratify survival within that disease stage

``` r
if (!requireNamespace("survminer", quietly = TRUE)) {
  message("Package 'survminer' not installed; skipping stage-specific KM plots.")
} else if (!requireNamespace("patchwork", quietly = TRUE)) {
  message("Package 'patchwork' not installed; skipping stage-specific KM plots.")
} else {
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(patchwork))

  plots <- list()
  stage_levels <- names(pal_stage)

  for (st in stage_levels) {
    dat_st <- subset(dat_stage, stage_group == st)
    if (nrow(dat_st) < 10) next

    # Align expression to these samples
    keep_ids <- intersect(rownames(bulk_mat), dat_st$geo_accession)
    if (length(keep_ids) < 10) next

    bulk_st <- bulk_mat[keep_ids, , drop = FALSE] # subset to stage samples (rows = samples if bulk_mat is samples×genes)
    pheno_st <- data.frame(
      sample_id = keep_ids,
      survival_time = dat_st$time_months[match(keep_ids, dat_st$geo_accession)],
      survival_event = dat_st$OS_event_raw[match(keep_ids, dat_st$geo_accession)],
      stringsAsFactors = FALSE
    )
    pheno_st <- pheno_st[!is.na(pheno_st$survival_time) & !is.na(pheno_st$survival_event), , drop = FALSE]
    bulk_st <- bulk_st[pheno_st$sample_id, , drop = FALSE]
    if (nrow(pheno_st) < 10) next

    ref_st <- suppressWarnings(
      derive_reference_from_bulk(
        bulk_expression = bulk_st,
        phenotype = pheno_st,
        sample_id_column = "sample_id",
        phenotype_type = "survival",
        survival_time = "survival_time",
        survival_event = "survival_event",
        verbose = FALSE
      )
    )

    # Score only within this stage using its own signature
    scores_st <- PhenoMap(
      expression = t(bulk_st), # PhenoMap expects genes × samples
      reference = ref_st,
      z_score_cutoff = 2,
      verbose = FALSE
    )
    col_sc <- grep("survival_z|weighted_sum", colnames(scores_st), value = TRUE)[1]
    sc <- scores_st[, col_sc]

    dat_km <- pheno_st
    dat_km$score_stage <- sc[match(dat_km$sample_id, rownames(scores_st))]
    dat_km <- dat_km[!is.na(dat_km$score_stage), , drop = FALSE]
    if (nrow(dat_km) < 10) next

    med <- median(dat_km$score_stage, na.rm = TRUE)
    dat_km$grp <- factor(ifelse(dat_km$score_stage >= med, "High", "Low"), levels = c("Low", "High"))
    fit <- survfit(Surv(survival_time, survival_event) ~ grp, data = dat_km)
    lr <- survdiff(Surv(survival_time, survival_event) ~ grp, data = dat_km)
    p_lr <- 1 - pchisq(lr$chisq, df = 1)
    cox <- coxph(Surv(survival_time, survival_event) ~ grp, data = dat_km)
    hr <- exp(coef(cox))[1]
    ci <- as.vector(exp(confint(cox)))[1:2]
    label <- sprintf("p = %s\nHR = %.1f (%.1f-%.1f)",
                     format.pval(p_lr, digits = 2, eps = 0.001), hr, ci[1], ci[2])
    max_t <- max(dat_km$survival_time, na.rm = TRUE)

    # Stage color for High; lighter shade for Low
    col_high <- pal_stage[[st]]
    col_low <- grDevices::adjustcolor(col_high, alpha.f = 0.55)

    p <- ggsurvplot(
      fit,
      data = dat_km,
      palette = unname(c(col_low, col_high)),
      risk.table = FALSE,
      legend.title = "Signature score",
      legend.labs = c("Low", "High"),
      title = st,
      xlab = "Time (months)",
      ylab = "Overall survival probability",
      pval = label,
      pval.coord = c(max_t * 0.35, 0.95),
      pval.size = 4
    )

    # Differentiate High vs Low with linetype so we can keep the stage color
    p$plot <- p$plot + ggplot2::aes(linetype = grp) +
      ggplot2::scale_linetype_manual(values = c("Low" = "solid", "High" = "dashed")) +
      ggplot2::guides(color = "none") + 
      theme(plot.title = element_text(hjust = 0.5))

    plots[[st]] <- p$plot
  }

  if (length(plots) == 0) {
    message("Not enough samples per stage with OS to derive and plot stage-specific signatures.")
  } else {
    ordered <- plots[intersect(stage_levels, names(plots))]
    p_surv_split <- wrap_plots(ordered, ncol = 4)
    p_surv_split
  }
}
```

![](custom-reference_files/figure-html/stage-specific-custom-1.png) It
appears that stage-specific **`PhenoMapR`** signature generation and
scoring significantly separates all but the metastatic patients.
Interestingly, it seems that the robustness of survival stratification
decreases with increasing disease stage severity, possibly suggesting
that applying **`PhenoMapR`** to traditionally lower risk populations
might identify patients who are actually higher risk.

## 9. Conclusions

In this vignette, we demonstrated how to generate a custom phenotype
signature using bulk expression and sample annotations using
**`PhenoMapR`**. The resulting signature successfully stratified a
cohort of 317 pancreatic cancer patients and did not enrich for specific
disease stage patients, suggesting that the **`PhenoMapR`** score
captured stage agnostic signatures. Additionally, we applied
stage-specific signature generation and scoring of patients using
**`PhenoMapR`** and show that the results further stratify all but the
metastatic stage patients.

## 10. References

**\[1\]** Nicolle, R. et al. Predictive genomic and transcriptomic
analysis on endoscopic ultrasound-guided fine needle aspiration
materials from primary pancreatic adenocarcinoma: A prospective
multicentre study. EBioMedicine 108, 105339. (2024).

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggsignif_0.6.4        survminer_0.5.2       ggpubr_0.6.3         
    ##  [4] survival_3.8-6        circlize_0.4.18       ComplexHeatmap_2.26.1
    ##  [7] patchwork_1.3.2       dplyr_1.2.1           ggplot2_4.0.2        
    ## [10] PhenoMapR_0.1.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1        farver_2.1.2            S7_0.2.1               
    ##  [4] fastmap_1.2.0           digest_0.6.39           lifecycle_1.0.5        
    ##  [7] cluster_2.1.8.2         magrittr_2.0.5          compiler_4.5.3         
    ## [10] rlang_1.2.0             sass_0.4.10             tools_4.5.3            
    ## [13] yaml_2.3.12             data.table_1.18.2.1     knitr_1.51             
    ## [16] labeling_0.4.3          htmlwidgets_1.6.4       curl_7.0.0             
    ## [19] splitstackshape_1.4.8.1 RColorBrewer_1.1-3      abind_1.4-8            
    ## [22] withr_3.0.2             purrr_1.2.2             BiocGenerics_0.56.0    
    ## [25] desc_1.4.3              stats4_4.5.3            googledrive_2.1.2      
    ## [28] colorspace_2.1-2        scales_1.4.0            iterators_1.0.14       
    ## [31] cli_3.6.6               rmarkdown_2.31          crayon_1.5.3           
    ## [34] ragg_1.5.2              generics_0.1.4          otel_0.2.0             
    ## [37] httr_1.4.8              rjson_0.2.23            cachem_1.1.0           
    ## [40] splines_4.5.3           parallel_4.5.3          matrixStats_1.5.0      
    ## [43] vctrs_0.7.3             Matrix_1.7-4            carData_3.0-6          
    ## [46] jsonlite_2.0.0          car_3.1-5               IRanges_2.44.0         
    ## [49] GetoptLong_1.1.1        S4Vectors_0.49.1-1      rstatix_0.7.3          
    ## [52] Formula_1.2-5           clue_0.3-68             systemfonts_1.3.2      
    ## [55] magick_2.9.1            foreach_1.5.2           jquerylib_0.1.4        
    ## [58] tidyr_1.3.2             glue_1.8.0              pkgdown_2.2.0          
    ## [61] codetools_0.2-20        gtable_0.3.6            shape_1.4.6.1          
    ## [64] tibble_3.3.1            pillar_1.11.1           htmltools_0.5.9        
    ## [67] R6_2.6.1                textshaping_1.0.5       doParallel_1.0.17      
    ## [70] evaluate_1.0.5          lattice_0.22-9          png_0.1-9              
    ## [73] backports_1.5.1         broom_1.0.12            gargle_1.6.1           
    ## [76] bslib_0.10.0            Rcpp_1.1.1              gridExtra_2.3          
    ## [79] HGNChelper_0.8.15       xfun_0.57               fs_2.0.1               
    ## [82] pkgconfig_2.0.3         GlobalOptions_0.1.4
