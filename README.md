# PhenoMap <a href='https://brooksbenard.github.io/PhenoMap'><img src='inst/figures/PhenoMap_logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/brooksbenard/PhenoMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brooksbenard/PhenoMap/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

PhenoMap is a semi-supervised method to map phenotypes associated with bulk gene expression data onto bulk, single cell, and spatial transcriptomics data. PhenoMap nominates and rank-orders samples, cells, and spatial locations most associated with gene expression signatures correlated with a phenotype of interest (e.g. overall survival).

![PhenoMap schematic](inst/figures/PhenoMap_schematic.png)

<details>
<summary><b>Installation</b></summary>

Install directly from GitHub (requires the `remotes` package):

```r
# Install from GitHub
remotes::install_github("brooksbenard/PhenoMap")
```

Or with devtools:

```r
devtools::install_github("brooksbenard/PhenoMap")
```

Dependencies (e.g. `dplyr`, `Matrix`, `glue`, `progress`) will be installed automatically. For Seurat/SCE support, install suggested packages as needed.

</details>

<details>
<summary><b>Quick Start</b></summary>

```r
library(PhenoMap)

# Score a bulk expression matrix
scores <- score_expression(
  expression = bulk_matrix,  # genes x samples
  reference = "precog",
  cancer_type = "BRCA"
)

# Score single-cell data
scores <- score_expression(
  expression = seurat_obj,
  reference = "tcga",
  cancer_type = "LUAD",
  assay = "RNA",
  slot = "data"
)

# Score with pseudobulk aggregation
scores <- score_expression(
  expression = seurat_obj,
  reference = "ici_precog",
  cancer_type = "MELANOMA_Metastatic",
  pseudobulk = TRUE,
  group_by = "patient_id"
)
```

</details>

## Features

- **Multiple Input Formats**: Supports matrices, data.frames, Seurat objects, SingleCellExperiment, SpatialExperiment, and AnnData
- **Multiple References**: PRECOG, TCGA, Pediatric, and ICI prognostic datasets
- **Flexible Scoring**: Cell-level or pseudobulk scoring
- **Custom Signatures**: Use your own z-score references
- **Efficient**: Vectorized operations for fast scoring

<details>
<summary><b>Supported Input Types</b></summary>

### 1. Matrix/Data.frame
```r
# Expression matrix: genes (rows) x samples/cells (columns)
expression_matrix <- matrix(...)
rownames(expression_matrix) <- gene_names
colnames(expression_matrix) <- cell_names

scores <- score_expression(expression_matrix, reference = "precog", cancer_type = "BRCA")
```

### 2. Seurat Objects
```r
library(Seurat)

# Single-cell
scores <- score_expression(
  seurat_obj,
  reference = "tcga",
  cancer_type = "LUAD",
  assay = "RNA",
  slot = "data"
)

# Add scores back to Seurat object
seurat_obj <- add_scores_to_seurat(seurat_obj, scores)

# Spatial
scores <- score_expression(
  spatial_seurat,
  reference = "precog",
  cancer_type = "BRCA",
  assay = "Spatial",
  slot = "counts"
)
```

### 3. SingleCellExperiment
```r
library(SingleCellExperiment)

scores <- score_expression(
  sce_obj,
  reference = "pediatric_precog",
  cancer_type = "Neuroblastoma",
  assay = "logcounts"
)

# Add scores to colData
sce_obj <- add_scores_to_sce(sce_obj, scores)
```

### 4. AnnData (Python)
```r
library(reticulate)

adata <- import("scanpy")$read_h5ad("data.h5ad")

scores <- score_expression(
  adata,
  reference = "precog",
  cancer_type = "BRCA"
)
```

</details>

<details>
<summary><b>Reference Datasets</b></summary>

Prognostic meta-z scores and cancer-type labels in PhenoMap are sourced from **PRECOG 2.0** ([Benard et al., *Nucleic Acids Research* 2026](https://academic.oup.com/nar/article/54/D1/D1579/8324954)). Additional citations for the underlying data and methods:

- **PRECOG / TCGA (pan-cancer meta-z and TCGA z-scores)**: [Gentles et al., *Nature Medicine* 2015](https://www.nature.com/articles/nm.3909) — The prognostic landscape of genes and infiltrating immune cells across human cancers.
- **Pediatric PRECOG**: [Stahl et al., *Cancers* 2021](https://www.mdpi.com/2072-6694/13/4/854).

### PRECOG
Pan-cancer prognostic meta-analysis (meta-z scores from PRECOG 2.0 / Gentles et al.).
```r
list_cancer_types("precog")
# Examples: "BRCA", "LUAD", "COAD", "PRAD", etc.
```

### TCGA
TCGA survival analysis z-scores.
```r
list_cancer_types("tcga")
# Examples: "BRCA", "LUAD", "UCEC", "KIRC", etc.
```

### Pediatric
Pediatric cancer prognostic signatures.
```r
list_cancer_types("pediatric_precog")
# Examples: "Neuroblastoma", "Medulloblastoma", etc.
```

### ICI (Immune Checkpoint Inhibitor)
Response-related signatures for immunotherapy patients.
```r
list_cancer_types("ici_precog")
# Format: "CANCER" or "CANCER_Metastatic"
# Examples: "MELANOMA", "MELANOMA_Metastatic", "NSCLC", etc.
```

</details>

<details>
<summary><b>Vignette: PAAD single-cell analysis</b></summary>

This example walks through loading the included **PAAD (pancreatic adenocarcinoma) GSE111672** Seurat object, scoring cells with a prognostic reference, and inspecting results. The file `PAAD_GSE111672_seurat.rds` is provided in the repository.

```r
library(PhenoMap)
library(Seurat)

# Load the PAAD single-cell dataset (place PAAD_GSE111672_seurat.rds in your working directory)
seurat <- readRDS("PAAD_GSE111672_seurat.rds")

# Score each cell with PRECOG pancreatic adenocarcinoma prognostic z-scores
scores <- score_expression(
  expression = seurat,
  reference = "precog",
  cancer_type = "PAAD",
  assay = "RNA",
  slot = "data",
  verbose = TRUE
)

# Attach scores to Seurat metadata for downstream plotting and subsetting
seurat <- add_scores_to_seurat(seurat, scores)

# Quick look at score distribution
head(scores)
summary(scores[, 1])

# Optional: plot score distribution
plot_score_distribution(scores, main = "PAAD GSE111672 — PRECOG prognostic score")
```

**Defining prognostic groups and marker genes (optional)**  
Identify cells in the top and bottom 5% by score (adverse vs. favorable prognosis) and find marker genes with Seurat:

```r
# Define adverse (top 5%) and favorable (bottom 5%) prognostic groups
groups <- define_prognostic_groups(scores, percentile = 0.05)

# Inspect group counts
score_col <- names(scores)[1]
table(groups[[paste0("prognostic_group_", score_col)]])

# Find marker genes (requires Seurat)
markers <- find_prognostic_markers(
  seurat,
  group_labels = groups,
  group_column = paste0("prognostic_group_", score_col),
  cell_id_column = "cell_id"
)
head(markers$adverse_markers)
head(markers$favorable_markers)
```

Interpretation: higher scores correspond to worse prognosis (adverse); lower scores to better prognosis (favorable). The derived groups and markers can be used for visualization (e.g. UMAP colored by score or group) and biological interpretation.

</details>

<details>
<summary><b>Advanced Usage</b></summary>

### Custom Reference Data
```r
# Create custom z-score reference
custom_ref <- data.frame(
  row.names = c("TP53", "MYC", "EGFR", "BRCA1"),
  my_signature = c(3.2, -2.5, 2.8, 1.9)
)

scores <- score_expression(
  expression = my_data,
  reference = custom_ref,
  z_score_cutoff = 1.5
)
```

### Derive reference from bulk expression and phenotype
If you have bulk expression (samples × genes) and a phenotype (e.g. response R/NR, survival, or continuous), you can derive gene-level z-scores and use them as the reference for scoring single-cell or spatial data. The function cleans gene names to HUGO symbols, checks/normalizes expression, then computes association z-scores (Cox for survival, logistic regression for binary, correlation for continuous).
```r
# Bulk: samples in rows, genes in columns
bulk_expr <- matrix(...)   # e.g. 50 samples × 5000 genes
pheno <- data.frame(
  sample_id = rownames(bulk_expr),
  response = c("R", "NR", ...)  # or time + event for survival
)

# Derive reference z-scores
ref <- derive_reference_from_bulk(
  bulk_expression = bulk_expr,
  phenotype = pheno,
  sample_id_column = "sample_id",
  phenotype_column = "response",
  phenotype_type = "binary"   # or "survival", "continuous", "auto"
)

# Score single-cell or spatial data with the derived reference
scores <- score_expression(expression = my_seurat, reference = ref)
```

For survival, provide `survival_time` and `survival_event` column names and set `phenotype_type = "survival"`. Optional: install `HGNChelper` for HUGO symbol cleaning and `survival` for Cox models.

### Pseudobulk Aggregation
```r
# Aggregate single cells by patient before scoring
scores <- score_expression(
  seurat_obj,
  reference = "tcga",
  cancer_type = "LUAD",
  pseudobulk = TRUE,
  group_by = "patient_id"
)
```

### Normalize Scores
```r
scores_df <- score_expression(...)

# Convert to z-scores
scores_df$score_zscore <- normalize_scores(scores_df[, 1])
```

### Prognostic Groups and Marker Genes
Define the top and bottom 5% prognostic cells per dataset (adverse = highest scores, favorable = lowest) and find unique marker genes using **Seurat's FindMarkers** (requires Seurat):
```r
# Score and define groups (top 5% = adverse, bottom 5% = favorable)
scores <- score_expression(seurat_obj, reference = "precog", cancer_type = "BRCA")
groups <- define_prognostic_groups(scores, percentile = 0.05)

# One group column per score; values: "adverse", "favorable", "middle"
table(groups$prognostic_group_weighted_sum_score_precog_BRCA)

# Find marker genes via Seurat::FindMarkers (adverse vs rest, favorable vs rest)
markers <- find_prognostic_markers(
  seurat_obj,
  group_labels = groups,
  group_column = "prognostic_group_weighted_sum_score_precog_BRCA",
  cell_id_column = "cell_id"
)
head(markers$adverse_markers)   # genes enriched in top 5% (worst prognosis)
head(markers$favorable_markers) # genes enriched in bottom 5% (best prognosis)
```

</details>

<details>
<summary><b>Utility Functions</b></summary>

### Check Gene Coverage
```r
my_genes <- rownames(expression_matrix)
coverage <- get_gene_coverage(my_genes)
print(coverage)
```

### Get Top Prognostic Genes
```r
# Top 100 genes for breast cancer
top_genes <- get_top_prognostic_genes(
  reference = "precog",
  cancer_type = "BRCA",
  n = 100,
  direction = "both"  # or "positive" or "negative"
)
```

### Visualize Scores
```r
plot_score_distribution(scores, main = "BRCA Prognostic Scores")
```

</details>

<details>
<summary><b>Parameters</b></summary>

### `score_expression()`

- **expression**: Expression data (matrix, Seurat, SCE, etc.)
- **reference**: Reference dataset name or custom data.frame
- **cancer_type**: Cancer type label (required for built-in references)
- **z_score_cutoff**: Absolute z-score threshold (default: 2)
- **pseudobulk**: Aggregate before scoring? (default: FALSE)
- **group_by**: Grouping variable for pseudobulk
- **assay**: Assay name for Seurat/SCE objects
- **slot**: Seurat slot ("data", "counts", "scale.data")
- **verbose**: Print progress messages

</details>

<details>
<summary><b>Output Format</b></summary>

Returns a data.frame with:
- Rows: samples/cells
- Columns: scores named `weighted_sum_score_{reference}_{cancer_type}`

```r
scores <- score_expression(...)

head(scores)
#                    weighted_sum_score_precog_BRCA
# Cell_1                                   234.5
# Cell_2                                   189.3
# Cell_3                                   -42.1
```

</details>

<details>
<summary><b>Citation</b></summary>

If you use PhenoMap, please cite the package and the reference datasets:

- **PhenoMap / PRECOG 2.0 (meta-z scores)**: Benard B et al. PRECOG 2.0: an updated resource of pan-cancer gene-level prognostic meta-z scores. *Nucleic Acids Research* (2026). [https://academic.oup.com/nar/article/54/D1/D1579/8324954](https://academic.oup.com/nar/article/54/D1/D1579/8324954)
- **PRECOG / TCGA**: Gentles AJ et al. The prognostic landscape of genes and infiltrating immune cells across human cancers. *Nature Medicine* 21, 938–945 (2015). [https://www.nature.com/articles/nm.3909](https://www.nature.com/articles/nm.3909)
- **Pediatric PRECOG**: Stahl et al. *Cancers* 13(4), 854 (2021). [https://www.mdpi.com/2072-6694/13/4/854](https://www.mdpi.com/2072-6694/13/4/854)

</details>

<details>
<summary><b>Session info (example)</b></summary>

Reproducibility summary from a typical session after `library(PhenoMap)`:

```r
sessionInfo()
```

Example output:

```
R version 4.4.x (YYYY-MM-DD)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: ...

Matrix products: default
BLAS:   ...
LAPACK: ...

locale:
 [1] LC_CTYPE=en_US.UTF-8    LC_NUMERIC=C            LC_TIME=en_US.UTF-8
 [4] LC_COLLATE=en_US.UTF-8  LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8    LC_NAME=C                LC_ADDRESS=C
[10] LC_TELEPHONE=C         LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] PhenoMap_0.1.0   dplyr_xxx         Matrix_xxx        glue_xxx
[5] progress_xxx

loaded via a namespace (and not attached):
[1] ...
```

To capture your own session info (e.g. for issues or reports):

```r
# Minimal
sessionInfo()

# With devtools (packages and versions)
devtools::session_info()
```

</details>

## License

MIT

## Contributing

Issues and pull requests welcome at https://github.com/brooksbenard/PhenoMap
