# scIMPEL

**S**ingle-**C**ell **IMP**act of **E**xpression on **L**ong-term outcomes — weighted sum prognostic scoring for bulk, single-cell, and spatial transcriptomics using gene expression and prognostic z-scores from PRECOG, TCGA, Pediatric, and ICI reference datasets.

## Installation

Install directly from GitHub (requires the `remotes` package):

```r
# Install from GitHub
remotes::install_github("bbenard/scIMPEL")
```

Or with devtools:

```r
devtools::install_github("bbenard/scIMPEL")
```

Dependencies (e.g. `dplyr`, `Matrix`, `glue`, `progress`) will be installed automatically. For Seurat/SCE support, install suggested packages as needed.

## Quick Start

```r
library(scIMPEL)

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

## Features

- **Multiple Input Formats**: Supports matrices, data.frames, Seurat objects, SingleCellExperiment, SpatialExperiment, and AnnData
- **Multiple References**: PRECOG, TCGA, Pediatric, and ICI prognostic datasets
- **Flexible Scoring**: Cell-level or pseudobulk scoring
- **Custom Signatures**: Use your own z-score references
- **Efficient**: Vectorized operations for fast scoring

## Supported Input Types

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

## Reference Datasets

### PRECOG
Pan-cancer prognostic meta-analysis from Gentles et al. (2015).
```r
# List available cancer types
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
Response predictions for immunotherapy patients.
```r
list_cancer_types("ici_precog")

# Format: "CANCER" or "CANCER_Metastatic"
# Examples: "MELANOMA", "MELANOMA_Metastatic", "NSCLC", etc.
```

## Advanced Usage

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

### Automatic Dataset Mapping
If you have the `datasets_info` table configured:
```r
scores <- score_expression(
  expression = my_seurat,
  reference = "precog",
  use_dataset_info = TRUE,
  dataset = "my_dataset_name"
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

## Utility Functions

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

## Parameters

### `score_expression()`

- **expression**: Expression data (matrix, Seurat, SCE, etc.)
- **reference**: Reference dataset name or custom data.frame
- **cancer_type**: Cancer type label (required for built-in references)
- **z_score_cutoff**: Absolute z-score threshold (default: 2)
- **pseudobulk**: Aggregate before scoring? (default: FALSE)
- **group_by**: Grouping variable for pseudobulk
- **assay**: Assay name for Seurat/SCE objects
- **slot**: Seurat slot ("data", "counts", "scale.data")
- **use_dataset_info**: Use automatic cancer type mapping
- **dataset**: Dataset name for mapping
- **verbose**: Print progress messages

## Output Format

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

## Citation

If you use this package, please cite the relevant reference datasets:

**PRECOG**: Gentles et al. (2015). The prognostic landscape of genes and infiltrating immune cells across human cancers. *Cell* 163(5):1193-1205.

**TCGA**: The Cancer Genome Atlas Research Network

## License

MIT

## Contributing

Issues and pull requests welcome at https://github.com/bbenard/scIMPEL
