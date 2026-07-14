# Score Expression Data with Prognostic Gene Signatures

Calculate weighted sum scores for expression data using prognostic
z-scores from various reference datasets (PRECOG, TCGA, Pediatric, ICI).

## Usage

``` r
PhenoMap(
  expression,
  reference,
  cancer_type = NULL,
  z_score_cutoff = 2,
  pseudobulk = FALSE,
  group_by = NULL,
  assay = NULL,
  slot = "data",
  verbose = TRUE,
  reference_sign = 1L,
  score_mode = c("weighted_sum", "activity_adjusted"),
  vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
  sample_metadata = NULL,
  cell_id_column = "cell_id",
  permutation_n = NULL,
  permutation_seed = 42L
)
```

## Arguments

- expression:

  Expression data. Can be:

  - Matrix or data.frame (genes x samples/cells)

  - Seurat object

  - SingleCellExperiment object

  - SpatialExperiment object

  - AnnData object (via reticulate)

- reference:

  Reference dataset name ("precog", "tcga", "pediatric_precog",
  "ici_precog") or a custom data.frame with genes as rownames and
  z-scores as values

- cancer_type:

  Cancer type label. Required if using built-in reference datasets.
  Should match column names in reference data. For ICI PRECOG, use
  format "ABBREV" or "ABBREV_Metastatic" (e.g., "MELANOMA",
  "MELANOMA_Metastatic")

- z_score_cutoff:

  Absolute z-score threshold for filtering genes (default: 2)

- pseudobulk:

  Logical. If TRUE, aggregate expression before scoring (default: FALSE)

- group_by:

  Column name for pseudobulk grouping. Required if pseudobulk = TRUE.
  For Seurat/SCE objects, should be a metadata column.

- assay:

  Assay name for Seurat/SCE objects (default: "RNA" for sc, "Spatial"
  for spatial)

- slot:

  Data layer for Seurat objects: "data" for normalized, "counts" for
  raw, "scale.data" for scaled (default: "data"). In Seurat v5+ this
  maps to the layer parameter; in Seurat v4, slot. The function handles
  both automatically.

- verbose:

  Logical. Print progress messages (default: TRUE)

- reference_sign:

  Multiply all reference gene z-scores by this value before the weighted
  sum. Use `1` (default) or `-1` to flip the sign of resulting PhenoMap
  scores without changing the reference matrix (e.g. if you want
  positive scores to mean the opposite association). For custom
  references from
  [`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md),
  prefer setting `binary_positive_reference` so the reference itself
  matches your factor levels (first level = positive association when
  `"first"`).

- score_mode:

  Scoring mode: `"weighted_sum"` (default) uses the standard dot product
  of expression and reference z-scores; `"activity_adjusted"` regresses
  technical covariates per gene (via Seurat `ScaleData`) before
  computing the weighted sum on the regressed layer. Requires Seurat for
  matrix inputs.

- vars_to_regress:

  Metadata columns to regress when `score_mode = "activity_adjusted"`.
  Default: `c("S.Score", "G2M.Score", "nCount_RNA")`. For matrix input,
  pass covariates via `sample_metadata` (see below). `nCount_RNA` is
  computed from column sums when absent.

- sample_metadata:

  Optional data.frame of per-cell/sample covariates for
  `score_mode = "activity_adjusted"` when `expression` is a matrix. Must
  include a cell identifier column (default `"cell_id"`).

- cell_id_column:

  Column in `sample_metadata` matching matrix colnames (default
  `"cell_id"`).

- permutation_n:

  If a positive integer, compute gene-shuffle permutation null scores
  and add `empirical_p_*` columns alongside raw score columns.

- permutation_seed:

  Random seed for permutation null (default `42`).

## Value

A data.frame with samples/cells as rows and score columns. Column names
follow pattern: `PhenoMapR_{reference}_{cancer_type}`. **Directionality
(built-in PRECOG/TCGA/ICI references)**: higher score = worse prognosis
(adverse); lower score = better prognosis (favorable), matching positive
reference z = worse survival. **Custom references** from
[`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
follow the sign convention you chose there (see
`binary_positive_reference`) times `reference_sign`. For Seurat/SCE
objects, scoring may use only reference genes internally for memory
efficiency. Always add scores to the **same (full) object** using
`add_scores_to_seurat` or `add_scores_to_sce` so that all genes are
retained for downstream analyses (e.g. cell type marker genes).

## Examples

``` r
if (FALSE) { # \dontrun{
# Bulk expression matrix
scores <- PhenoMap(
  expression = bulk_matrix,
  reference = "precog",
  cancer_type = "BRCA"
)

# Single cell Seurat object
scores <- PhenoMap(
  expression = seurat_obj,
  reference = "tcga",
  cancer_type = "LUAD",
  assay = "RNA",
  slot = "data"
)

# Spatial with pseudobulk
scores <- PhenoMap(
  expression = spatial_seurat,
  reference = "ici_precog",
  cancer_type = "MELANOMA_Metastatic",
  pseudobulk = TRUE,
  group_by = "sample_id"
)

# Custom reference data
custom_ref <- data.frame(
  row.names = c("TP53", "MYC", "EGFR"),
  my_signature = c(3.2, -2.5, 2.8)
)
scores <- PhenoMap(
  expression = my_data,
  reference = custom_ref
)
} # }
```
