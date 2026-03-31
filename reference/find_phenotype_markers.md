# Find Unique Marker Genes for Adverse and Favorable Phenotype Groups

Performs differential expression between the top (adverse) and bottom
(favorable) phenotype groups versus the rest. For `Seurat` input, uses
[`Seurat::FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html).
For matrix, `Matrix`, or `SingleCellExperiment` input, uses a
Wilcoxon-based path (presto if available, else base R) and does not
require Seurat.

## Usage

``` r
find_phenotype_markers(
  expression,
  group_labels,
  group_column = NULL,
  cell_id_column = "cell_id",
  marker_scope = c("phenotype_groups", "cell_type_specific"),
  cell_type_column = NULL,
  assay = NULL,
  slot = "data",
  test.use = c("wilcox", "t", "roc", "negbinom", "poisson", "LR"),
  min.pct = 0.1,
  logfc.threshold = 0.25,
  pval_threshold = 0.05,
  verbose = TRUE,
  max_cells_per_ident = 5000L,
  validate_expression_axes = TRUE,
  celltype_unique_genes = TRUE,
  celltype_contrast = c("within_cell_type", "vs_cohort_rest"),
  ...
)
```

## Arguments

- expression:

  Expression matrix (genes x cells), a `Matrix` (e.g. `dgCMatrix`), a
  Seurat object, or a SingleCellExperiment. Must match cells in
  `group_labels`.

- group_labels:

  Either a character vector of group labels (`"Most Adverse"`,
  `"Most Favorable"`, `"Other"`) in the same order as columns of
  `expression`, or a data.frame from
  [`define_phenotype_groups()`](https://brooksbenard.github.io/PhenoMapR/reference/define_phenotype_groups.md)
  (see `group_column`).

- group_column:

  If `group_labels` is a data.frame, the name of the column containing
  `"Most Adverse"` / `"Most Favorable"` / `"Other"`.

- cell_id_column:

  If `group_labels` is a data.frame, the column name for cell/sample IDs
  (default `"cell_id"`).

- marker_scope:

  Either:

  - `"phenotype_groups"`: find markers for the adverse and favorable
    phenotype groups globally (cell type agnostic; default).

  - `"cell_type_specific"`: for each cell type, find markers for cells
    of that type in the adverse (or favorable) tail vs a reference group
    controlled by `celltype_contrast` (default: same cell type, not in
    that tail; use `"vs_cohort_rest"` for the original cohort-wide
    contrast).

- cell_type_column:

  When `marker_scope = "cell_type_specific"`, the column in
  `group_labels` that contains cell type labels.

- assay:

  Assay name for Seurat/SCE (e.g. `"RNA"`).

- slot:

  Layer for Seurat: `"data"`, `"counts"`, or `"scale.data"` (default
  `"data"`).

- test.use:

  Seurat `test.use`: `"wilcox"` (default), `"t"`, `"roc"`, `"negbinom"`,
  `"poisson"`, or `"LR"`.

- min.pct:

  Minimum fraction of cells expressing the gene in either group (default
  0.1). Passed to `FindMarkers`.

- logfc.threshold:

  Minimum absolute log2 fold change (default 0.25). Passed to
  `FindMarkers`.

- pval_threshold:

  Maximum unadjusted p-value to include (default 0.05).

- verbose:

  Print progress messages (default TRUE).

- max_cells_per_ident:

  When any phenotype group exceeds this many cells, subsample to this
  limit before FindMarkers (default 5000). Reduces memory for large
  objects. Set to `Inf` to disable.

- validate_expression_axes:

  If TRUE (default), check that the expression matrix is genes
  \\\times\\ cells (transpose when a clear samples \\\times\\ genes
  layout is detected) and normalize colnames for ID matching. Use
  log-normalized (e.g. log1p CPM) `data` for matrix input when possible.

- celltype_unique_genes:

  When `marker_scope = "cell_type_specific"` and `TRUE` (default), each
  gene is kept in at most one row across adverse and favorable marker
  tables: the row with the largest `avg_log2FC` (tie-break by `p_adj`).
  This avoids the same gene (e.g. highly expressed housekeeping
  transcripts) appearing as a top marker in many cell-type blocks. Set
  to `FALSE` to retain all significant hits per block.

- celltype_contrast:

  When `marker_scope = "cell_type_specific"`: `"within_cell_type"`
  (default) compares each phenotype tail within a cell type to other
  phenotype groups in the **same** cell type only. `"vs_cohort_rest"`
  restores the original behaviour: (cell type \\\cap\\ tail) vs **all
  other cells** in the dataset (other cell types and other phenotype
  groups).

- ...:

  Additional arguments passed to
  [`Seurat::FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  when input is a Seurat object (ignored for matrix/SCE/Matrix input).

## Value

A list with:

- `adverse_markers`: data.frame of genes that are markers of the adverse
  (top score) group (vs rest).

- `favorable_markers`: data.frame of genes that are markers of the
  favorable (bottom score) group (vs rest).

Each data.frame has columns: `gene`, `avg_log2FC`, `pct_in_group`,
`pct_rest`, `p_val`, `p_adj`. When
`marker_scope = "cell_type_specific"`, each row is one gene for one cell
type; the `cell_type` column identifies which type the in-group was
anchored on. If `celltype_contrast = "within_cell_type"` (default), the
reference is other phenotype groups within that cell type; if
`"vs_cohort_rest"`, the reference is all other cells with a phenotype
label (original behaviour).

## Details

Phenotype tails (e.g. top/bottom 5\\ For
`marker_scope = "cell_type_specific"`, the contrast is set by
`celltype_contrast` (see above). `group_labels` cell IDs must match
`colnames(expression)` (after trimming); use the same identifiers you
used when building the score table and phenotype groups.

## See also

[`plot_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/plot_phenotype_markers.md)
for heatmaps of marker expression.

## Examples

``` r
if (FALSE) { # \dontrun{
scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
groups <- define_phenotype_groups(scores, percentile = 0.05)
markers <- find_phenotype_markers(
  seurat_obj,
  group_labels = groups,
  group_column = "phenotype_group_weighted_sum_score_precog_BRCA",
  cell_id_column = "cell_id"
)
head(markers$adverse_markers)
head(markers$favorable_markers)
} # }
```
