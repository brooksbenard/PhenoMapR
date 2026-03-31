# Run marker detection on a matrix by cell type

For each cell type, compute markers for adverse and favorable tails; the
reference group is set by `celltype_contrast` (see
[`find_phenotype_markers`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)).

## Usage

``` r
run_markers_on_matrix_by_celltype(
  mat,
  group_vec,
  cell_type_vec,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  pval_threshold = 0.05,
  max_cells_per_ident = 5000L,
  verbose = TRUE,
  unique_genes_across_celltypes = TRUE,
  celltype_contrast = c("within_cell_type", "vs_cohort_rest")
)
```
