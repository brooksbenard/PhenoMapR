# Regress technical covariates and return scaled expression for scoring

Uses Seurat normalization, optional cell-cycle scoring, and `ScaleData`
with `vars.to.regress`. Intended for activity-adjusted PhenoMap scoring.

## Usage

``` r
regress_expression_for_scoring(
  expression_matrix,
  signature_genes,
  cell_metadata = NULL,
  cell_id_column = "cell_id",
  vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
  use_counts = TRUE,
  verbose = TRUE
)
```

## Arguments

- expression_matrix:

  Genes-by-cells count or normalized matrix.

- signature_genes:

  Gene symbols to scale (typically reference signature genes).

- cell_metadata:

  Optional data.frame with cell covariates; must include a column
  matching `colnames(expression_matrix)` via `cell_id_column`.

- cell_id_column:

  Metadata column with cell barcodes (default `"cell_id"`).

- vars_to_regress:

  Character vector of metadata columns to regress out.

- use_counts:

  Logical. If `TRUE` (default), treat `expression_matrix` as counts and
  run `NormalizeData`; if `FALSE`, use the matrix as normalized input
  without re-normalizing.

- verbose:

  Print progress messages.

## Value

A genes-by-cells matrix from the Seurat `scale.data` layer.
