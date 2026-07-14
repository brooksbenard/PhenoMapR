# Define Phenotype Groups (Top and Bottom Percentile)

For each score column (dataset), labels cells as adverse (top
percentile, highest scores), favorable (bottom percentile, lowest
scores), or middle.

## Usage

``` r
define_phenotype_groups(scores, percentile = 0.05, score_columns = NULL)
```

## Arguments

- scores:

  Data.frame of phenotype scores from `PhenoMap`. Rows = cells/samples,
  columns = score variables (e.g. PhenoMapR_precog_BRCA).

- percentile:

  Fraction of cells in each tail (default 0.05 for top and bottom 5%).

- score_columns:

  Character vector of score column names to use. If NULL, all numeric
  columns in `scores` are used.

## Value

A data.frame with:

- `cell_id`: same as row names of `scores`

- For each score column: a `phenotype_group_<name>` column with values
  `"Most Adverse"` (top percentile), `"Most Favorable"` (bottom
  percentile), or `"Other"`

## Examples

``` r
if (FALSE) { # \dontrun{
scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
groups <- define_phenotype_groups(scores, percentile = 0.05)
table(groups$phenotype_group_PhenoMapR_precog_BRCA)
} # }
```
