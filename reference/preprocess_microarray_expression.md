# PRECOG-style microarray preprocessing (genes x samples)

PRECOG-style microarray preprocessing (genes x samples)

## Usage

``` r
preprocess_microarray_expression(
  expr_genes_by_sample,
  min_gene_sd = 1e-05,
  max_missing_frac = 0.8,
  verbose = TRUE
)
```

## Arguments

- expr_genes_by_sample:

  Numeric matrix, genes as rows, samples as columns.

- min_gene_sd:

  Minimum per-gene standard deviation to retain a gene.

- max_missing_frac:

  Maximum allowed missing fraction per gene or sample.

- verbose:

  Logical.

## Value

Preprocessed matrix (genes x samples).
