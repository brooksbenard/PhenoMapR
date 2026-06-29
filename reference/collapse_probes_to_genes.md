# Collapse probe-level microarray expression to gene symbols

Maps platform probe IDs to gene symbols using a GPL-style annotation
table, then averages expression for probes that map to the same symbol.
Follows the probe summarization strategy used in the original PRECOG
resource.

## Usage

``` r
collapse_probes_to_genes(expr_probe_by_sample, annot_df)
```

## Arguments

- expr_probe_by_sample:

  Matrix with probes as rows and samples as columns.

- annot_df:

  Annotation table (e.g. from GEO GPL) with probe and gene columns.

## Value

Matrix with gene symbols as row names and the same sample columns.
