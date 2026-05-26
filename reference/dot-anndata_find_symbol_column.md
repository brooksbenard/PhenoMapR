# Find a HUGO-symbol-style column in AnnData.var that overlaps the reference

Scans `AnnData.var` for the most common gene-symbol column conventions
(`gene_symbol`, `gene_symbols`, `feature_name`, `Symbol`,
`gene_short_name`, `hgnc_symbol`, `gene_name`) and returns the first
column with any overlap against `gene_subset`.

## Usage

``` r
.anndata_find_symbol_column(obj, gene_subset)
```

## Value

A list with `column`, `var_names` (the AnnData `var_names` whose symbol
matched, suitable for passing to `.anndata_subset_var`) and `symbols`
(the matched HUGO symbols in the same order). Returns `NULL` when no
usable column is found.
