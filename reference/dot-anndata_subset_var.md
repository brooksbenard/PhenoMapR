# Subset an AnnData on the var (gene) axis in Python

Uses a tiny Python helper defined in `__main__` to mask AnnData on
`var_names` and return the subset's `.X` plus the new `var_names` /
`obs_names`. All return values are kept as Python objects so the caller
can apply the most memory-efficient conversion.

## Usage

``` r
.anndata_subset_var(obj, gene_keep)
```
