# Convert AnnData.X into a genes × cells Matrix

Returns a `dgCMatrix` (sparse) when `adata.X` is a scipy sparse matrix,
or a regular dense matrix otherwise. Always genes × cells (i.e. the
transpose of AnnData's native cells × genes layout) with
`rownames = var_names` and `colnames = obs_names`.

## Usage

``` r
.anndata_X_to_genes_cells(obj, gene_subset = NULL)
```

## Details

Two memory optimisations versus the naive `as.matrix(adata$X)` approach:

1.  When `gene_subset` is supplied, the AnnData is sliced on the var
    axis *in Python* (`adata[:, mask].X`) before any data is copied
    into R. This is the dominant cost reduction for multi-GB AnnData
    objects: only the genes that actually contribute to the score are
    transferred.

2.  For scipy-sparse `.X`, the AnnData native CSR storage of a (n_obs ×
    n_vars) matrix is the same as the CSC storage of the transposed
    (n_vars × n_obs) matrix. We reuse `indices`, `indptr` and `data`
    arrays directly to build a `dgCMatrix` in genes × cells orientation,
    with no extra allocation for a transpose pass.
