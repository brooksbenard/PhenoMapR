# Process AnnData Object

Convert a Python anndata object into a genes × cells expression matrix
that the rest of PhenoMapR understands. Optimized for very large
objects:

- When `genes_to_extract` is supplied (e.g. the reference genes that
  pass the z-score cutoff in
  [`PhenoMap()`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md)),
  the AnnData is subset *on the Python side* before any data is copied
  into R. This is the single biggest memory win for multi-GB `.h5ad`
  files: only a few hundred to a few thousand genes typically pass the
  cutoff, so we transfer ~1-3\\ it.

- scipy-sparse `.X` is reinterpreted directly as a `dgCMatrix` in genes
  × cells orientation by treating the native CSR-of-(cells×genes)
  storage as CSC-of-(genes×cells); this avoids `Matrix::t()` and the
  doubling of memory it would otherwise cause.

## Usage

``` r
process_anndata(obj, pseudobulk, group_by, genes_to_extract = NULL)
```
