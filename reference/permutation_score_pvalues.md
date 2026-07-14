# Gene-shuffle permutation p-values for weighted-sum scores.

Gene-shuffle permutation p-values for weighted-sum scores.

## Usage

``` r
permutation_score_pvalues(
  expression_matrix,
  meta_z,
  n_perm = 100L,
  seed = 42L,
  pseudobulk = FALSE,
  verbose = TRUE
)
```

## Arguments

- expression_matrix:

  Genes-by-cells matrix aligned to `meta_z`.

- meta_z:

  Named numeric vector of reference z-scores.

- n_perm:

  Number of label shuffles (default 100).

- seed:

  Random seed.

- pseudobulk:

  Passed to
  [`compute_scores`](https://brooksbenard.github.io/PhenoMapR/reference/compute_scores.md).

- verbose:

  Print progress.

## Value

List with `observed` scores and `empirical_p` per cell.
