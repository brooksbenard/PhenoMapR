# Derive a pan-cohort meta-z reference from multiple bulk studies

Computes per-study gene z-scores with
[`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
and combines them with a Stouffer weighted meta-analysis (weights =
`sqrt(n_samples)` per study), analogous to PRECOG / TCGA meta-z
signatures.

## Usage

``` r
derive_meta_z_from_bulk_studies(
  studies,
  meta_z_label = "meta_z",
  hugo_species = c("human", "mouse"),
  binary_positive_reference = c("second", "first"),
  verbose = TRUE
)
```

## Arguments

- studies:

  A named list. Each element must be a list with at least
  `bulk_expression` and `phenotype`. Optional per-study fields are
  forwarded to
  [`derive_reference_from_bulk()`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
  (e.g. `phenotype_type`, `platform`, `probe_annotation`).

- meta_z_label:

  Character label for the output z-score column (default `"meta_z"`).

- hugo_species:

  Species for HGNC symbol validation (`"human"` or `"mouse"`); passed to
  each per-study call.

- binary_positive_reference:

  For binary phenotypes, which level is the positive reference
  (`"second"` or `"first"`); passed to each per-study call.

- verbose:

  Logical.

## Value

A data.frame of combined meta-z scores suitable for
[`PhenoMap`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md).
