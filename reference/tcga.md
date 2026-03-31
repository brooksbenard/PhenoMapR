# TCGA Prognostic Z-Scores

Gene-level prognostic z-scores derived from The Cancer Genome Atlas
(TCGA) survival analyses.

## Usage

``` r
tcga
```

## Format

A matrix with genes as rows and cancer types as columns. Values
represent z-scores from Cox proportional hazards models. The object
shipped in the package is **sparse**: for package size,
`data-raw/shrink_reference_data.R` sets entries with \|z\| \\\leq\\ 2
(the default `z_score_cutoff` in
[`PhenoMap()`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md))
to `NA` and drops genes that are `NA` in every column. The upstream
table `TCGA_metaz.csv` (gitignored) contains full-density columns (e.g.
~20k finite values per TCGA cancer) before this step.

## Source

The Cancer Genome Atlas Research Network
