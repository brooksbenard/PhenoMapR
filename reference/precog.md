# PRECOG Prognostic Z-Scores

Gene-level prognostic z-scores from the PRECOG (PREdiction of Clinical
Outcomes from Genomic Profiles) meta-analysis.

## Usage

``` r
precog
```

## Format

A matrix with genes as rows and cancer types as columns. Values
represent meta-z-scores indicating prognostic association. The object
shipped in the package is **sparse**: for package size,
`data-raw/shrink_reference_data.R` sets entries with \|z\| \\\leq\\ 2
(the default `z_score_cutoff` in
[`PhenoMap()`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md))
to `NA` and drops genes that are `NA` in every column. Raw exports have
many more finite values per column.

## Source

Gentles et al. (2015) Cell. DOI: 10.1016/j.cell.2015.11.025
