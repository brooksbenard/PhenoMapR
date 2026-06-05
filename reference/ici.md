# Immune Checkpoint Inhibitor (ICI) Prognostic Z-Scores

Gene-level prognostic z-scores for patients treated with immune
checkpoint inhibitors, separated by primary and metastatic disease.

## Usage

``` r
ici
```

## Format

A matrix with genes as rows and columns in format
`CANCER_THERAPY_STAGE_ENDPOINT` (e.g. `"SKCM_PD1_Primary_Naive"`).
Unlike the other built-in references, the ICI table is shipped at a
**relaxed** cutoff: `data-raw/shrink_reference_data.R` keeps every entry
with `|z| > 1` and only sets `|z| \eqn{\leq} 1` to `NA`. PRECOG / TCGA /
pediatric ship at `|z| > 2`, but several ICI columns drop to **zero**
finite values once trimmed at 2; relaxing to 1 keeps all 47 ICI columns
usable while only adding ~1 MB to the install.
