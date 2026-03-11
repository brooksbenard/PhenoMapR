# PhenoMapR

**PhenoMapR** is a semi-supervised method to map phenotypes associated
with bulk gene expression onto single cell, spatial, and bulk
transcriptomics data. PhenoMapR nominates and rank-orders samples,
cells, and spatial locations associated with gene expression signatures
that are correlated with a phenotype of interest (e.g. overall
survival).

![PhenoMapR
visualization](reference/figures/PhenoMapR_visualization.png)

## Installation

``` r
# Download PhenoMapR using the following:
if (!require(devtools)) install.packages("devtools")
devtools::install_github("brooksbenard/PhenoMapR")
```
