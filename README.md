# <a href='https://brooksbenard.github.io/PhenoMapR'><img src='man/figures/PhenoMapR_logo.png' alt='' class='logo' style='float:right; height:120px; width:auto; margin-top:0.5em;' /></a> PhenoMapR

<!-- badges: start -->

[![R-CMD-check](https://github.com/brooksbenard/PhenoMapR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brooksbenard/PhenoMapR/actions)
[![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-teal.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: Stanford](https://img.shields.io/badge/License-Stanford-yellow.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://brooksbenard.github.io/PhenoMapR)

<!-- badges: end -->

**PhenoMapR** is a semi-supervised method to map phenotypes associated with bulk gene expression onto bulk, single cell, and spatial transcriptomics data. PhenoMapR nominates and rank-orders samples, cells, and spatial locations associated with gene expression signatures that are correlated with a phenotype of interest (e.g. overall survival).

<p align="center">
  <img src="man/figures/PhenoMapR_visualization.png" alt="PhenoMapR visualization" width="420" />
</p>

## Installation

```r
# Download PhenoMapR using the following:
if (!require(devtools)) install.packages("devtools")
devtools::install_github("brooksbenard/PhenoMapR")
```
