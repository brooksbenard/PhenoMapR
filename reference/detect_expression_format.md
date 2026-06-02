# Detect the format of an expression matrix

Inspect a numeric expression matrix (genes x samples) and return a
classification of (1) the gene-ID style (HUGO symbols vs Ensembl IDs vs
mixed), (2) the likely expression format (raw counts / CPM or TPM /
log-normalized / z-scaled), and (3) whether the data looks like
single-cell or bulk based on sparsity and the number of samples. The
result is intended to power user-facing diagnostics (e.g., the PhenoMapR
Shiny app's 1. Data sidebar) and to drive optional cleanup with
[`clean_matrix_input`](https://brooksbenard.github.io/PhenoMapR/reference/clean_matrix_input.md).

## Usage

``` r
detect_expression_format(x, sample_cap = 10000L, verbose = FALSE)
```

## Arguments

- x:

  Matrix, Matrix sparse matrix, or data.frame with genes as rows and
  samples as columns. A data.frame is first coerced to a numeric matrix.

- sample_cap:

  Maximum number of matrix values to sample for the numeric heuristics
  (default 10000). Larger matrices are subsampled so the detector
  remains fast on very large objects.

- verbose:

  If `TRUE`, prints a one-line message summarising the detection result.

## Value

A named list with elements: `gene_id_kind`
("hugo"/"ensembl"/"mixed"/"unknown"), `n_genes`, `n_samples`,
`n_ensembl`, `n_hugo_like`, `prop_ensembl`, `prop_hugo_like`, `n_dup`
(duplicate gene IDs), `dup_examples` (up to 5), `format` (one of
`"raw_counts"`, `"cpm_or_tpm"`, `"log_normalized"`, `"z_scaled"`,
`"unknown"`), `format_label` (human-readable), `format_confidence`
("high"/"medium"/"low"), `sparsity`, `sc_or_bulk`
("single_cell"/"bulk"/"unclear"), `sc_or_bulk_label`, `stats` (numeric
vector of summary statistics), and `recommendations` (character vector
of actionable suggestions, empty when nothing is recommended).

## Details

All heuristics are intentionally conservative: when in doubt the
detector returns `"unknown"` with a low confidence rating so the caller
can prompt the user for clarification.

## Detection heuristics

- Gene-ID kind:

  Fraction of rownames matching `^ENSG\d` (human Ensembl gene IDs) vs
  the fraction matching a generic HUGO-symbol pattern (uppercase letters
  / digits / hyphens). \>50\\ \>50\\ Mouse symbols (mixed-case) are also
  counted as HUGO-like.

- Format:

  Subsamples up to `sample_cap` values: \* `"z_scaled"` if a substantial
  fraction (\>=10\\ values are negative AND the per-column mean is close
  to 0. \* `"raw_counts"` if \>=99\\ integers. \* `"cpm_or_tpm"` if
  column sums cluster near 1e6 and the matrix max is large (\>50) – TPM
  and CPM are indistinguishable without gene-length info. \*
  `"log_normalized"` if values are non-negative, mostly non-integer, and
  max is moderate (typically \<20). \* `"unknown"` otherwise.

- Single-cell vs bulk:

  Combines sparsity (fraction of zeros across the sampled values) and
  the number of samples (columns). Sparsity \>=0.5 OR (\>=200 samples
  AND sparsity \>=0.3) -\> single cell; \<200 samples AND sparsity \<0.3
  -\> bulk; everything else -\> unclear.

## Examples

``` r
set.seed(1)
counts <- matrix(rpois(2000, 5), nrow = 100,
                 dimnames = list(paste0("G", 1:100), paste0("S", 1:20)))
detect_expression_format(counts)
#> $gene_id_kind
#> [1] "hugo"
#> 
#> $n_genes
#> [1] 100
#> 
#> $n_samples
#> [1] 20
#> 
#> $n_ensembl
#> [1] 0
#> 
#> $n_hugo_like
#> [1] 100
#> 
#> $prop_ensembl
#> [1] 0
#> 
#> $prop_hugo_like
#> [1] 1
#> 
#> $n_dup
#> [1] 0
#> 
#> $dup_examples
#> character(0)
#> 
#> $format
#> [1] "raw_counts"
#> 
#> $format_label
#> [1] "Raw counts (integers)"
#> 
#> $format_confidence
#> [1] "high"
#> 
#> $sparsity
#> [1] 0.0065
#> 
#> $sc_or_bulk
#> [1] "bulk"
#> 
#> $sc_or_bulk_label
#> [1] "Bulk-like (1% zeros, 20 samples)"
#> 
#> $stats
#>           min           max          mean        median  frac_integer 
#>        0.0000       15.0000        4.9745        5.0000        1.0000 
#> frac_negative     frac_zero 
#>        0.0000        0.0065 
#> 
#> $recommendations
#> [1] "Raw counts detected. PhenoMapR's weighted-sum scoring expects log-normalized expression. Click \"Clean & normalize\" to apply log2(CPM+1) (bulk-style log normalization)."
#> 
```
