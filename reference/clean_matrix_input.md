# Clean and normalize an expression matrix

Optionally clean gene IDs to approved HUGO symbols (collapsing duplicate
rows by mean) and apply a log-normalization that matches the chosen
analysis mode (single-cell or bulk). The function is a standalone
counterpart to the cleanup that
[`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
performs internally on the supplied bulk expression matrix; it can be
used to pre-process the expression matrix you pass to
[`PhenoMap`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md)
so that gene IDs and value scale match what the reference signatures
expect.

## Usage

``` r
clean_matrix_input(
  x,
  do_hugo = TRUE,
  do_collapse_dups = TRUE,
  do_normalize = TRUE,
  mode = c("auto", "single_cell", "bulk"),
  hugo_species = c("human", "mouse"),
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix / Matrix / data.frame with genes as rows and samples as
  columns.

- do_hugo:

  Logical. Clean gene IDs to approved HUGO symbols via HGNChelper
  (default `TRUE`). A no-op when HGNChelper is not installed; the
  function returns the matrix unchanged on the gene-ID axis and adds a
  note to `steps`.

- do_collapse_dups:

  Logical. Collapse duplicate gene rows by mean per sample (default
  `TRUE`). Always honoured even when `do_hugo = FALSE`.

- do_normalize:

  Logical. Apply log-normalization when the detected format is raw
  counts or CPM/TPM (default `TRUE`).

- mode:

  `"auto"` / `"single_cell"` / `"bulk"`. When `"auto"` the helper calls
  [`detect_expression_format`](https://brooksbenard.github.io/PhenoMapR/reference/detect_expression_format.md)
  and uses its `sc_or_bulk` guess; ties default to bulk.

- hugo_species:

  `"human"` or `"mouse"` (passed to
  [`HGNChelper::checkGeneSymbols`](https://waldronlab.io/HGNChelper/reference/checkGeneSymbols.html)).

- verbose:

  Logical. Print progress messages.

## Value

A list with elements `matrix` (the cleaned matrix), `steps` (character
vector of operations performed), `n_collapsed` (number of duplicate gene
IDs that were collapsed by mean), `detection` (the
[`detect_expression_format`](https://brooksbenard.github.io/PhenoMapR/reference/detect_expression_format.md)
result computed before any cleanup), and `mode` (the resolved mode).

## Normalization

Mode is resolved (when `mode = "auto"`) from a fresh call to
[`detect_expression_format`](https://brooksbenard.github.io/PhenoMapR/reference/detect_expression_format.md),
so callers don't have to duplicate the detection logic. The resulting
transformations are:

- single cell, raw counts: per-column library-size scaling to 10,000
  followed by `log1p`, mirroring Seurat's `LogNormalize`.

- bulk, raw counts: per-column scaling to 1e6 (CPM) followed by
  `log2(x + 1)`.

- raw CPM / TPM (any mode): `log2(x + 1)`.

- already log-normalized or z-scaled: no transformation; the function
  still cleans gene IDs if requested.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
counts <- matrix(rpois(2000, 5), nrow = 100,
                 dimnames = list(paste0("G", 1:100), paste0("S", 1:20)))
cleaned <- clean_matrix_input(counts, mode = "bulk")
str(cleaned$matrix); cleaned$steps
} # }
```
