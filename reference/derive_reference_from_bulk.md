# Derive Reference Z-Scores from Bulk Expression and Phenotype

When you have bulk expression (**genes as rows, samples as columns**)
and a phenotype (binary, continuous, or survival), this function
computes gene-level association z-scores that can be used as a custom
reference in
[`PhenoMap`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md).

## Usage

``` r
derive_reference_from_bulk(
  bulk_expression,
  phenotype,
  sample_id_column = NULL,
  phenotype_column = NULL,
  phenotype_type = c("auto", "survival", "binary", "continuous"),
  survival_time = NULL,
  survival_event = NULL,
  normalize = TRUE,
  hugo_species = c("human", "mouse"),
  binary_positive_reference = c("second", "first"),
  verbose = TRUE
)
```

## Arguments

- bulk_expression:

  Matrix or data.frame. Bulk expression with **genes as rows** and
  **samples as columns**. If the matrix appears to be samples x genes
  (e.g. fewer rows than columns), the function will transpose and
  message the user.

- phenotype:

  Data.frame with sample identifiers and phenotype column(s). Must align
  with `bulk_expression` by sample ID (see `sample_id_column`); sample
  IDs must match the **column names** of `bulk_expression`.

- sample_id_column:

  Character. Column name in `phenotype` that matches the column names of
  `bulk_expression` (sample IDs). If `NULL`, the first column of
  `phenotype` is used.

- phenotype_column:

  Character. Column name in `phenotype` for the outcome. For
  `phenotype_type = "survival"` this is ignored; use `survival_time` and
  `survival_event` instead.

- phenotype_type:

  One of `"auto"`, `"survival"`, `"binary"`, `"continuous"`. If
  `"auto"`, inferred from the phenotype column (numeric with \>2 unique
  -\> continuous; 2 unique -\> binary; or use survival if
  `survival_time` and `survival_event` are provided).

- survival_time:

  Character. Column name in `phenotype` for time-to-event (e.g. overall
  survival time). Required when `phenotype_type = "survival"`.

- survival_event:

  Character. Column name in `phenotype` for event indicator (0/1 or
  FALSE/TRUE). Required when `phenotype_type = "survival"`.

- normalize:

  Logical. If `TRUE`, run normalization when expression looks like
  counts (default `TRUE`). Set `FALSE` to skip.

- hugo_species:

  Character. Species for HUGO symbol cleaning: `"human"` or `"mouse"`
  (default `"human"`).

- binary_positive_reference:

  For `phenotype_type` `"binary"` only: which level of the binary factor
  should correspond to the **positive** outcome in logistic regression
  (`y = 1`), so that genes with **positive** z-scores are those whose
  **higher** expression is associated with that level. Use `"first"`
  when the first level of `factor(..., levels = c(...))` is your
  phenotype of interest (e.g. `mutated` vs `wt`); use `"second"` for the
  legacy convention (second level coded as `y = 1`). Default `"second"`
  preserves behaviour from previous versions when levels were implicit
  (e.g. alphabetical).

- verbose:

  Logical. Print progress messages (default `TRUE`).

## Value

A data.frame with genes as rownames and a single column of
phenotype-association z-scores, suitable for `reference` in
[`PhenoMap`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md).
When scoring with
[`PhenoMap`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md),
a **positive** weighted sum means higher expression of genes with
**positive** reference z is associated with the level you chose via
`binary_positive_reference` (for binary phenotypes). Use
`PhenoMap(..., reference_sign = -1)` if you need to flip the sign of the
entire reference after the fact.

## Details

Steps: (1) ensure genes x samples format (transpose with message if
heuristic suggests the matrix was provided as samples x genes), (2)
clean gene names to approved HUGO symbols, (3) check if expression is
normalized and normalize if needed, (4) compute phenotype association
z-scores per gene (Cox for survival, logistic regression for binary,
correlation for continuous).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated bulk: genes (rows) x samples (columns)
set.seed(1)
bulk <- matrix(rnorm(100 * 50), 100, 50,
  dimnames = list(paste0("G", 1:100), paste0("S", 1:50)))
pheno <- data.frame(
  sample_id = paste0("S", 1:50),
  response = sample(c("R", "NR"), 50, replace = TRUE))

ref <- derive_reference_from_bulk(
  bulk_expression = bulk,
  phenotype = pheno,
  sample_id_column = "sample_id",
  phenotype_column = "response",
  phenotype_type = "binary")

# Use in scoring
scores <- PhenoMap(expression = my_single_cell_data, reference = ref)
} # }
```
