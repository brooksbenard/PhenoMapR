# Plot reference signature z-scores as a 1xN heatmap

Takes the gene z-score data.frame returned by
[`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)
and draws a single-row heatmap with genes ordered by z-score, top/bottom
gene labels, and a strip indicating `|z| > 2`. Requires the
ComplexHeatmap and circlize packages.

## Usage

``` r
plot_reference_signature(
  reference,
  n_label = 15L,
  z_cutoff = 2,
  row_title = NULL,
  legend_side = c("bottom", "right"),
  ...
)
```

## Arguments

- reference:

  Data.frame with genes as rownames and a single numeric column of
  z-scores (e.g. the return value of
  [`derive_reference_from_bulk`](https://brooksbenard.github.io/PhenoMapR/reference/derive_reference_from_bulk.md)).
  Can also be a numeric vector with names (gene IDs); it will be
  converted to the expected format.

- n_label:

  Integer. Number of genes to label at the low and high end of the
  z-score axis (default 15). Total labels will be at most `2 * n_label`.

- z_cutoff:

  Numeric. Absolute z-score threshold for the significance strip
  (default 2). Genes with `|z| > z_cutoff` are shown as a colored bar
  (blue for z \< -z_cutoff, red for z \> z_cutoff).

- row_title:

  Character. Title for the heatmap row / legend (e.g. "Survival
  z-score"). If `NULL`, the column name of `reference` is used when
  `reference` is a data.frame; otherwise `"z-score"`.

- legend_side:

  Character. Where to draw the heatmap legend: `"bottom"` or `"right"`
  (default `"bottom"`).

- ...:

  Optional arguments passed to
  [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  (e.g. `heatmap_legend_param`).

## Value

The return value of
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
(a `HeatmapList` object), invisibly. The plot is drawn in the current
device.

## Examples

``` r
if (FALSE) { # \dontrun{
ref <- derive_reference_from_bulk(
  bulk_expression = bulk, phenotype = pheno,
  sample_id_column = "sample_id", phenotype_column = "response",
  phenotype_type = "binary")
plot_reference_signature(ref, row_title = "Survival z-score")
} # }
```
