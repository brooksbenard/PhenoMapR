# Top/right/bottom/left draw() padding (mm) for anno_mark label headroom. Kept modest so draw() margins do not blow up the layout; prefer `heatmap_height` / knitr `fig.height` for tall matrices.

Top/right/bottom/left draw() padding (mm) for anno_mark label headroom.
Kept modest so draw() margins do not blow up the layout; prefer
`heatmap_height` / knitr `fig.height` for tall matrices.

## Usage

``` r
.phenomap_draw_padding_mm(
  n_rows,
  marks_at,
  fontsize = 7,
  mark_anno_width_mm = 18,
  heatmap_type = c("global", "cell_type_specific")
)
```
