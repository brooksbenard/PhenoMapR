# Build `labels_gp` for ComplexHeatmap::anno_mark().

When `color_by_celltype` is TRUE and `cell_types` is parallel to the
labels, text is colored from `pal_celltype`; otherwise labels use
`default_col`.

## Usage

``` r
.anno_mark_labels_gp(
  cell_types,
  color_by_celltype,
  pal_celltype,
  fontsize = 7,
  default_col = "black"
)
```
