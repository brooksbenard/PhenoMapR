# List built-in marker heatmap color palette options

Returns RColorBrewer ([ColorBrewer2](https://colorbrewer2.org/)) palette
names grouped by type and viridis palette names from the
[viridis](https://github.com/sjmgarnier/viridis) family.

## Usage

``` r
list_marker_heatmap_color_palettes()
```

## Value

A list with elements `brewer` (list of `sequential`, `diverging`,
`qualitative` character vectors) and `viridis` (character vector of
palette names).
