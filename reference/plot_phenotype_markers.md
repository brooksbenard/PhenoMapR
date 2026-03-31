# Plot phenotype marker gene heatmaps (ComplexHeatmap)

Draws a **cell-type-agnostic** (global) or **cell-type-specific**
heatmap from the output of
[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md).
Expression is subset to the selected marker genes and ordered cells
**before** row-wise scaling, to keep transpose/scale costs small on
large matrices.

## Usage

``` r
plot_phenotype_markers(
  markers,
  expr_mat,
  meta,
  cell_id_col = "Cell",
  group_col,
  score_col,
  celltype_col = "celltype_original",
  celltype_palette = NULL,
  heatmap_type = c("global", "cell_type_specific"),
  top_n_markers = 20L,
  n_mark_labels = 5L,
  p_adj_threshold = 0.05,
  scale_clip = NULL,
  column_title = NULL,
  draw = TRUE,
  use_raster = FALSE
)
```

## Arguments

- markers:

  List returned by
  [`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
  with elements `adverse_markers` and `favorable_markers`.

- expr_mat:

  Numeric matrix, genes \\\times\\ cells (column names = cell IDs).

- meta:

  Data frame of cell metadata with at least `cell_id_col`, `group_col`,
  `score_col`, and `celltype_col`.

- cell_id_col:

  Column in `meta` with cell IDs matching `colnames(expr_mat)`.

- group_col:

  Column with phenotype groups (`Most Favorable`, `Other`,
  `Most Adverse`).

- score_col:

  Column with continuous phenotype scores for the top color bar.

- celltype_col:

  Column with cell type labels.

- celltype_palette:

  Named vector of colors for cell types. If `NULL`,
  [`get_celltype_palette()`](https://brooksbenard.github.io/PhenoMapR/reference/get_celltype_palette.md)
  is used.

- heatmap_type:

  `"global"` (cell-type agnostic markers) or `"cell_type_specific"`
  (markers per cell type from `marker_scope = "cell_type_specific"`).

- top_n_markers:

  Maximum number of genes to keep per contrast block (per tail for
  global; per phenotype bin \\\times\\ cell type for
  cell-type-specific).

- n_mark_labels:

  Number of row labels to draw per block via
  [`ComplexHeatmap::anno_mark`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_mark.html)
  (top genes by `avg_log2FC` within each block).

- p_adj_threshold:

  Only genes with `p_adj` below this value (and positive `avg_log2FC`)
  are candidates before ordering by log fold change.

- scale_clip:

  Length-2 numeric vector `c(lo, hi)` applied after row scaling (values
  outside are clipped). If `NULL`, uses `c(-3, 3)` for
  `heatmap_type = "global"` and `c(-5, 5)` for cell-type-specific.

- column_title:

  Optional title above the heatmap.

- draw:

  If `TRUE` (default), calls
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).
  If `FALSE`, returns the `Heatmap` object invisibly.

- use_raster:

  Passed to `Heatmap()` (default `FALSE`).

## Value

Invisibly, the
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object (or `NULL` if nothing was plotted).

## Details

Requires suggested packages **ComplexHeatmap** and **circlize**.

For `heatmap_type = "global"`, rows are favorable-marker genes then
adverse-only genes (no duplicate genes in the adverse block). Columns
are ordered by PhenoMapR score (`score_col`) from **low to high** across
all cells in `expr_mat`. For `cell_type_specific`, columns follow
phenotype group, then cell type, then score within each block (see
[`find_phenotype_markers`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
grouping logic); rows follow phenotype bin then cell type, matching that
column order.

Column annotations (top to bottom): phenotype group, cell type,
PhenoMapR score (nearest the heatmap first); annotation names on the
right; legends are drawn explicitly (`Legend()` +
`annotation_legend_list`; `merge_legends = TRUE` and correct parameter
name `merge_legends`). PhenoMapR score colors are diverging
blue–white–red with **white at 0**: negative scores are blue, positive
scores are red (using the min and max of the score column across **all**
cells in `meta`). Row annotations: for `left_annotation`, the first
track is leftmost (farthest from the matrix). For
`heatmap_type = "cell_type_specific"`, the order is `anno_mark` (outer),
**cell type**, then **phenotype bin** (adjacent to the heatmap). For
`right_annotation`, strips are first (next to the heatmap): phenotype
then cell type, then `anno_mark` outermost. Favorable-tail strips and
gene marks on the **left**, adverse on the **right**. For
`heatmap_type = "cell_type_specific"`, after drawing, white outline
boxes are added around each marker-gene block (EcoTyper-style
`decorate_heatmap_body` + `grid.rect`), spanning only the columns for
that block\\s phenotype bin and cell type (not the full matrix width).
Row-split slice titles are suppressed. Heatmap fill uses ColorBrewer
**RdGy** (11-class): **high** scaled expression = red, **low** = black.
Heatmap and column annotation legends merge on the right
(`merge_legends = TRUE`; extra right `padding` for PDFs). `row_gap = 0`
between split blocks.

## See also

[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
