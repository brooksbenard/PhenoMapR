# Scoring spatial transcriptomics data with PhenoMapR

### Overview

This vignette demonstrates applying
[PhenoMapR](https://brooksbenard.github.io/PhenoMapR) to spatial
transcriptomics data in two parts. **Part 1** uses **spot-level** 10X
Visium data from an HTAN PDAC sample to show broad tissue regions of
phenotypic interest. **Part 2** uses the same sample after we mapped
paired single-cell data to the spots using
[**CytoSPACE**](https://www.nature.com/articles/s41587-023-01697-9); we
then mirror the single-cell workflow (score by cell type, identify
marker genes, and plot heatmaps) at **cell** resolution and assess
**co-localization** of prognostic cell types.

The sample (`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test`) is from a pancreatic
cancer dataset available from [HTAN](https://humantumoratlas.org/).

### Download data and vizualize H&E

``` r

suppressPackageStartupMessages({
  library(PhenoMapR)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggchicklet2)
})

knitr::opts_chunk$set(fig.width = 12, out.width = "100%", warning = FALSE)
theme_set(theme_minimal(base_size = 14))

options(googledrive_quiet = TRUE)
googledrive::drive_deauth()

gd_id_spot <- "1OkIr7ksAWxKVjtdlGqYHMidvHZZsySEE"
rds_spot <- "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds"

googledrive::drive_download(googledrive::as_id(gd_id_spot), rds_spot, overwrite = TRUE)

seurat_spot <- readRDS("HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds")

h_and_e <- SpatialPlot(object = seurat_spot,
  features = NULL,
  image.alpha = 1,
  pt.size.factor = 0) +
  guides(fill = "none") 

print(h_and_e)
```

![](spatial-transcriptomics_files/figure-html/load-setup-1.png)

### Score spots with PhenoMapR

``` r

# Score spots
scores_spot <- PhenoMapR::PhenoMap(
    expression = seurat_spot,
    reference = "precog",
    cancer_type = "Pancreatic",
    assay = "Spatial",
    slot = "data",
    verbose = FALSE
  )

# Add scores to Seurat object
seurat_spot <- PhenoMapR::add_scores_to_seurat(seurat_spot, scores_spot)

score_col_spot <- grep("weighted_sum_score", names(scores_spot), value = TRUE)[1]
## Z-scaled PhenoMapR score on spots (for hex / metadata plots)
sc_sp <- suppressWarnings(as.numeric(seurat_spot@meta.data[[score_col_spot]]))
seurat_spot$phenomapr_scaled_spot <- as.numeric(scale(sc_sp))

plot_score_distribution(
  seurat_spot$phenomapr_scaled_spot,
  main = "PRECOG Pancreatic score distribution (Visium spots)",
  base_size = 14
)
```

![](spatial-transcriptomics_files/figure-html/score-spots-part1-1.png)

### PhenoMapR score distribution across spots

``` r

p <- SpatialPlot(
  object = seurat_spot,
  features = "weighted_sum_score_Pancreatic", image.alpha = 0
)
cols_pheno <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

df <- p$data  # extract data from SpatialPlot

spot_hex_bins <- max(10, min(50L, as.integer(round(sqrt(nrow(df))))))

hex_phenomapr <- ggplot(df, aes(
  x = x,
  y = -y,
  z = scale(weighted_sum_score_Pancreatic)
)) +
  stat_summary_hex(fun = mean, bins = spot_hex_bins) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    trans = scales::pseudo_log_trans(sigma = 0.2),
    # limits = c(-4, 4),
    name = "PhenoMapR\nz-score"
  )+
  coord_fixed() +
  theme_void()


(h_and_e | hex_phenomapr) + 
  plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/plot-spot-scores-1.png) We
can see that there are clear morphological associations with where
adverse and favorable signals are being mapped to.

Next, let’s define the most adverse and favorable groups by selecting
the top and bottom 5th percentile of spots based on PhenoMapR scores.

``` r

df <- SpatialPlot(
  object = seurat_spot,
  features = "weighted_sum_score_Pancreatic", image.alpha = 0
)$data %>%
  as.data.frame()

groups_spot <- PhenoMapR::define_phenotype_groups(df, score_columns = "weighted_sum_score_Pancreatic", percentile = 0.05)

df <- left_join(df, groups_spot, by = c("cell" = "cell_id"))

# Helper function for mode
mode_val <- function(z) {
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  as.numeric(names(sort(table(z), decreasing = TRUE)[1L]))
}

hex_phenomapr_groups <- ggplot(df, aes(
  x = x,
  y = -y,
  z = as.numeric(factor(
    phenotype_group_weighted_sum_score_Pancreatic,
    levels = c("Most Favorable", "Other", "Most Adverse")
  ))
)) +
  stat_summary_hex(
    aes(
      fill = after_stat(factor(
        value,
        levels = 1:3,
        labels = c("Most Favorable", "Other", "Most Adverse")
      ))
    ),
    bins = 50,
    fun = mode_val,
    colour = NA
  ) +
  scale_fill_manual(
    values = c(
      "Most Favorable" = "#2166AC",
      "Other"            = "#F7F7F7",
      "Most Adverse"     = "#B2182B"
    ),
    drop = FALSE,
    name = "PhenoMapR\nGroup"
  ) +
  coord_fixed() +
  theme_void()


(hex_phenomapr | hex_phenomapr_groups) + 
  plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/unnamed-chunk-1-1.png)
Here, we can still see that the most favorable spots seem to cluster
together, while a subset of the most adverse also tend to co-localize.

## Part 2: CytoSPACE-mapped cells

This sample from HTAN contains a paired single-cell sample from the same
tissue block. In order to leverage the increased resolution of the
single-cell data in a spatial context, we use CytoSPACE to map the
single-cells to their spots.

### Load and score sample with PhenoMapR (CytoSPACE)

Here, we load the pre-processed Seurat object with the cells already
maped.

``` r

# Part 1 spot object is no longer needed; free memory before loading CytoSPACE cells.
for (nm in c("seurat_spot", "scores_spot", "h_and_e", "hex_phenomapr",
             "hex_phenomapr_groups", "df", "groups_spot", "p")) {
  if (exists(nm, inherits = FALSE)) rm(list = nm, inherits = FALSE)
}
gc(verbose = FALSE)
```

    ##            used  (Mb) gc trigger  (Mb)  max used  (Mb)
    ## Ncells  4390071 234.5    7859087 419.8   7859087 419.8
    ## Vcells 11270289  86.0  126615171 966.0 126265721 963.4

``` r

# Load the **CytoSPACE** object (single cells placed on Visium coordinates).
rds_cyto <- "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds"
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

googledrive::drive_download(googledrive::as_id(gd_id_cyto), rds_cyto, overwrite = TRUE)

seurat <- readRDS(rds_cyto)

  scores_spatial <- PhenoMapR::PhenoMap(
    expression = seurat,
    reference = "precog",
    cancer_type = "Pancreatic",
    assay = "Spatial",
    slot = "data",
    verbose = FALSE
  )
  seurat <- PhenoMapR::add_scores_to_seurat(seurat, scores_spatial)
score_col <- grep("weighted_sum_score", names(scores_spatial), value = TRUE)[1]

plot_score_distribution(
  seurat@meta.data[[score_col]],
  main = "PRECOG Pancreatic score distribution (CytoSPACE cells)",
  base_size = 14
)
```

![](spatial-transcriptomics_files/figure-html/score-spots-1.png)

### Score by cell type

Before looking at any spatial information, we plot the **PhenoMapR**
score by cell type to see cell types most enriched in the adverse and
favorable prognostic groups.

``` r

spatial_celltype_pal <- NULL
spatial_celltype_col <- NULL
meta_names <- names(seurat@meta.data)
celltype_col <- "CellType"

df <- seurat@meta.data
n_meta <- nrow(df)
ann_vec <- NULL
if (!is.null(celltype_col)) {
  raw <- df[[celltype_col]]
  if (is.list(raw)) {
    ann_vec <- vapply(raw, function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)[1], character(1))
  } else {
    ann_vec <- as.vector(raw)
  }
  if (length(ann_vec) != n_meta) ann_vec <- NULL
}

  df$annotation <- factor(ann_vec, exclude = NULL)

  pal <- PhenoMapR::get_celltype_palette(levels(df$annotation))
print(ggplot(df, aes(
    x = reorder(.data$annotation, .data[[score_col]], FUN = median),
    y = .data[[score_col]],
    fill = .data$annotation
  )) +
        # geom_boxplot(outlier.alpha = 0.5, median.linewidth = 0.5,
        #          outlier.fill = NULL,
        #          outlier.color = NULL, outlier.shape = 21) +
      ggchicklet2::geom_chicklet_boxplot(radius = grid::unit(1, "pt"),
                                     outlier.alpha = 0.5, median.linewidth = 0.5,
                 outlier.fill = NULL,
                 outlier.color = NULL, outlier.shape = 21) +
    scale_fill_manual(values = pal, name = celltype_col) +
    coord_cartesian(ylim = c(-10000, 15000)) +
    guides(fill = guide_legend(ncol = 2)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(y = "PRECOG Pancreatic score", x = celltype_col, title = "PhenoMapR Score by Cell Type"))
```

![](spatial-transcriptomics_files/figure-html/score-by-annotation-1.png)
As seen in the other pancreatic cancer vignettes, the ductal cell type
is the most associated with the adverse prognostic signal, while the
beta and delta cells are most associated with favorable signal.
Interestingly, plasma cells show quite a wide range of score and contain
some of the most favorably prognostic cells across the sample.

### Where are the different cell types?

Here, we use the single-cell coordinates from the spatial seurat object
and pair them with the cell level metadata in order to plot the results
in spatial context.

``` r

cell_locations <- seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::select("Cell", "row", "col", 
                "CellType", "weighted_sum_score_Pancreatic") %>%
    group_by(.data$row, .data$col) %>%
    mutate(points_per_location = n()) %>%
    ungroup()

cell_locations$CellType <- as.factor(cell_locations$CellType)

  # Shared jitter and point size so all three spatial plots match (multiple cells per spot)
  rng_row <- diff(range(cell_locations$row, na.rm = TRUE))
  rng_col <- diff(range(-cell_locations$col, na.rm = TRUE))
  spatial_jitter_w <- max(0.25, if (rng_row > 0) rng_row * 0.025 else 0.15)
  spatial_jitter_h <- max(0.25, if (rng_col > 0) rng_col * 0.025 else 0.15)
  spatial_point_range <- c(0.5, 1.6)

  spatial_celltype_pal <- PhenoMapR::get_celltype_palette(levels(cell_locations$CellType))
  
    ct_freq <- sort(table(cell_locations$CellType, useNA = "no"), decreasing = TRUE)
  ct_order <- names(ct_freq)
  cell_locations$celltype_zorder <- as.numeric(factor(as.character(cell_locations$CellType), levels = ct_order))
  ct_pal <- if (!is.null(spatial_celltype_pal)) spatial_celltype_pal else PhenoMapR::get_celltype_palette(levels(cell_locations$CellType))
 
 cytospace_loc <-  ggplot(cell_locations, aes(x = .data$row, y = -.data$col, color = .data$CellType,
    size = points_per_location, zorder = .data$celltype_zorder)) +
    geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    scale_color_manual(values = ct_pal, name = "Cell Type", na.value = "grey90") +
    scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
    guides(
      color = guide_legend(override.aes = list(size = 4), ncol = 2)
      # size  = guide_legend(override.aes = list(size = 4))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank()
    ) +
       coord_fixed(ratio = 0.6) +
  theme_void()
  
 print(cytospace_loc)
```

![](spatial-transcriptomics_files/figure-html/unnamed-chunk-2-1.png)
Interestingly, it appears that several cell types tend to co-localize
with themselves, some of which seem to be localized to the areas where
we saw the greatest adverse and favorable scores mapping to at the spot
level.

#### Where raw PhenoMapR scores are

Spatial map of PhenoMapR score (z-scaled for color gradient). Blue =
more favorable, red = more adverse.

``` r

cell_locations <- cell_locations[order(abs(cell_locations$weighted_sum_score_Pancreatic)), ]

sc_phenomapr <- ggplot(cell_locations, aes(
  x = .data$row,
  y = -.data$col,
  color = scale(weighted_sum_score_Pancreatic),
  size = .data$points_per_location
)) +
  geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
  scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    trans = scales::pseudo_log_trans(sigma = 0.2),
    name = "PhenoMapR\nz-score"
  ) +
  coord_fixed(ratio = 0.6) +
  theme_void()

(cytospace_loc | sc_phenomapr) + 
  # plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/unnamed-chunk-3-1.png)

#### Where 5th percentile cells are

In order to make the visualization more apparent, we restrict the
spatial map of prognostic groups to: top 5% (Most Adverse), bottom 5%
(Most Favorable), and the rest (Other).

``` r

cell_locations <- cell_locations %>% 
  mutate(percentile = percent_rank(weighted_sum_score_Pancreatic)) %>%
  mutate(prognostic_group = case_when(
    percentile < 0.05 ~ "Favorable",
    percentile >= 0.95 ~ "Adverse",
    TRUE ~ "Other"
  ))

  df_other <- cell_locations %>%
    dplyr::filter(prognostic_group=="Other")
  
  df_extreme <- cell_locations %>%
    dplyr::filter(prognostic_group!="Other")
  
 sc_phenomapr_5 <- ggplot() +
    geom_jitter(data = df_other, aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
      size = points_per_location), alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    geom_jitter(data = df_extreme, aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
      size = points_per_location), alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    # ggtitle("5th percentile: Most Adverse vs Most Favorable") +
    scale_color_manual(
      values = c(`Adverse` = "#B2182B", Other = "#f7f7f7", `Favorable` = "#2166AC"),
      name = "Prognostic group",
      na.value = "grey90",
      drop = FALSE
    ) +
     guides(
      color = guide_legend(override.aes = list(size = 4))
      # size  = guide_legend(override.aes = list(size = 4))
    ) +
    scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
  coord_fixed(ratio = 0.6) +
  theme_void()

  (sc_phenomapr | sc_phenomapr_5) + 
  # plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
```

![](spatial-transcriptomics_files/figure-html/unnamed-chunk-4-1.png)

Visually, adverse and favorable cells can appear to segregate on the
tissue. The next section quantifies **neighborhood co-occurrence** of
prognostic groups with
**[`spatialCooccur`](https://github.com/juninamo/spatialCooccur)**
(Inamo *et al.*, [medRxiv
2025](https://doi.org/10.1101/2025.08.05.25332835)).

### Co-localization of prognostic cells (CytoSPACE)

We use the same **Visium row / column** coordinates as the spatial
plots, assembled as a **`data.frame`** with **`x`**, **`y`**, and a
**combined label** that includes **cell type and prognostic group**.
Each cell is labeled as `CellType_prognostic_group`
(e.g. `Ductal_Adverse`). This lets us ask whether, for example,
**adverse ductal cells** co-localize differently than **favorable ductal
cells**.

Install **`spatialCooccur`** from GitHub if needed:
`remotes::install_github("juninamo/spatialCooccur")`.

**`spatialCooccur::nhood_enrichment()`** builds a **kNN graph**
(**[`Seurat::FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)**
inside the package), optionally normalizes adjacency, aggregates
**co-occurrence** between cluster pairs, and compares the observed
matrix to **permutation** nulls to produce a **z-score** matrix. We use
**`n_jobs = 1`** so the vignette runs on a single core (set higher
locally if you install the package).
**`spatialCooccur::cooccur_local()`** scores each cell by whether its
**radius** neighborhood contains both **Adverse** and **Favorable**
labels, then applies a short diffusion-style step when
**`maxnsteps > 0`**. Duplicate spot coordinates (many CytoSPACE cells
per spot) get a small jitter.

The neighborhood **Z-score** matrix is shown below as a **pre-rendered
figure** (full CytoSPACE cell set, no subsampling). Stars use
Benjamini–Hochberg FDR on two-sided p(\|Z\|): \*\*\* q\<0.001, \*\*
q\<0.01, \* q\<0.05. To regenerate after changing labels or
`spatialCooccur` settings, run
`Rscript scripts/render_spatial_colocalization_heatmap.R` from the
package root.

![](../inst/figures/spatial_colocalization_nhood_enrichment.png)

### Prognostic markers

[`find_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/find_phenotype_markers.md)
supports the same **marker scopes** on spatial data (cells mapped to
spots) as in the single-cell vignette:

- **Cell type agnostic** — `marker_scope = "phenotype_groups"`: **Most
  Adverse** vs all other labeled cells and **Most Favorable** vs all
  other labeled cells, **ignoring** cell type.
- **Cell type × phenotype vs all opposite-tail cells** —
  `marker_scope = "cell_type_specific"`,
  `celltype_contrast = "vs_opposite_tail"`: for each cell type, that
  type in one tail vs **every cell in the opposite tail** (works when a
  type exists in only one tail).
- **Cell type specific (within cell type)** —
  `marker_scope = "cell_type_specific"`, default
  `celltype_contrast = "within_cell_type"`: that type in one tail vs
  **other phenotype bins within the same type**.

Below we run all three contrasts and draw
**[`plot_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/plot_phenotype_markers.md)**
heatmaps with **five genes per block** and **five `anno_mark` labels per
block** (`top_n_markers = 5`, `n_mark_labels = 5`). Gene labels are
colored by cell type (`color_mark_labels_by_celltype = TRUE`) on the
cell-type-specific heatmaps. On CI/pkgdown builds, only the **global**
contrast is rendered to stay within memory limits; run locally for all
three heatmaps.

#### Step 1: Phenotype groups and marker discovery

``` r

markers_global <- NULL
markers_ct_within <- NULL
markers_ct_opp <- NULL

cells <- colnames(seurat)
scores_df_markers <- seurat@meta.data[cells, score_col, drop = FALSE]
rownames(scores_df_markers) <- cells
groups_markers <- PhenoMapR::define_phenotype_groups(
  scores_df_markers,
  percentile = 0.05,
  score_columns = score_col
)
group_col <- grep("phenotype_group", names(groups_markers), value = TRUE)[1]
seurat@meta.data[[group_col]] <-
  groups_markers[rownames(seurat@meta.data), group_col, drop = TRUE]

ct_col_markers <- NULL
if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col) &&
    spatial_df_celltype_col %in% names(seurat@meta.data)) {
  ct_col_markers <- spatial_df_celltype_col
} else if (exists("celltype_col") && !is.null(celltype_col) &&
           celltype_col %in% names(seurat@meta.data)) {
  ct_col_markers <- celltype_col
}

group_vec <- seurat@meta.data[cells, group_col]
group_df <- data.frame(
  cell_id = cells,
  phenotype_group = as.character(group_vec),
  stringsAsFactors = FALSE
)
if (!is.null(ct_col_markers)) {
  group_df$cell_type <- as.character(seurat@meta.data[cells, ct_col_markers])
}

max_cells_markers <- if (.is_vignette_ci()) 2500L else 5000L

markers_global <- PhenoMapR::find_phenotype_markers(
  seurat,
  group_labels = group_df,
  group_column = "phenotype_group",
  cell_id_column = "cell_id",
  marker_scope = "phenotype_groups",
  assay = "Spatial",
  slot = "data",
  max_cells_per_ident = max_cells_markers,
  verbose = FALSE
)

if (.spatial_run_full_markers && !is.null(ct_col_markers)) {
  markers_ct_opp <- PhenoMapR::find_phenotype_markers(
    seurat,
    group_labels = group_df,
    group_column = "phenotype_group",
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "cell_type_specific",
    celltype_contrast = "vs_opposite_tail",
    assay = "Spatial",
    slot = "data",
    max_cells_per_ident = max_cells_markers,
    verbose = FALSE
  )
  markers_ct_within <- PhenoMapR::find_phenotype_markers(
    seurat,
    group_labels = group_df,
    group_column = "phenotype_group",
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "cell_type_specific",
    celltype_contrast = "within_cell_type",
    assay = "Spatial",
    slot = "data",
    max_cells_per_ident = max_cells_markers,
    verbose = FALSE
  )
}
```

#### Step 2: Heatmap prep (expression matrix + metadata)

``` r

expr_pm <- NULL
meta_pm <- NULL
group_col_hm <- NULL
ct_col_hm <- NULL
pal_ct <- NULL

.get_seurat_layer <- function(obj, assay, layer) {
  mat <- tryCatch(
    Seurat::GetAssayData(obj, layer = layer, assay = assay),
    error = function(e) tryCatch(Seurat::GetAssayData(obj, slot = layer, assay = assay), error = function(e2) NULL)
  )
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    mat <- tryCatch(
      SeuratObject::LayerData(obj, layer = layer, assay = assay),
      error = function(e) NULL
    )
  }
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) NULL else mat
}

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  genes_hm <- character(0)
  for (m in list(markers_global, markers_ct_opp, markers_ct_within)) {
    if (is.null(m)) next
    for (df in list(m$adverse_markers, m$favorable_markers)) {
      if (!is.null(df) && nrow(df) > 0L && "gene" %in% names(df)) {
        genes_hm <- union(genes_hm, df$gene)
      }
    }
  }

  expr <- NULL
  layer_used <- "data"
  assay_use_hm <- if (exists("assay_use") && !is.null(assay_use)) as.character(assay_use)[1] else "Spatial"
  assay_order <- unique(c(assay_use_hm, "Spatial", "RNA", "SCT"))
  for (a in assay_order) {
    if (!a %in% names(seurat@assays)) next
    expr <- .get_seurat_layer(seurat, a, "data")
    layer_used <- "data"
    if (is.null(expr)) {
      expr <- .get_seurat_layer(seurat, a, "counts")
      layer_used <- "counts"
    }
    if (!is.null(expr)) break
  }

  if (!is.null(expr) && ncol(expr) > 0) {
    cells_expr <- colnames(expr)
    if (is.null(cells_expr)) cells_expr <- character(0)
    cells_obj <- colnames(seurat)
    cells_meta <- rownames(seurat@meta.data)
    cells_use <- intersect(cells_expr, cells_obj)
    if (length(cells_use) == 0 && length(cells_meta) > 0) {
      cells_use <- intersect(cells_expr, cells_meta)
    }
    if (length(cells_use) == 0 && length(cells_expr) > 0 && any(grepl("_", cells_expr, fixed = TRUE))) {
      stripped <- sub("^[^_]+_", "", cells_expr)
      cells_use <- intersect(stripped, cells_obj)
      if (length(cells_use) == 0) cells_use <- intersect(stripped, cells_meta)
      if (length(cells_use) > 0) {
        expr_cols <- vapply(cells_use, function(c) cells_expr[which(stripped == c)[1]], character(1))
        expr <- expr[, expr_cols, drop = FALSE]
        colnames(expr) <- cells_use
      }
    }
    if (length(cells_use) == 0) {
      message(
        "No overlapping cells between expression and metadata; skipping marker heatmaps. ",
        "expr ncol=", length(cells_expr), ", obj ncol=", length(cells_obj)
      )
    } else {
      max_cells_hm <- if (.is_vignette_ci()) 2000L else 4000L
      if (length(cells_use) > max_cells_hm) {
        set.seed(5L)
        cells_use <- sample(cells_use, max_cells_hm)
      }
      expr_sub <- expr[, cells_use, drop = FALSE]
      if (length(genes_hm) > 0) {
        genes_use <- intersect(genes_hm, rownames(expr_sub))
        if (length(genes_use) > 0) {
          expr_sub <- expr_sub[genes_use, , drop = FALSE]
        }
      }
      if (identical(layer_used, "counts")) {
        if (requireNamespace("Matrix", quietly = TRUE) &&
            (inherits(expr_sub, "Matrix") || inherits(expr_sub, "sparseMatrix"))) {
          cs <- Matrix::colSums(expr_sub)
          expr_sub <- Matrix::t(Matrix::t(expr_sub) / pmax(cs, 1) * 10000)
          expr_sub <- log1p(expr_sub)
        } else {
          cs <- colSums(expr_sub)
          expr_sub <- t(t(expr_sub) / pmax(cs, 1) * 10000)
          expr_sub <- log1p(expr_sub)
        }
      }
      expr_pm <- as.matrix(expr_sub)
      if (exists("expr", inherits = FALSE)) rm(expr)
      gc(verbose = FALSE)
      meta_pm <- seurat@meta.data[cells_use, , drop = FALSE]
      meta_pm$cell_id_plot <- cells_use
      group_col_hm <- paste0("phenotype_group_", score_col)
      if (!group_col_hm %in% names(meta_pm)) group_col_hm <- group_col
      ct_col_hm <- ct_col_markers
      if (is.null(ct_col_hm) || !ct_col_hm %in% names(meta_pm)) {
        meta_pm$celltype_plot <- factor(rep("Cell", nrow(meta_pm)))
        ct_col_hm <- "celltype_plot"
      } else {
        meta_pm[[ct_col_hm]] <- factor(as.character(meta_pm[[ct_col_hm]]))
      }
      pal_ct <- PhenoMapR::get_celltype_palette(levels(meta_pm[[ct_col_hm]]))
    }
  } else {
    message("Could not extract a non-empty expression matrix for marker heatmaps.")
  }
}
```

#### Step 3: Cell-type-agnostic marker heatmap

``` r

if (!is.null(expr_pm) && !is.null(markers_global)) {
  adv_ok <- !is.null(markers_global$adverse_markers) && nrow(markers_global$adverse_markers) > 0L
  fav_ok <- !is.null(markers_global$favorable_markers) && nrow(markers_global$favorable_markers) > 0L
  if (adv_ok && fav_ok) {
    PhenoMapR::plot_phenotype_markers(
      markers = markers_global,
      expr_mat = expr_pm,
      meta = meta_pm,
      cell_id_col = "cell_id_plot",
      group_col = group_col_hm,
      score_col = score_col,
      celltype_col = ct_col_hm,
      celltype_palette = pal_ct,
      heatmap_type = "global",
      top_n_markers = 5L,
      n_mark_labels = 5L,
      p_adj_threshold = 0.05,
      column_title = "Global phenotype marker genes (CytoSPACE cells)"
    )
  }
}
```

![](spatial-transcriptomics_files/figure-html/heatmap-markers-global-spatial-1.png)

#### Step 4: Cell-type × phenotype vs all opposite-tail cells

``` r

if (.spatial_run_full_markers && !is.null(expr_pm) && !is.null(markers_ct_opp)) {
  PhenoMapR::plot_phenotype_markers(
    markers = markers_ct_opp,
    expr_mat = expr_pm,
    meta = meta_pm,
    cell_id_col = "cell_id_plot",
    group_col = group_col_hm,
    score_col = score_col,
    celltype_col = ct_col_hm,
    celltype_palette = pal_ct,
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L,
    n_mark_labels = 5L,
    color_mark_labels_by_celltype = TRUE,
    outline_marker_blocks = TRUE,
    block_outline_color = "black",
    p_adj_threshold = 0.05,
    column_title = "Cell-type markers vs all opposite-tail cells (CytoSPACE cells)"
  )
}
```

#### Step 5: Cell-type-specific markers (within cell type)

``` r

if (.spatial_run_full_markers && !is.null(expr_pm) && !is.null(markers_ct_within)) {
  PhenoMapR::plot_phenotype_markers(
    markers = markers_ct_within,
    expr_mat = expr_pm,
    meta = meta_pm,
    cell_id_col = "cell_id_plot",
    group_col = group_col_hm,
    score_col = score_col,
    celltype_col = ct_col_hm,
    celltype_palette = pal_ct,
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L,
    n_mark_labels = 5L,
    color_mark_labels_by_celltype = TRUE,
    outline_marker_blocks = TRUE,
    block_outline_color = "black",
    p_adj_threshold = 0.05,
    column_title = "Cell-type-specific phenotype marker genes (CytoSPACE cells)"
  )
}
```

#### Stacked barplot: adverse vs. favorable by cell type

Each column is a cell type; the height is filled by the number of
adverse (5th percentile) or favorable (5th percentile) cells.

``` r

if (!is.null(group_col)) {
  ct_col_bar <- NULL
  if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col) && spatial_df_celltype_col %in% names(seurat@meta.data)) {
    ct_col_bar <- spatial_df_celltype_col
  } else if (exists("celltype_col") && !is.null(celltype_col) && celltype_col %in% names(seurat@meta.data)) {
    ct_col_bar <- celltype_col
  } else {
    for (c in c("CellType", "Celltype..major.lineage.", "cell_type", "celltype", "annotation", "seurat_clusters")) {
      if (c %in% names(seurat@meta.data) && length(unique(na.omit(seurat@meta.data[[c]]))) >= 2) {
        ct_col_bar <- c
        break
      }
    }
  }
  if (!is.null(ct_col_bar)) {
    meta <- seurat@meta.data
    meta$pg <- meta[[group_col]]
    meta$ct <- as.character(meta[[ct_col_bar]])
    meta$ct_ok <- !is.na(meta$ct) & nzchar(trimws(meta$ct))
    idx_extreme <- meta$pg %in% c("Most Adverse", "Most Favorable") & meta$ct_ok
    df_bar <- as.data.frame(table(
      CellType = meta$ct[idx_extreme],
      Prognostic_group = meta$pg[idx_extreme],
      useNA = "no"
    ))
    if (nrow(df_bar) > 0) {
      df_bar$Prognostic_group <- factor(df_bar$Prognostic_group, levels = c("Most Favorable", "Most Adverse"))
      pal_bar <- c(`Most Adverse` = "#B2182B", `Most Favorable` = "#2166AC")
      adverse <- df_bar[df_bar$Prognostic_group == "Most Adverse", c("CellType", "Freq")]
      fav <- df_bar[df_bar$Prognostic_group == "Most Favorable", c("CellType", "Freq")]
      names(adverse)[2] <- "n_adverse"
      names(fav)[2] <- "n_fav"
      df_labels <- merge(adverse, fav, by = "CellType", all = TRUE)
      df_labels$n_adverse[is.na(df_labels$n_adverse)] <- 0
      df_labels$n_fav[is.na(df_labels$n_fav)] <- 0
      df_labels$label <- paste0(df_labels$n_adverse, "/", df_labels$n_fav)
      df_labels$total <- df_labels$n_adverse + df_labels$n_fav
      ct_ord <- levels(reorder(df_bar$CellType, df_bar$Freq, function(x) -sum(x)))
      df_labels$CellType <- factor(df_labels$CellType, levels = ct_ord)
      p_bar <- ggplot(df_bar, aes(x = reorder(.data$CellType, .data$Freq, function(x) -sum(x)), y = .data$Freq, fill = .data$Prognostic_group)) +
        geom_col(position = "stack") +
        # ggchicklet2::geom_chicklet_bar(radius = grid::unit(1.5, "pt"), position = "stack") +
        geom_text(data = df_labels, aes(x = .data$CellType, y = .data$total, label = .data$n_adverse),
                  inherit.aes = FALSE, vjust = -0.3, size = 3.5, color = "#B2182B", hjust = 1.05) +
        geom_text(data = df_labels, aes(x = .data$CellType, y = .data$total, label = paste0("/", .data$n_fav)),
                  inherit.aes = FALSE, vjust = -0.3, size = 3.5, color = "#2166AC", hjust = -0.05) +
        scale_fill_manual(values = pal_bar) +
        labs(x = "Cell type", y = "Number of cells", fill = "Prognostic group",
             title = "Adverse vs. favorable cells (5th percentile) by cell type") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = "top")
      print(p_bar)
    }
  }
}
```

![](spatial-transcriptomics_files/figure-html/stacked-bar-adverse-favorable-by-celltype-1.png)
In line with the other Pancreatic cancer vignettes, PhenoMapR identifies
ductal cells (and some fibroblasts) as the most associated cell type
with an adverse PhenoMapR prognostic score. Alpha, plasma, and beta
cells are the most associated with the more favorable prognostic signal.
Interestingly, fibroblasts seems to comprise both adverse and favorable
phenotypes.

### Summary

- **Part 1 (spots)**: Spot-level object
  `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds`; **H&E** via
  `SpatialPlot`; PhenoMapR scoring; **hex** summary of z-scaled scores
  plus **point** map of 5th-percentile groups (**same `coord_fixed` /
  limits as H&E**; no jitter).
- **Part 2 (CytoSPACE)**:
  `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`; same PRECOG
  **Pancreatic** reference at **cell** resolution; score distribution,
  score by cell type, prognostic groups, and spatial maps (jittered
  cells per spot).
- **Co-localization**: pre-rendered
  **`spatialCooccur::nhood_enrichment()`** heatmap on full CytoSPACE
  cells (`CellType` × prognostic group labels); regenerate with
  `scripts/render_spatial_colocalization_heatmap.R`.
- **Markers**: three
  **[`plot_phenotype_markers()`](https://brooksbenard.github.io/PhenoMapR/reference/plot_phenotype_markers.md)**
  heatmaps on CytoSPACE cells — (1) **cell-type agnostic**
  (`heatmap_type = "global"`), (2) **cell type × phenotype vs all
  opposite-tail cells** (`celltype_contrast = "vs_opposite_tail"`), (3)
  **within cell type** (`celltype_contrast = "within_cell_type"`); each
  uses **five genes and five `anno_mark` labels per block**.

### References

**\[1\]** Benard, B. A. et al. PRECOG update: an augmented resource of
clinical outcome associations with gene expression for adult, pediatric,
and immunotherapy cohorts. Nucleic Acids Res. 54, D1579–D1589 (2026).

### Session Info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggchicklet2_0.7.0  patchwork_1.3.2    dplyr_1.2.1        ggplot2_4.0.3     
    ## [5] Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-1           PhenoMapR_0.1.0   
    ## [9] testthat_3.3.2    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     shape_1.4.6.1          jsonlite_2.0.0        
    ##   [4] magrittr_2.0.5         magick_2.9.1           spatstat.utils_3.2-3  
    ##   [7] farver_2.1.2           rmarkdown_2.31         GlobalOptions_0.1.4   
    ##  [10] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
    ##  [13] ROCR_1.0-12            spatstat.explore_3.8-1 htmltools_0.5.9       
    ##  [16] progress_1.2.3         curl_7.1.0             sass_0.4.10           
    ##  [19] sctransform_0.4.3      parallelly_1.48.0      KernSmooth_2.23-26    
    ##  [22] bslib_0.11.0           htmlwidgets_1.6.4      desc_1.4.3            
    ##  [25] ica_1.0-3              plyr_1.8.9             plotly_4.12.0         
    ##  [28] zoo_1.8-15             cachem_1.1.0           igraph_2.3.3          
    ##  [31] iterators_1.0.14       mime_0.13              lifecycle_1.0.5       
    ##  [34] pkgconfig_2.0.3        Matrix_1.7-5           R6_2.6.1              
    ##  [37] fastmap_1.2.0          clue_0.3-68            fitdistrplus_1.2-6    
    ##  [40] future_1.70.0          shiny_1.14.0           digest_0.6.39         
    ##  [43] colorspace_2.1-2       S4Vectors_0.50.1       rprojroot_2.1.1       
    ##  [46] tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
    ##  [49] pkgload_1.5.3          textshaping_1.0.5      labeling_0.4.3        
    ##  [52] progressr_0.19.0       spatstat.sparse_3.2-0  httr_1.4.8            
    ##  [55] polyclip_1.10-7        abind_1.4-8            compiler_4.6.1        
    ##  [58] gargle_1.6.1           doParallel_1.0.17      withr_3.0.3           
    ##  [61] S7_0.2.2               fastDummies_1.7.6      hexbin_1.28.5         
    ##  [64] pkgbuild_1.4.8         MASS_7.3-65            rjson_0.2.23          
    ##  [67] tools_4.6.1            lmtest_0.9-40          otel_0.2.0            
    ##  [70] googledrive_2.1.2      httpuv_1.6.17          future.apply_1.20.2   
    ##  [73] goftest_1.2-3          glue_1.8.1             nlme_3.1-169          
    ##  [76] promises_1.5.0         grid_4.6.1             Rtsne_0.17            
    ##  [79] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4        
    ##  [82] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
    ##  [85] data.table_1.18.4      hms_1.1.4              BiocGenerics_0.58.1   
    ##  [88] spatstat.geom_3.8-1    RcppAnnoy_0.0.23       foreach_1.5.2         
    ##  [91] ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1         
    ##  [94] stringr_1.6.0          limma_3.68.4           spam_2.11-4           
    ##  [97] RcppHNSW_0.7.0         later_1.4.8            circlize_0.4.18       
    ## [100] splines_4.6.1          lattice_0.22-9         survival_3.8-6        
    ## [103] deldir_2.0-4           tidyselect_1.2.1       ComplexHeatmap_2.28.0 
    ## [106] miniUI_0.1.2           pbapply_1.7-4          knitr_1.51            
    ## [109] gridExtra_2.3.1        IRanges_2.46.0         scattermore_1.2       
    ## [112] stats4_4.6.1           xfun_0.59              statmod_1.5.2         
    ## [115] brio_1.1.5             matrixStats_1.5.0      stringi_1.8.7         
    ## [118] lazyeval_0.2.3         yaml_2.3.12            evaluate_1.0.5        
    ## [121] codetools_0.2-20       tibble_3.3.1           cli_3.6.6             
    ## [124] uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0     
    ## [127] systemfonts_1.3.2      jquerylib_0.1.4        Rcpp_1.1.1-1.1        
    ## [130] globals_0.19.1         spatstat.random_3.5-0  png_0.1-9             
    ## [133] spatstat.univar_3.2-0  parallel_4.6.1         pkgdown_2.2.0         
    ## [136] presto_1.0.0           prettyunits_1.2.0      dotCall64_1.2         
    ## [139] listenv_1.0.0          viridisLite_0.4.3      scales_1.4.0          
    ## [142] ggridges_0.5.7         purrr_1.2.2            crayon_1.5.3          
    ## [145] GetoptLong_1.1.1       rlang_1.2.0            cowplot_1.2.0
