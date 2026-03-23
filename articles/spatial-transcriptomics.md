# Scoring spatial transcriptomics data with PhenoMapR

## Overview

This vignette demonstrates using **{PhenoMapR}** on spatial
transcriptomics data. We use a preprocessed 10X Visium spatial object
(`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`) to score spots with
the PRECOG **Pancreatic** reference, define adverse vs. favorable
prognostic groups, and visualize both score distributions and **where**
PhenoMapR scores localize on the tissue image. The workflow mirrors the
single-cell vignette but adds spatial plots to show score and prognostic
group across the slide.

## Load data

Here, we use a 10X Visium spatial transcriptomics sample from a
pancreatic cancer dataset available from
[HTAN](https://humantumoratlas.org/). The sample has been preprocessed
and has paired single cell RNAseq data mapped to the Visium spots using
[**CytoSPACE**](https://www.nature.com/articles/s41587-023-01697-9).

``` r
suppressPackageStartupMessages(library(PhenoMapR))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(googledrive))
suppressPackageStartupMessages(library(dplyr))

knitr::opts_chunk$set(fig.width = 12, out.width = "100%", warning = FALSE)
theme_set(theme_minimal(base_size = 14))

googledrive::drive_deauth()
googledrive::drive_download(googledrive::as_id("1HM0dBrQnaNsdm5mnq23aaQ2ILofJ0_vj"), "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds", overwrite = TRUE)

  suppressMessages({
    seurat <- PhenoMapR::load_rds_fast("HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
    update_fn <- tryCatch(
      get("UpdateSeuratObject", envir = asNamespace("SeuratObject")),
      error = function(e) get0("UpdateSeuratObject", envir = asNamespace("Seurat"), mode = "function")
    )
    if (is.function(update_fn)) {
      seurat <- tryCatch(update_fn(seurat), error = function(e) seurat)
    }
  })
  assay_use <- if ("Spatial" %in% names(seurat@assays)) "Spatial" else "RNA"
  n_genes <- nrow(seurat)
  n_spots <- ncol(seurat)
  n_images <- length(seurat@images)
```

## Score sample with PhenoMapR

Apply **{PhenoMapR}** using the built-in PRECOG primary **Pancreatic**
reference [\[1\]](#ref1), then add scores to the Seurat object metadata.

``` r
suppressMessages({
  scores_spatial <- PhenoMapR::PhenoMap(
    expression = seurat,
    reference = "precog",
    cancer_type = "Pancreatic",
    assay = assay_use,
    slot = "data",
    verbose = FALSE
  )
  seurat <- PhenoMapR::add_scores_to_seurat(seurat, scores_spatial)
})
score_col <- grep("weighted_sum_score", names(scores_spatial), value = TRUE)[1]
```

## Score distribution

Distribution of PhenoMapR (PRECOG Pancreatic) scores across spots.

``` r
plot_score_distribution(
  seurat@meta.data[[score_col]],
  main = "PRECOG Pancreatic score distribution",
  base_size = 14
)
```

![](spatial-transcriptomics_files/figure-html/score-distribution-1.png)

## Score by cell type

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
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = pal, name = celltype_col) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "right"
    ) +
    labs(y = "PRECOG Pancreatic score", x = celltype_col, title = "Score by annotation (spatial)"))
```

![](spatial-transcriptomics_files/figure-html/score-by-annotation-1.png)
As see in the other Pancreatic caner vignettes, the ductal cell type is
the most associated with the adverse prognostic signal, while the beta
and delta cells are most associated with favorable signal.
Interestingly, plasma cells show quite a wide range of score and contain
some of the most favorably prognostic cells across the sample.

## Define phenotype groups

Label cells as Most Adverse (top 5%), Most Favorable (bottom 5%), or
Other, then attach to metadata for spatial plotting.

``` r
scores_df <- seurat@meta.data[, grep("weighted_sum_score", names(seurat@meta.data)), drop = FALSE]
groups <- PhenoMapR::define_phenotype_groups(scores_df, percentile = 0.05)
suppressMessages(seurat <- AddMetaData(seurat, groups))

group_col <- paste0("phenotype_group_", score_col)
if (!group_col %in% names(seurat@meta.data)) group_col <- grep("phenotype_group", names(seurat@meta.data), value = TRUE)[1]
table(seurat@meta.data[[group_col]], useNA = "ifany")
```

    ## 
    ##   Most Adverse Most Favorable          Other 
    ##           1073           1073          19312

## Build spatial locations dataframe for plotting

Here, we use the coordinates information from the spatial seurat object
and pair them with the cell level metadata in order to plot the results
in spatial contect.

``` r
coords <- tryCatch(
  suppressMessages(Seurat::GetTissueCoordinates(seurat)),
  error = function(e) NULL
)
if (is.null(coords) && length(seurat@images) > 0L) {
  img <- seurat@images[[1]]
  coords <- tryCatch(img@coordinates, error = function(e) NULL)
}
# Fallback: look for row/col or imagerow/imagecol in metadata
if (is.null(coords) || nrow(coords) == 0) {
  meta <- seurat@meta.data
  for (pair in list(
    c("row", "col"),
    c("imagerow", "imagecol"),
    c("x", "y"),
    c("tile", "cell")
  )) {
    if (pair[1] %in% names(meta) && pair[2] %in% names(meta)) {
      coords <- meta[, pair, drop = FALSE]
      names(coords) <- c("row", "col")
      break
    }
  }
}
if (!is.null(coords) && nrow(coords) > 0) {
  if (!"row" %in% names(coords)) coords$row <- coords[[grep("row|imagerow|x", names(coords), ignore.case = TRUE)[1]]]
  if (!"col" %in% names(coords)) coords$col <- coords[[grep("col|imagecol|y", names(coords), ignore.case = TRUE)[1]]]
  spatial_df <- as.data.frame(seurat@meta.data)
  spatial_df$cell_id <- rownames(spatial_df)
  # Align coordinates by cell id (coords rownames match cell barcodes)
  ids <- spatial_df$cell_id
  spatial_df$row <- coords[match(ids, rownames(coords)), "row"]
  spatial_df$col <- coords[match(ids, rownames(coords)), "col"]
  if (!is.null(group_col)) spatial_df$prognostic_group <- spatial_df[[group_col]]
  # Cell type for spatial: use spatial_celltype_col if set, else find CellType or any 2--100 level column
  celltype_plot_col <- NULL
  if (!is.null(spatial_celltype_col) && spatial_celltype_col %in% names(spatial_df)) {
    celltype_plot_col <- spatial_celltype_col
  } else if ("CellType" %in% names(spatial_df)) {
    celltype_plot_col <- "CellType"
  } else {
    idx <- grep("celltype|cell_type|annotation", names(spatial_df), ignore.case = TRUE)
    if (length(idx) > 0) celltype_plot_col <- names(spatial_df)[idx[1]]
  }
  if (is.null(celltype_plot_col)) {
    for (c in setdiff(names(spatial_df), c(score_col, group_col, "cell_id", "row", "col", "prognostic_group"))) {
      raw <- spatial_df[[c]]
      if (is.list(raw)) raw <- vapply(raw, function(x) as.character(x)[1], character(1))
      else raw <- as.character(raw)
      if (length(unique(na.omit(raw))) >= 2L && length(unique(na.omit(raw))) <= 100L) {
        celltype_plot_col <- c
        break
      }
    }
  }
  if (!is.null(celltype_plot_col)) {
    spatial_df$celltype_plot <- factor(spatial_df[[celltype_plot_col]], exclude = NULL)
    spatial_df_celltype_col <- celltype_plot_col
    n_lev <- length(levels(spatial_df$celltype_plot))
    if (n_lev >= 2L && (is.null(spatial_celltype_pal) || length(spatial_celltype_pal) != n_lev)) {
      spatial_celltype_pal <- PhenoMapR::get_celltype_palette(levels(spatial_df$celltype_plot))
      spatial_celltype_col <- celltype_plot_col
    }
  } else {
    spatial_df$celltype_plot <- NULL
    spatial_df_celltype_col <- NULL
  }
  spatial_df <- spatial_df %>%
    group_by(.data$row, .data$col) %>%
    mutate(points_per_location = n()) %>%
    ungroup()
  # Scaled PhenoMapR score (z-score) for gradient plot so colors show variation
  if (score_col %in% names(spatial_df)) {
    raw_score <- as.numeric(spatial_df[[score_col]])
    spatial_df$phenomapr_scaled <- as.numeric(scale(raw_score))
  } else {
    spatial_df$phenomapr_scaled <- NA_real_
  }
  # Shared jitter and point size so all three spatial plots match (multiple cells per spot)
  rng_row <- diff(range(spatial_df$row, na.rm = TRUE))
  rng_col <- diff(range(-spatial_df$col, na.rm = TRUE))
  spatial_jitter_w <- max(0.15, if (rng_row > 0) rng_row * 0.025 else 0.15)
  spatial_jitter_h <- max(0.15, if (rng_col > 0) rng_col * 0.025 else 0.15)
  spatial_point_range <- c(0.5, 1.6)
  has_spatial_coords <- TRUE
} else {
  has_spatial_coords <- FALSE
  message("No spatial coordinates found; skipping spatial plots.")
}
```

### Where cell types are

Here, we plot the location of the different cell types across the
sample. Points are sized based on the number of cells mapped to each
spot (CytoSPACE) and jittered in order to see all points when multiple
cells map to the same spot.

``` r
if (exists("has_spatial_coords") && has_spatial_coords && !is.null(spatial_df$celltype_plot) &&
    length(levels(spatial_df$celltype_plot)) >= 2L) {
  ct_col <- if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col)) spatial_df_celltype_col else spatial_celltype_col
  ct_freq <- sort(table(spatial_df$celltype_plot, useNA = "no"), decreasing = TRUE)
  ct_order <- names(ct_freq)
  spatial_df$celltype_zorder <- as.numeric(factor(as.character(spatial_df$celltype_plot), levels = ct_order))
  ct_pal <- if (!is.null(spatial_celltype_pal)) spatial_celltype_pal else PhenoMapR::get_celltype_palette(levels(spatial_df$celltype_plot))
  p3 <- ggplot(spatial_df, aes(x = .data$row, y = -.data$col, color = .data$celltype_plot,
    size = points_per_location, zorder = .data$celltype_zorder)) +
    geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    ggtitle(paste("Cell type (", ct_col, ")", sep = "")) +
    scale_color_manual(values = ct_pal, name = ct_col, na.value = "grey90") +
    scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  print(p3)
}
```

![](spatial-transcriptomics_files/figure-html/spatial-celltype-plot-1.png)

### Where raw PhenoMapR scores are

Spatial map of PhenoMapR score (z-scaled for color gradient). Blue =
more favorable, red = more adverse.

``` r
if (exists("has_spatial_coords") && has_spatial_coords) {
  score_vals <- spatial_df$phenomapr_scaled
  p1 <- ggplot(spatial_df, aes(x = .data$row, y = -.data$col, color = score_vals,
    size = points_per_location, zorder = score_vals)) +
    geom_jitter(alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    ggtitle("PhenoMapR score (PRECOG Pancreatic, z-scaled)") +
    scale_color_gradientn(
      colors = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
      name = "PhenoMapR\nscore (z)",
      limits = c(-2, 2),
      na.value = "grey80"
    ) +
    scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  print(p1)
}
```

![](spatial-transcriptomics_files/figure-html/spatial-score-1.png)

### Where 5th percentile cells are

In order to make the visualization more apparent, we restrict the
spatial map of prognostic groups to: top 5% (Most Adverse), bottom 5%
(Most Favorable), and the rest (Other).

``` r
if (exists("has_spatial_coords") && has_spatial_coords && !is.null(group_col)) {
  pg <- as.character(spatial_df$prognostic_group)
  df_other <- spatial_df[pg == "Other" | is.na(pg), ]
  df_extreme <- spatial_df[pg %in% c("Most Adverse", "Most Favorable"), ]
  p2 <- ggplot() +
    geom_jitter(data = df_other, aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
      size = points_per_location), alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    geom_jitter(data = df_extreme, aes(x = .data$row, y = -.data$col, color = .data$prognostic_group,
      size = points_per_location), alpha = 0.8, width = spatial_jitter_w, height = spatial_jitter_h, shape = 16) +
    ggtitle("5th percentile: Most Adverse vs Most Favorable") +
    scale_color_manual(
      values = c(`Most Adverse` = "#B2182B", Other = "#f7f7f7", `Most Favorable` = "#2166AC"),
      name = "Prognostic group",
      na.value = "grey90",
      drop = FALSE
    ) +
    scale_size_continuous(range = spatial_point_range, trans = "reverse", name = "Cells per spot") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  print(p2)
}
```

![](spatial-transcriptomics_files/figure-html/spatial-group-1.png)

It looks like the adverse and favorable cells tend to co-localize
independently (adverse with adverse, and favorable with favorable).

## Prognostic markers

We first find marker genes for adverse vs. favorable spots, visualize
them in a heatmap (similar to the single-cell vignette), then define
**unique markers for each adverse/favorable cell type combination**
(e.g. adverse ductal, favorable B cell) and visualize those in a second
heatmap.

### Step 1: Markers for adverse vs. favorable spots

``` r
markers <- NULL
assay_markers <- NULL
if (!is.null(group_col)) {
  cells <- colnames(seurat)
  group_vec <- seurat@meta.data[cells, group_col]
  group_df <- data.frame(
    cell_id = cells,
    phenotype_group = as.character(group_vec),
    stringsAsFactors = FALSE
  )
  for (a in unique(c(assay_use, "RNA", "SCT"))) {
    if (!a %in% names(seurat@assays)) next
    markers <- tryCatch(
      PhenoMapR::find_phenotype_markers(
        seurat,
        group_labels = group_df,
        group_column = "phenotype_group",
        cell_id_column = "cell_id",
        assay = a,
        slot = "data",
        max_cells_per_ident = 5000L,
        verbose = FALSE
      ),
      error = function(e) {
        message("find_phenotype_markers (assay ", a, ") error: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(markers) && (nrow(markers$adverse_markers) > 0 || nrow(markers$favorable_markers) > 0)) {
      assay_markers <- a
      break
    }
  }
  if (!is.null(markers)) {
    message("Adverse markers (top 5):")
    print(head(markers$adverse_markers, 5))
    message("Favorable markers (top 5):")
    print(head(markers$favorable_markers, 5))
  } else {
    message("Marker analysis returned no results. Check that Most Adverse and Most Favorable groups have enough cells.")
  }
}
```

    ##   p_val avg_log2FC pct_in_group pct_rest     gene p_adj
    ## 1     0   3.221594        0.829    0.100      MET     0
    ## 2     0   2.425253        0.870    0.159    MECOM     0
    ## 3     0   3.033991        0.916    0.210    ITGA2     0
    ## 4     0   2.679157        0.851    0.152 BAIAP2L1     0
    ## 5     0   2.662331        0.801    0.105    RASEF     0

    ##   p_val avg_log2FC pct_in_group pct_rest     gene p_adj
    ## 1     0   4.271152        0.627    0.064    NRXN1     0
    ## 2     0   5.161995        0.626    0.066    RIMS2     0
    ## 3     0   5.324439        0.640    0.081   KCNMB2     0
    ## 4     0   5.275982        0.595    0.041    ABCC8     0
    ## 5     0   5.057609        0.589    0.036 TMEM132D     0

### Step 2: Heatmap of adverse vs. favorable markers

``` r
if (exists("markers") && !is.null(markers)) {
  n_top <- 15
  expr <- NULL
  assay_order <- unique(c(
    if (exists("assay_markers") && !is.null(assay_markers)) assay_markers else character(0),
    assay_use, "RNA", "SCT"
  ))
  for (a in assay_order) {
    if (!a %in% names(seurat@assays)) next
    expr <- tryCatch(
      Seurat::GetAssayData(seurat, layer = "data", assay = a),
      error = function(e) tryCatch(Seurat::GetAssayData(seurat, slot = "data", assay = a), error = function(e2) NULL)
    )
    if (!is.null(expr) && nrow(expr) > 0 && ncol(expr) > 0) break
    expr <- tryCatch(
      SeuratObject::LayerData(seurat, layer = "data", assay = a),
      error = function(e) NULL
    )
    if (is.null(expr) || ncol(expr) == 0) {
      expr <- tryCatch(
        SeuratObject::LayerData(seurat, layer = "counts", assay = a),
        error = function(e) tryCatch(SeuratObject::LayerData(seurat, layer = NULL, assay = a), error = function(e2) NULL)
      )
    }
    if (!is.null(expr) && nrow(expr) > 0 && ncol(expr) > 0) break
  }
  if (!is.null(expr) && ncol(expr) > 0) {
  adverse_pos <- markers$adverse_markers
  if (nrow(adverse_pos) > 0 && "avg_log2FC" %in% names(adverse_pos)) adverse_pos <- adverse_pos[adverse_pos$avg_log2FC > 0, ]
  favorable_pos <- markers$favorable_markers
  if (nrow(favorable_pos) > 0 && "avg_log2FC" %in% names(favorable_pos)) favorable_pos <- favorable_pos[favorable_pos$avg_log2FC > 0, ]
  pcol <- if ("p_adj" %in% names(markers$adverse_markers)) "p_adj" else if ("p_val_adj" %in% names(markers$adverse_markers)) "p_val_adj" else "p_val"
  top_genes <- unique(c(
    if (nrow(adverse_pos) > 0 && "gene" %in% names(adverse_pos)) head(adverse_pos$gene[order(adverse_pos[[pcol]])], n_top) else character(0),
    if (nrow(favorable_pos) > 0 && "gene" %in% names(favorable_pos)) head(favorable_pos$gene[order(favorable_pos[[pcol]])], n_top) else character(0)
  ))
  top_genes <- top_genes[top_genes %in% rownames(expr)]
  if (length(top_genes) == 0) top_genes <- head(rownames(expr), 20)

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
    message("No overlapping cells between expression and metadata; skipping heatmap. ",
            "expr ncol=", length(cells_expr), ", obj ncol=", length(cells_obj),
            if (length(cells_expr) > 0) paste0("; expr sample: ", head(cells_expr, 2)) else "")
  } else {
  mat <- as.matrix(expr[top_genes, cells_use, drop = FALSE])
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[mat_scaled < -3] <- -3
  mat_scaled[mat_scaled > 3] <- 3
  meta <- seurat@meta.data[cells_use, , drop = FALSE]
  ord <- order(meta[[score_col]])
  mat_scaled <- mat_scaled[, ord, drop = FALSE]
  meta_ord <- meta[ord, , drop = FALSE]

  group_col_heatmap <- paste0("phenotype_group_", score_col)
  if (!group_col_heatmap %in% names(meta_ord)) group_col_heatmap <- group_col

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann_col <- data.frame(
      `PhenoMapR Score` = meta_ord[[score_col]],
      `Prognostic group` = factor(meta_ord[[group_col_heatmap]], levels = c("Most Adverse", "Other", "Most Favorable")),
      check.names = FALSE
    )
    ct_col_ann <- if (exists("spatial_df_celltype_col") && !is.null(spatial_df_celltype_col)) spatial_df_celltype_col else if (exists("celltype_col") && !is.null(celltype_col)) celltype_col else NULL
    if (!is.null(ct_col_ann) && ct_col_ann %in% names(meta_ord)) {
      ann_col$`Cell type` <- factor(meta_ord[[ct_col_ann]])
    }
    rownames(ann_col) <- colnames(mat_scaled)
    pal_score <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    pal_group <- c(`Most Adverse` = "#B2182B", Other = "#f7f7f7", `Most Favorable` = "#2166AC")
    pal_celltype <- if ("Cell type" %in% names(ann_col)) {
      PhenoMapR::get_celltype_palette(levels(ann_col$`Cell type`))
    } else list()
    ann_colors <- list(
      `PhenoMapR Score` = pal_score,
      `Prognostic group` = pal_group
    )
    if (length(pal_celltype) > 0) ann_colors$`Cell type` <- pal_celltype
    heatmap_colors <- if (requireNamespace("paletteer", quietly = TRUE)) {
      colorRampPalette(paletteer::paletteer_d("MexBrewer::Vendedora"))(100)
    } else {
      colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    }
    pheatmap::pheatmap(
      mat_scaled,
      scale = "none",
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      color = heatmap_colors,
      breaks = seq(-3, 3, length.out = 101),
      main = "Top adverse & favorable marker genes (spatial spots)",
      fontsize = 12,
      fontsize_row = 10,
      treeheight_row = 25
    )
  }
  }
  } else {
    msg <- "Could not retrieve expression data for heatmap"
    if (!is.null(expr)) {
      msg <- paste0(msg, " (expr nrow=", nrow(expr), ", ncol=", ncol(expr), "; try LayerData for Seurat v5)")
    }
    message(msg, ".")
  }
}
```

![](spatial-transcriptomics_files/figure-html/heatmap-markers-spatial-1.png)

### Stacked barplot: adverse vs. favorable by cell type

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

## Summary

- **Data**: Processed spatial object
  (`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`).
- **Scoring**: PhenoMapR with PRECOG **Pancreatic**; scores added to
  Seurat metadata.
- **Plots**: Score distribution, score by annotation (if present),
  prognostic groups, and **spatial maps** of cell type, continuous
  score, and prognostic group (points ordered so less frequent groups
  are drawn on top for visibility).
- **Markers**: (1) Adverse vs. favorable markers with heatmap; (2)
  unique markers for each adverse/favorable cell type combination
  (e.g. adverse ductal, favorable B cell) with a second heatmap.

## References

**\[1\]** Benard, B. A. et al. PRECOG update: an augmented resource of
clinical outcome associations with gene expression for adult, pediatric,
and immunotherapy cohorts. Nucleic Acids Res. 54, D1579–D1589 (2026).

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
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
    ## [1] dplyr_1.2.0        googledrive_2.1.2  ggplot2_4.0.2      Seurat_5.4.0      
    ## [5] SeuratObject_5.3.0 sp_2.2-1           PhenoMapR_0.1.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.4        
    ##   [4] spatstat.utils_3.2-2   farver_2.1.2           rmarkdown_2.30        
    ##   [7] fs_2.0.0               ragg_1.5.1             vctrs_0.7.2           
    ##  [10] ROCR_1.0-12            spatstat.explore_3.8-0 paletteer_1.7.0       
    ##  [13] htmltools_0.5.9        curl_7.0.0             sass_0.4.10           
    ##  [16] sctransform_0.4.3      parallelly_1.46.1      KernSmooth_2.23-26    
    ##  [19] bslib_0.10.0           htmlwidgets_1.6.4      desc_1.4.3            
    ##  [22] ica_1.0-3              plyr_1.8.9             plotly_4.12.0         
    ##  [25] zoo_1.8-15             cachem_1.1.0           igraph_2.2.2          
    ##  [28] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
    ##  [31] Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0         
    ##  [34] fitdistrplus_1.2-6     future_1.70.0          shiny_1.13.0          
    ##  [37] digest_0.6.39          rematch2_2.1.2         patchwork_1.3.2       
    ##  [40] tensor_1.5.1           prismatic_1.1.2        RSpectra_0.16-2       
    ##  [43] irlba_2.3.7            textshaping_1.0.5      labeling_0.4.3        
    ##  [46] progressr_0.18.0       spatstat.sparse_3.1-0  httr_1.4.8            
    ##  [49] polyclip_1.10-7        abind_1.4-8            compiler_4.5.3        
    ##  [52] gargle_1.6.1           withr_3.0.2            S7_0.2.1              
    ##  [55] fastSave_0.1.0         fastDummies_1.7.5      MASS_7.3-65           
    ##  [58] tools_4.5.3            lmtest_0.9-40          otel_0.2.0            
    ##  [61] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
    ##  [64] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
    ##  [67] grid_4.5.3             Rtsne_0.17             cluster_2.1.8.2       
    ##  [70] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
    ##  [73] spatstat.data_3.1-9    tidyr_1.3.2            data.table_1.18.2.1   
    ##  [76] spatstat.geom_3.7-2    RcppAnnoy_0.0.23       ggrepel_0.9.8         
    ##  [79] RANN_2.6.2             pillar_1.11.1          stringr_1.6.0         
    ##  [82] spam_2.11-3            RcppHNSW_0.6.0         limma_3.66.0          
    ##  [85] later_1.4.8            splines_4.5.3          lattice_0.22-9        
    ##  [88] survival_3.8-6         deldir_2.0-4           tidyselect_1.2.1      
    ##  [91] miniUI_0.1.2           pbapply_1.7-4          knitr_1.51            
    ##  [94] gridExtra_2.3          scattermore_1.2        xfun_0.57             
    ##  [97] statmod_1.5.1          matrixStats_1.5.0      pheatmap_1.0.13       
    ## [100] stringi_1.8.7          lazyeval_0.2.2         yaml_2.3.12           
    ## [103] evaluate_1.0.5         codetools_0.2-20       tibble_3.3.1          
    ## [106] cli_3.6.5              uwot_0.2.4             xtable_1.8-8          
    ## [109] reticulate_1.45.0      systemfonts_1.3.2      jquerylib_0.1.4       
    ## [112] Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-5 
    ## [115] png_0.1-9              spatstat.univar_3.1-7  parallel_4.5.3        
    ## [118] pkgdown_2.2.0          presto_1.0.0           dotCall64_1.2         
    ## [121] listenv_0.10.1         viridisLite_0.4.3      scales_1.4.0          
    ## [124] ggridges_0.5.7         purrr_1.2.1            rlang_1.1.7           
    ## [127] cowplot_1.2.0
