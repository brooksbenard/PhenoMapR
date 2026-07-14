# Shared ComplexHeatmap helpers for spatial co-localization + CellChat figures.

.spatial_hm_cell_mm <- 2.8
.spatial_hm_cell_mm_compact <- 1.55
.spatial_anno_mm <- 2.5
.spatial_anno_mm_compact <- 1.8
.spatial_bar_mm <- 10
.spatial_bar_mm_compact <- 6

.spatial_two_line_title <- function(line1, line2) {
  paste0(line1, "\n", line2)
}

.spatial_title_gp <- function(fontsize = 11) {
  grid::gpar(fontsize = fontsize, fontface = "bold", lineheight = 0.95)
}

.spatial_hm_wh <- function(mat, cell_mm = .spatial_hm_cell_mm) {
  list(
    width = grid::unit(ncol(mat) * cell_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_mm, "mm")
  )
}

.spatial_label_ncells_path <- function(root) {
  file.path(root, "inst", "data", "spatial_coloc_label_ncells.rds")
}

.spatial_label_ncells_from_tab <- function(labels, tab) {
  out <- rep(0L, length(labels))
  names(out) <- labels
  if (length(tab)) {
    hit <- labels %in% names(tab)
    out[hit] <- as.integer(tab[labels[hit]])
  }
  out
}

.spatial_save_label_ncells <- function(root, tab) {
  data_dir <- file.path(root, "inst", "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(as.integer(tab), file.path(data_dir, "spatial_coloc_label_ncells.rds"))
  invisible(tab)
}

.spatial_load_label_ncells <- function(root, labels) {
  rds <- .spatial_label_ncells_path(root)
  if (!file.exists(rds)) {
    return(.spatial_label_ncells_from_tab(labels, integer(0)))
  }
  .spatial_label_ncells_from_tab(labels, readRDS(rds))
}

.spatial_png_inches <- function(mat, pad_w = 4.8, pad_h = 3, cell_mm = .spatial_hm_cell_mm) {
  wh <- .spatial_hm_wh(mat, cell_mm = cell_mm)
  anno_mm <- if (identical(cell_mm, .spatial_hm_cell_mm_compact)) {
    .spatial_anno_mm_compact
  } else {
    .spatial_anno_mm
  }
  bar_mm <- if (identical(cell_mm, .spatial_hm_cell_mm_compact)) {
    .spatial_bar_mm_compact
  } else {
    .spatial_bar_mm
  }
  top_anno_mm <- bar_mm + 2 * anno_mm + 2
  # Column annotation names + barplot y-axis need extra left margin.
  left_anno_mm <- 2 * anno_mm + bar_mm * 0.45 + 3
  legend_pad <- 1.6
  c(
    width = as.numeric(grid::convertWidth(wh$width, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertWidth(grid::unit(left_anno_mm, "mm"), "in", valueOnly = TRUE)) +
      pad_w + legend_pad,
    height = as.numeric(grid::convertHeight(wh$height, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertHeight(grid::unit(top_anno_mm, "mm"), "in", valueOnly = TRUE)) +
      pad_h + 0.8
  )
}

.spatial_bar_axis_at <- function(x) {
  mx <- max(x, 1L, na.rm = TRUE)
  if (mx <= 1000) {
    pretty(c(0, mx), n = 3)
  } else {
    c(0, round(mx / 2), mx)
  }
}

.spatial_spearman_palette <- function() {
  c(
    "#283C93", "#4960C1", "#6987D6", "#91ABE4", "#C5D5F1", "#FDF9FA",
    "#E9A5B8", "#D5708F", "#BC5270", "#A83E5D", "#8D2842"
  )
}

.spatial_spearman_col_fun <- function() {
  cols <- .spatial_spearman_palette()
  breaks <- seq(-1, 1, length.out = length(cols))
  circlize::colorRamp2(breaks, cols)
}

.spatial_cor_dist <- function(x) {
  x[!is.finite(x)] <- 0
  stats::dist(x)
}

.spatial_make_ha <- function(
    mat,
    meta_labs,
    col_ncells,
    pal_ct,
    pal_pg,
    compact = FALSE,
    anno_legend_ncol = 1L,
    legend_title_size = 9
) {
  legend_title_gp <- grid::gpar(fontsize = legend_title_size, fontface = "bold")
  row_ct <- as.character(meta_labs$CellType[match(rownames(mat), meta_labs$label)])
  row_pg <- as.character(meta_labs$Prognostic[match(rownames(mat), meta_labs$label)])
  col_ct <- as.character(meta_labs$CellType[match(colnames(mat), meta_labs$label)])
  col_pg <- as.character(meta_labs$Prognostic[match(colnames(mat), meta_labs$label)])
  col_nc <- col_ncells[colnames(mat)]
  col_nc[is.na(col_nc)] <- 0L
  anno_mm <- if (compact) .spatial_anno_mm_compact else .spatial_anno_mm
  bar_mm <- if (compact) .spatial_bar_mm_compact else .spatial_bar_mm
  anno_u <- grid::unit(anno_mm, "mm")
  bar_u <- grid::unit(bar_mm, "mm")
  bar_axis_at <- .spatial_bar_axis_at(col_nc)
  list(
    row = ComplexHeatmap::rowAnnotation(
      CellType = row_ct,
      `Prognostic Group` = row_pg,
      col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
      show_annotation_name = FALSE,
      simple_anno_size = anno_u,
      show_legend = c(FALSE)
    ),
    col = ComplexHeatmap::HeatmapAnnotation(
      `N cells` = ComplexHeatmap::anno_barplot(
        col_nc,
        gp = grid::gpar(fill = "#666666", col = NA),
        border = FALSE,
        height = bar_u,
        ylim = c(0, max(col_nc, 1L, na.rm = TRUE)),
        axis_param = list(
          at = bar_axis_at,
          gp = grid::gpar(fontsize = if (compact) 5 else 6),
          side = "left",
          labels_rot = 0
        )
      ),
      CellType = col_ct,
      `Prognostic Group` = col_pg,
      col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
      show_annotation_name = TRUE,
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(
        fontsize = if (compact) 6 else 7,
        fontface = "plain"
      ),
      annotation_name_offset = grid::unit(1.2, "mm"),
      gap = grid::unit(1.2, "mm"),
      annotation_legend_param = list(
        CellType = list(
          title = "Cell type",
          ncol = anno_legend_ncol,
          title_gp = legend_title_gp,
          labels_gp = grid::gpar(fontsize = legend_title_size - 1)
        ),
        `Prognostic Group` = list(
          title = "Prognostic group",
          ncol = 1L,
          title_gp = legend_title_gp,
          labels_gp = grid::gpar(fontsize = legend_title_size - 1)
        )
      ),
      annotation_height = grid::unit(c(bar_mm, anno_mm, anno_mm), "mm")
    )
  )
}

.spatial_heatmap_legend_param <- function(title, legend_title_size = 9) {
  list(
    title = title,
    title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = legend_title_size - 1)
  )
}

.spatial_legend_draw_args <- function(mode = c("right_1col", "right_2col", "bottom")) {
  mode <- match.arg(mode)
  switch(
    mode,
    right_1col = list(
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(2, 18, 2, 24), "mm")
    ),
    right_2col = list(
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(2, 16, 2, 24), "mm")
    ),
    bottom = list(
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      padding = grid::unit(c(10, 2, 2, 2), "mm")
    )
  )
}

.spatial_draw_heatmap <- function(
    ht,
    column_title = NULL,
    column_title_gp = NULL,
    merge_legend = FALSE,
    ht_gap = grid::unit(2, "mm"),
    legend_mode = "right_2col"
) {
  draw_args <- .spatial_legend_draw_args(legend_mode)
  do.call(
    ComplexHeatmap::draw,
    c(
      list(
        ht,
        column_title = column_title,
        column_title_gp = column_title_gp,
        merge_legend = merge_legend,
        ht_gap = ht_gap
      ),
      draw_args
    )
  )
}

.build_pair_matrix <- function(labs, pair_df, value_col, fill = NA_real_) {
  mat <- matrix(fill, length(labs), length(labs), dimnames = list(labs, labs))
  if (!nrow(pair_df)) return(mat)
  for (i in seq_len(nrow(pair_df))) {
    r <- as.character(pair_df$source[i])
    c <- as.character(pair_df$target[i])
    if (r %in% labs && c %in% labs) {
      mat[r, c] <- pair_df[[value_col]][i]
    }
  }
  mat
}

.cluster_heatmap_order <- function(mat, col_fun, ha_row, ha_col) {
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "tmp",
    col = col_fun,
    na_col = "#eeeeee",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = .spatial_cor_dist,
    clustering_distance_columns = .spatial_cor_dist,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ha_row,
    top_annotation = ha_col,
    border = TRUE
  )
  grDevices::pdf(NULL)
  drawn <- ComplexHeatmap::draw(ht)
  ord <- list(
    row = ComplexHeatmap::row_order(drawn),
    col = ComplexHeatmap::column_order(drawn)
  )
  grDevices::dev.off()
  ord
}

.spatial_draw_ht_png <- function(
    ht,
    out_path,
    title_line1,
    title_line2 = NULL,
    width_in = NULL,
    height_in = NULL,
    legend_mode = "right_2col"
) {
  title <- if (is.null(title_line2) || !nzchar(title_line2)) {
    title_line1
  } else {
    .spatial_two_line_title(title_line1, title_line2)
  }
  message("Writing ", basename(out_path), " ...")
  if (is.null(width_in) || is.null(height_in)) {
    stop("width_in and height_in required for .spatial_draw_ht_png()")
  }
  grDevices::png(
    out_path,
    width = width_in,
    height = height_in,
    units = "in",
    res = 150
  )
  .spatial_draw_heatmap(
    ht,
    column_title = title,
    column_title_gp = .spatial_title_gp(10),
    legend_mode = legend_mode
  )
  grDevices::dev.off()
}

.build_compact_integrated_ht <- function(
    mat,
    name,
    col_fun,
    hm_ord,
    ha_row = NULL,
    ha_col = NULL,
    cell_fun = NULL,
    title_line1,
    title_line2 = NULL,
    show_legend = TRUE,
    legend_param = list()
) {
  wh <- .spatial_hm_wh(mat, cell_mm = .spatial_hm_cell_mm_compact)
  args <- list(
    mat,
    name = name,
    col = col_fun,
    width = wh$width,
    height = wh$height,
    na_col = "#eeeeee",
    row_order = hm_ord$row,
    column_order = hm_ord$col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    border = TRUE,
    column_title = if (is.null(title_line2) || !nzchar(title_line2)) {
      title_line1
    } else {
      .spatial_two_line_title(title_line1, title_line2)
    },
    column_title_gp = .spatial_title_gp(9),
    show_heatmap_legend = show_legend
  )
  if (!is.null(ha_row)) args$left_annotation <- ha_row
  if (!is.null(ha_col)) args$top_annotation <- ha_col
  if (!is.null(cell_fun)) args$cell_fun <- cell_fun
  if (length(legend_param)) args$heatmap_legend_param <- legend_param
  do.call(ComplexHeatmap::Heatmap, args)
}
