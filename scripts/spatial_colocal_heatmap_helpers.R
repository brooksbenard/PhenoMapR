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
  # as.integer(table) drops names; keep a named integer vector for reload.
  vals <- as.integer(tab)
  names(vals) <- names(tab)
  saveRDS(vals, file.path(data_dir, "spatial_coloc_label_ncells.rds"))
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
  # Empty title spacer + bar + CellType + Prognostic strips.
  top_anno_mm <- bar_mm + 2 * anno_mm + (if (identical(cell_mm, .spatial_hm_cell_mm_compact)) 2.5 else 4) + 2
  # Left margin for bar axis ticks + "N cells" name offset.
  left_anno_mm <- 2 * anno_mm + bar_mm * 0.45 + 18
  legend_pad <- 1.6
  # Room for multi-line column_title + gap above the annotations.
  title_pad <- 0.9
  c(
    width = as.numeric(grid::convertWidth(wh$width, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertWidth(grid::unit(left_anno_mm, "mm"), "in", valueOnly = TRUE)) +
      pad_w + legend_pad,
    height = as.numeric(grid::convertHeight(wh$height, "in", valueOnly = TRUE)) +
      as.numeric(grid::convertHeight(grid::unit(top_anno_mm, "mm"), "in", valueOnly = TRUE)) +
      pad_h + 0.8 + title_pad
  )
}

.spatial_bar_axis_at <- function(x) {
  mx <- max(x, 1L, na.rm = TRUE)
  if (mx <= 1000) {
    # Prefer few ticks so labels stay short and leave room for "N cells".
    if (mx <= 200) {
      c(0, mx)
    } else if (mx <= 500) {
      c(0, round(mx / 2), mx)
    } else {
      c(0, 500, max(1000, mx))
    }
  } else {
    c(0, round(mx / 2), mx)
  }
}

.spatial_bar_axis_labels <- function(at) {
  vapply(at, function(v) {
    if (!is.finite(v)) return("")
    if (v >= 1000) {
      paste0(format(v / 1000, trim = TRUE, digits = 2), "k")
    } else {
      as.character(as.integer(round(v)))
    }
  }, character(1))
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
  rn <- rownames(mat)
  cn <- colnames(mat)
  # Named vectors so ComplexHeatmap maps annotations by label under
  # row_order / column_order reordering (not fragile positional match).
  row_ct <- stats::setNames(
    as.character(meta_labs$CellType[match(rn, meta_labs$label)]),
    rn
  )
  row_pg <- stats::setNames(
    as.character(meta_labs$Prognostic[match(rn, meta_labs$label)]),
    rn
  )
  col_ct <- stats::setNames(
    as.character(meta_labs$CellType[match(cn, meta_labs$label)]),
    cn
  )
  col_pg <- stats::setNames(
    as.character(meta_labs$Prognostic[match(cn, meta_labs$label)]),
    cn
  )
  col_nc <- stats::setNames(as.integer(col_ncells[cn]), cn)
  col_nc[is.na(col_nc)] <- 0L
  if (anyNA(row_ct) || anyNA(col_ct) || anyNA(row_pg) || anyNA(col_pg)) {
    stop("Annotation labels failed to match matrix row/column names")
  }
  tick_candidates <- .spatial_bar_axis_at(col_nc)
  bar_ylim_max <- max(c(as.numeric(col_nc), tick_candidates, 1), na.rm = TRUE) * 1.02
  anno_mm <- if (compact) .spatial_anno_mm_compact else .spatial_anno_mm
  bar_mm <- if (compact) .spatial_bar_mm_compact else .spatial_bar_mm
  spacer_mm <- if (compact) 2.5 else 5
  anno_u <- grid::unit(anno_mm, "mm")
  bar_u <- grid::unit(bar_mm, "mm")
  # Only "N cells" needs a large offset to clear left-side tick labels;
  # CellType / Prognostic Group stay close to their strips.
  ncells_name_offset <- if (compact) 8 else 12
  name_offsets <- grid::unit(c(0, ncells_name_offset, 1.2, 1.2), "mm")
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
      title_spacer = ComplexHeatmap::anno_empty(
        border = FALSE,
        height = grid::unit(spacer_mm, "mm")
      ),
      `N cells` = ComplexHeatmap::anno_barplot(
        col_nc,
        gp = grid::gpar(fill = "#666666", col = NA),
        border = FALSE,
        height = bar_u,
        ylim = c(0, bar_ylim_max),
        axis_param = list(
          at = tick_candidates,
          labels = as.character(tick_candidates),
          gp = grid::gpar(fontsize = if (compact) 5 else 6),
          side = "left",
          labels_rot = 0
        )
      ),
      CellType = col_ct,
      `Prognostic Group` = col_pg,
      col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
      show_annotation_name = c(FALSE, TRUE, TRUE, TRUE),
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(
        fontsize = if (compact) 6 else 7,
        fontface = "plain"
      ),
      annotation_name_offset = name_offsets,
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
      annotation_height = grid::unit(
        c(spacer_mm, bar_mm, anno_mm, anno_mm),
        "mm"
      )
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
      # top / left padding: title gap + room for bar ticks and "N cells"
      padding = grid::unit(c(10, 18, 4, 28), "mm")
    ),
    right_2col = list(
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(10, 16, 4, 28), "mm")
    ),
    bottom = list(
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      padding = grid::unit(c(10, 2, 10, 2), "mm")
    )
  )
}

.spatial_draw_heatmap <- function(
    ht,
    column_title = NULL,
    column_title_gp = NULL,
    merge_legend = FALSE,
    ht_gap = grid::unit(2, "mm"),
    legend_mode = "right_2col",
    annotation_legend_list = NULL,
    heatmap_legend_list = NULL
) {
  draw_args <- .spatial_legend_draw_args(legend_mode)
  extra <- list()
  if (!is.null(annotation_legend_list)) {
    extra$annotation_legend_list <- annotation_legend_list
  }
  if (!is.null(heatmap_legend_list)) {
    extra$heatmap_legend_list <- heatmap_legend_list
  }
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
      draw_args,
      extra
    )
  )
}

.spatial_outline_legend <- function(
    dual_col = "#762A83",
    score_label = "score \u2265 0.75",
    dual_label = "Spatial + CellChat",
    title = "Outlines",
    legend_title_size = 9
) {
  ComplexHeatmap::Legend(
    title = title,
    labels = c(score_label, dual_label),
    type = "points",
    pch = 0,
    size = grid::unit(c(4, 4.5), "mm"),
    legend_gp = grid::gpar(
      col = c("black", dual_col),
      lwd = c(1.5, 2.2),
      fill = NA
    ),
    background = "white",
    title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = legend_title_size - 1)
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

.cluster_heatmap_order <- function(
    mat,
    col_fun,
    ha_row,
    ha_col,
    distance = .spatial_cor_dist
) {
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "tmp",
    col = col_fun,
    na_col = "#eeeeee",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
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
    legend_mode = "right_2col",
    annotation_legend_list = NULL,
    heatmap_legend_list = NULL
) {
  # Trailing blank lines add space between the title and the top annotations.
  title <- if (is.null(title_line2) || !nzchar(title_line2)) {
    paste0(title_line1, "\n\n")
  } else {
    paste0(.spatial_two_line_title(title_line1, title_line2), "\n\n")
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
    legend_mode = legend_mode,
    annotation_legend_list = annotation_legend_list,
    heatmap_legend_list = heatmap_legend_list
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
