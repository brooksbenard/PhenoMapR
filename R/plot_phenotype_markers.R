#' Plot phenotype marker gene heatmaps (ComplexHeatmap)
#'
#' Draws a **cell-type-agnostic** (global) or **cell-type-specific** heatmap from
#' the output of \code{\link{find_phenotype_markers}()}. Expression is subset to
#' the selected marker genes and ordered cells **before** row-wise scaling, to keep
#' transpose/scale costs small on large matrices.
#'
#' Requires suggested packages \strong{ComplexHeatmap} and \strong{circlize}.
#'
#' @param markers List returned by \code{find_phenotype_markers()} with elements
#'   \code{adverse_markers} and \code{favorable_markers}.
#' @param expr_mat Numeric matrix, genes \eqn{\times} cells (column names = cell IDs).
#' @param meta Data frame of cell metadata with at least \code{cell_id_col},
#'   \code{group_col}, \code{score_col}, and \code{celltype_col}.
#' @param cell_id_col Column in \code{meta} with cell IDs matching \code{colnames(expr_mat)}.
#' @param group_col Column with phenotype groups (\code{Most Favorable}, \code{Other},
#'   \code{Most Adverse}).
#' @param score_col Column with continuous phenotype scores for the top color bar.
#' @param celltype_col Column with cell type labels.
#' @param celltype_palette Named vector of colors for cell types. If \code{NULL},
#'   \code{\link{get_celltype_palette}()} is used.
#' @param heatmap_type \code{"global"} (cell-type agnostic markers) or
#'   \code{"cell_type_specific"} (markers per cell type from
#'   \code{marker_scope = "cell_type_specific"}).
#' @param top_n_markers Maximum number of genes to keep per contrast block (per tail
#'   for global; per phenotype bin \eqn{\times} cell type for cell-type-specific).
#' @param n_mark_labels Number of row labels to draw per block via
#'   \code{ComplexHeatmap::anno_mark} (top genes by \code{avg_log2FC} within each block).
#' @param p_adj_threshold Only genes with \code{p_adj} below this value (and positive
#'   \code{avg_log2FC}) are candidates before ordering by log fold change.
#' @param scale_clip Length-2 numeric vector \code{c(lo, hi)} applied after row scaling
#'   (values outside are clipped). If \code{NULL}, uses \code{c(-3, 3)} for
#'   \code{heatmap_type = "global"} and \code{c(-5, 5)} for cell-type-specific.
#' @param column_title Optional title above the heatmap.
#' @param draw If \code{TRUE} (default), calls \code{ComplexHeatmap::draw()}. If
#'   \code{FALSE}, returns the \code{Heatmap} object invisibly.
#' @param use_raster Passed to \code{Heatmap()} (default \code{FALSE}).
#'
#' @return Invisibly, the \code{ComplexHeatmap::Heatmap} object (or \code{NULL} if
#'   nothing was plotted).
#'
#' @details
#' For \code{heatmap_type = "global"}, rows are favorable-marker genes then
#' adverse-only genes (no duplicate genes in the adverse block). For
#' \code{cell_type_specific}, rows follow phenotype bin then cell type, matching
#' column order.
#'
#' Column annotations (top to bottom): phenotype group, cell type, PhenoMapR
#' score; annotation names are drawn on the right. Row \code{anno_mark} gene
#' labels only (no colored row strips). Heatmap and annotation legends are
#' drawn on the right (\code{merge_legend = TRUE}). Scaled expression uses
#' ColorBrewer \strong{RdGy} (11-class): red (low) through white to gray (high).
#' No extra gap is inserted between row-split groups (\code{row_gap = 0}).
#'
#' @seealso \code{\link{find_phenotype_markers}()}
#' @export
plot_phenotype_markers <- function(markers,
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
                                   use_raster = FALSE) {
  heatmap_type <- match.arg(heatmap_type)

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    message("Install packages 'ComplexHeatmap' and 'circlize' to plot phenotype marker heatmaps.")
    return(invisible(NULL))
  }

  if (is.null(markers) || !is.list(markers)) {
    message("'markers' must be a non-NULL list from find_phenotype_markers().")
    return(invisible(NULL))
  }

  adverse_df <- markers$adverse_markers
  favorable_df <- markers$favorable_markers
  if (is.null(adverse_df) || is.null(favorable_df) ||
      nrow(adverse_df) == 0L || nrow(favorable_df) == 0L) {
    message("Marker tables are empty; skipping heatmap.")
    return(invisible(NULL))
  }

  if (is.null(scale_clip)) {
    scale_clip <- if (heatmap_type == "global") c(-3, 3) else c(-5, 5)
  }
  if (length(scale_clip) != 2L || !is.numeric(scale_clip)) {
    stop("'scale_clip' must be a numeric vector of length 2 (lower, upper clip).")
  }

  top_n_markers <- as.integer(top_n_markers)[1L]
  n_mark_labels <- as.integer(n_mark_labels)[1L]
  if (top_n_markers < 1L) stop("'top_n_markers' must be >= 1.")
  if (n_mark_labels < 1L) stop("'n_mark_labels' must be >= 1.")

  req_meta <- c(cell_id_col, group_col, score_col, celltype_col)
  if (!all(req_meta %in% names(meta))) {
    stop("meta must contain columns: ", paste(req_meta, collapse = ", "))
  }

  gene_ids <- rownames(expr_mat)
  if (is.null(gene_ids)) {
    stop("'expr_mat' must have row names (gene symbols).")
  }

  hm_group_levels <- c("Most Favorable", "Other", "Most Adverse")
  hm_celltype_levels <- levels(factor(meta[[celltype_col]]))

  ord <- .phenotype_heatmap_cell_order(
    meta = meta,
    expr_mat = expr_mat,
    cell_id_col = cell_id_col,
    group_col = group_col,
    celltype_col = celltype_col,
    score_col = score_col,
    hm_group_levels = hm_group_levels,
    hm_celltype_levels = hm_celltype_levels
  )
  cell_order_hm <- ord$cell_order
  meta_idx_hm <- ord$meta_idx

  pal_group <- c(`Most Adverse` = "#B2182B", Other = "#F7F7F7", `Most Favorable` = "#2166AC")

  if (is.null(celltype_palette)) {
    celltype_palette <- get_celltype_palette(hm_celltype_levels)
  }
  pal_celltype <- celltype_palette[hm_celltype_levels]
  pal_celltype[is.na(pal_celltype)] <- "#BBBBBB"

  score_ann <- as.numeric(meta[[score_col]][meta_idx_hm])
  score_ann[!is.finite(score_ann)] <- NA_real_
  if (all(is.na(score_ann))) score_ann <- rep(0, length(score_ann))
  smin <- suppressWarnings(min(score_ann, na.rm = TRUE))
  smax <- suppressWarnings(max(score_ann, na.rm = TRUE))
  if (!is.finite(smin) || !is.finite(smax) || smin == smax) {
    smin <- -1
    smax <- 1
  }
  score_col_fun <- circlize::colorRamp2(
    c(smin, (smin + smax) / 2, smax),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  if (heatmap_type == "global") {
    pick_global <- function(df, n_keep) {
      .pick_marker_genes(
        df, n_keep = n_keep,
        p_adj_threshold = p_adj_threshold,
        valid_genes = gene_ids
      )
    }
    fav_genes <- pick_global(favorable_df, top_n_markers)
    adv_genes <- pick_global(adverse_df, top_n_markers)
    adv_only <- adv_genes[!adv_genes %in% fav_genes]
    top_genes <- c(fav_genes, adv_only)

    if (length(top_genes) == 0L) {
      message("No global phenotype markers passed filters; skipping heatmap.")
      return(invisible(NULL))
    }

    mat_sub <- as.matrix(expr_mat[top_genes, cell_order_hm, drop = FALSE])
    mat_scaled <- t(scale(t(mat_sub)))
    mat_plot <- mat_scaled
    mat_plot[!is.finite(mat_plot)] <- NA_real_
    mat_plot[mat_plot > scale_clip[2]] <- scale_clip[2]
    mat_plot[mat_plot < scale_clip[1]] <- scale_clip[1]

    n_fav <- length(fav_genes)
    n_adv <- length(adv_only)
    marks_at_g <- c(
      seq_len(min(n_mark_labels, n_fav)),
      if (n_adv > 0L) n_fav + seq_len(min(n_mark_labels, n_adv)) else integer(0)
    )
    marks_lab_g <- c(
      fav_genes[seq_len(min(n_mark_labels, n_fav))],
      adv_only[seq_len(min(n_mark_labels, n_adv))]
    )
    marker_tail <- c(rep("Most Favorable", n_fav), rep("Most Adverse", n_adv))

    # Top stack (bottom to top): score → cell type → phenotype group
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      `PhenoMapR score` = ComplexHeatmap::anno_simple(
        score_ann,
        col = score_col_fun
      ),
      `Cell type` = ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        width = grid::unit(3, "mm")
      ),
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        width = grid::unit(3, "mm")
      ),
      annotation_name_side = "right",
      show_annotation_name = TRUE
    )

    ha_left <- ComplexHeatmap::rowAnnotation(
      marks = ComplexHeatmap::anno_mark(
        at = marks_at_g,
        labels = marks_lab_g,
        side = "left",
        labels_gp = grid::gpar(fontsize = 7),
        link_gp = grid::gpar(col = "grey50", lwd = 0.6),
        padding = grid::unit(1, "mm")
      ),
      show_annotation_name = FALSE
    )

    row_split_g <- factor(marker_tail, levels = c("Most Favorable", "Most Adverse"))
    hm_col_fun <- .scaled_expr_col_fun_rdgy11(scale_clip)

    ct <- column_title %||% "Global phenotype marker genes (favorable vs adverse)"

    ht <- ComplexHeatmap::Heatmap(
      mat_plot,
      name = "Scaled expr",
      col = hm_col_fun,
      use_raster = use_raster,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = row_split_g,
      cluster_row_slices = FALSE,
      row_gap = grid::unit(0, "mm"),
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = ha_top,
      left_annotation = ha_left,
      column_title = ct,
      heatmap_legend_param = list(
        title = "Scaled expr",
        direction = "vertical",
        title_position = "leftcenter-rot",
        legend_height = grid::unit(3, "cm")
      )
    )
  } else {
    # cell_type_specific
    gene_info <- list()
    for (g in hm_group_levels) {
      if (!g %in% c("Most Adverse", "Most Favorable")) next
      for (ct in hm_celltype_levels) {
        df_ct <- if (g == "Most Adverse") {
          adverse_df[adverse_df$cell_type == ct, , drop = FALSE]
        } else {
          favorable_df[favorable_df$cell_type == ct, , drop = FALSE]
        }
        if (nrow(df_ct) == 0L) next
        top_g <- .pick_marker_genes(
          df_ct,
          n_keep = top_n_markers,
          p_adj_threshold = p_adj_threshold,
          valid_genes = gene_ids
        )
        if (length(top_g) == 0L) next
        gene_info[[length(gene_info) + 1L]] <- data.frame(
          gene = top_g,
          cell_type = ct,
          phenotype_bin = g,
          row_id = paste(ct, g, top_g, sep = "__"),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(gene_info) == 0L) {
      message("No cell-type-specific markers passed filters; skipping heatmap.")
      return(invisible(NULL))
    }

    gene_info <- do.call(rbind, gene_info)
    block_key <- paste(
      as.character(gene_info$phenotype_bin),
      as.character(gene_info$cell_type),
      sep = "||"
    )

    mat_raw <- as.matrix(expr_mat[gene_info$gene, cell_order_hm, drop = FALSE])
    row_names_hm <- make.unique(as.character(gene_info$gene))
    rownames(mat_raw) <- row_names_hm
    mat_scaled <- t(scale(t(mat_raw)))
    mat_plot <- mat_scaled
    mat_plot[!is.finite(mat_plot)] <- NA_real_
    mat_plot[mat_plot > scale_clip[2]] <- scale_clip[2]
    mat_plot[mat_plot < scale_clip[1]] <- scale_clip[1]

    marks_at <- integer(0)
    marks_lab <- character(0)
    start <- 1L
    for (bk in unique(block_key)) {
      ii <- which(block_key == bk)
      n <- length(ii)
      if (n == 0L) next
      nm <- min(n_mark_labels, n)
      marks_at <- c(marks_at, start + seq_len(nm) - 1L)
      marks_lab <- c(marks_lab, as.character(gene_info$gene[ii[seq_len(nm)]]))
      start <- start + n
    }

    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      `PhenoMapR score` = ComplexHeatmap::anno_simple(
        score_ann,
        col = score_col_fun
      ),
      `Cell type` = ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        width = grid::unit(3, "mm")
      ),
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        width = grid::unit(3, "mm")
      ),
      annotation_name_side = "right",
      show_annotation_name = TRUE
    )

    row_split <- factor(block_key, levels = unique(block_key))

    ha_left <- ComplexHeatmap::rowAnnotation(
      marks = ComplexHeatmap::anno_mark(
        at = marks_at,
        labels = marks_lab,
        side = "left",
        labels_gp = grid::gpar(fontsize = 7),
        link_gp = grid::gpar(col = "grey50", lwd = 0.6),
        padding = grid::unit(1, "mm")
      ),
      show_annotation_name = FALSE
    )

    hm_col_fun <- .scaled_expr_col_fun_rdgy11(scale_clip)

    ct <- column_title %||% "Cell-type-specific phenotype marker genes"

    ht <- ComplexHeatmap::Heatmap(
      mat_plot,
      name = "Scaled expr",
      col = hm_col_fun,
      use_raster = use_raster,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = row_split,
      cluster_row_slices = FALSE,
      row_gap = grid::unit(0, "mm"),
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = ha_top,
      left_annotation = ha_left,
      column_title = ct,
      heatmap_legend_param = list(
        title = "Scaled expr",
        direction = "vertical",
        title_position = "leftcenter-rot",
        legend_height = grid::unit(3, "cm")
      )
    )
  }

  if (isTRUE(draw)) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legend = TRUE,
      legend_grouping = "original"
    )
  }
  invisible(ht)
}


#' ColorBrewer diverging \code{RdGy} palette, 11 classes (low = red, high = gray).
#' Hex values from \code{RColorBrewer::brewer.pal(11, "RdGy")} (no runtime dependency).
#'
#' @noRd
#' @keywords internal
.rdgy11_brewer <- c(
  "#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF",
  "#E0E0E0", "#BABABA", "#878787", "#4D4D4D", "#1A1A1A"
)

#' @noRd
#' @keywords internal
.scaled_expr_col_fun_rdgy11 <- function(scale_clip) {
  breaks <- seq(scale_clip[1], scale_clip[2], length.out = 11L)
  circlize::colorRamp2(breaks, .rdgy11_brewer)
}


#' @keywords internal
.pick_marker_genes <- function(df, n_keep, p_adj_threshold, valid_genes = NULL) {
  req <- c("gene", "avg_log2FC", "p_adj")
  if (is.null(df) || nrow(df) == 0L || !all(req %in% names(df))) {
    return(character(0))
  }
  df <- df[
    is.finite(df$avg_log2FC) & df$avg_log2FC > 0 &
      is.finite(df$p_adj) & df$p_adj < p_adj_threshold,
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0L) {
    return(character(0))
  }
  df <- df[order(-df$avg_log2FC, df$p_adj), , drop = FALSE]
  g <- head(df$gene, n_keep)
  if (!is.null(valid_genes)) {
    g <- g[g %in% valid_genes]
  }
  g
}


#' @keywords internal
.phenotype_heatmap_cell_order <- function(meta,
                                            expr_mat,
                                            cell_id_col,
                                            group_col,
                                            celltype_col,
                                            score_col,
                                            hm_group_levels,
                                            hm_celltype_levels) {
  score_vec <- meta[[score_col]]
  group_vec <- meta[[group_col]]
  ct_vec <- meta[[celltype_col]]
  cell_ids_hm <- meta[[cell_id_col]]
  cell_order_hm <- character(0)
  for (g in hm_group_levels) {
    for (ct in hm_celltype_levels) {
      idx <- which(group_vec == g & ct_vec == ct)
      if (length(idx) == 0L) next
      idx <- idx[order(score_vec[idx], na.last = TRUE)]
      cell_order_hm <- c(cell_order_hm, cell_ids_hm[idx])
    }
  }
  cell_order_hm <- c(cell_order_hm, setdiff(colnames(expr_mat), cell_order_hm))
  meta_idx_hm <- match(cell_order_hm, meta[[cell_id_col]])
  list(cell_order = cell_order_hm, meta_idx = meta_idx_hm)
}
