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
#' adverse-only genes (no duplicate genes in the adverse block). Columns are
#' ordered by PhenoMapR score (\code{score_col}) from \strong{low to high} across
#' all cells in \code{expr_mat}. For \code{cell_type_specific}, columns follow
#' phenotype group, then cell type, then score within each block (see
#' \code{\link{find_phenotype_markers}} grouping logic); rows follow phenotype
#' bin then cell type, matching that column order.
#'
#' Column annotations (top to bottom): phenotype group, cell type, PhenoMapR
#' score (nearest the heatmap first); annotation names on the right; legends
#' are drawn explicitly (\code{Legend()} + \code{annotation_legend_list};
#' \code{merge_legends = TRUE} and correct parameter name \code{merge_legends}).
#' PhenoMapR score colors use the min--max of the score column across
#' \strong{all} cells in \code{meta}. Row annotations: for \code{left_annotation},
#' the first track is leftmost (farthest from the matrix), so \code{anno_mark}
#' is listed first and phenotype/cell-type strips last (adjacent to the heatmap);
#' for \code{right_annotation}, strips are first (next to the heatmap) and marks
#' last. Adverse-tail strips and gene marks on the \strong{left}, favorable-tail
#' strips and gene marks on the \strong{right} (global and cell-type-specific;
#' favorable \code{anno_mark} uses \code{side = "right"}).
#' Row-split slice titles are suppressed. Heatmap fill uses ColorBrewer
#' \strong{RdGy} (11-class): \strong{high} scaled expression = red, \strong{low}
#' = black. Heatmap and column annotation legends merge on the right
#' (\code{merge_legends = TRUE}; extra right \code{padding} for PDFs).
#' \code{row_gap = 0} between split blocks.
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

  if (heatmap_type == "global") {
    ord <- .global_marker_heatmap_cell_order(
      meta = meta,
      expr_mat = expr_mat,
      cell_id_col = cell_id_col,
      score_col = score_col
    )
  } else {
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
  }
  cell_order_hm <- ord$cell_order
  meta_idx_hm <- ord$meta_idx

  pal_group <- c(`Most Adverse` = "#B2182B", Other = "#F7F7F7", `Most Favorable` = "#2166AC")

  if (is.null(celltype_palette)) {
    celltype_palette <- get_celltype_palette(hm_celltype_levels)
  }
  pal_celltype <- celltype_palette[hm_celltype_levels]
  pal_celltype[is.na(pal_celltype)] <- "#BBBBBB"

  # PhenoMapR score bar: color range from ALL cells in meta (not re-scaled per block)
  score_all <- as.numeric(meta[[score_col]])
  score_all <- score_all[is.finite(score_all)]
  if (length(score_all) == 0L) {
    smin <- -1
    smax <- 1
  } else {
    smin <- min(score_all)
    smax <- max(score_all)
  }
  if (!is.finite(smin) || !is.finite(smax) || smin == smax) {
    smin <- -1
    smax <- 1
  }
  score_ann <- as.numeric(meta[[score_col]][meta_idx_hm])
  score_ann[!is.finite(score_ann)] <- NA_real_
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
    marker_tail <- c(rep("Most Favorable", n_fav), rep("Most Adverse", n_adv))

    # Top stack (CH: first = nearest heatmap): phenotype → cell type → score
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        width = grid::unit(3, "mm")
      ),
      `Cell type` = ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        width = grid::unit(3, "mm")
      ),
      `PhenoMapR score` = ComplexHeatmap::anno_simple(
        score_ann,
        col = score_col_fun
      ),
      annotation_name_side = "right",
      show_annotation_name = TRUE,
      show_legend = TRUE,
      gap = grid::unit(0, "mm")
    )

    pal_marker_row <- c(`Most Favorable` = "#2166AC", `Most Adverse` = "#B2182B")
    # Left = adverse only; right = favorable only (user request).
    strip_l <- rep(NA_character_, nrow(mat_plot))
    if (n_adv > 0L) {
      strip_l[n_fav + seq_len(n_adv)] <- "Most Adverse"
    }
    strip_r <- rep(NA_character_, nrow(mat_plot))
    strip_r[seq_len(n_fav)] <- "Most Favorable"
    marks_at_fav <- seq_len(min(n_mark_labels, n_fav))
    marks_lab_fav <- fav_genes[seq_len(min(n_mark_labels, n_fav))]
    marks_at_adv <- if (n_adv > 0L) {
      n_fav + seq_len(min(n_mark_labels, n_adv))
    } else {
      integer(0)
    }
    marks_lab_adv <- adv_only[seq_len(min(n_mark_labels, n_adv))]

    # Left (adverse): first track leftmost; put marks first, strip last (next to heatmap).
    ha_left <- NULL
    if (n_adv > 0L) {
      ha_left <- ComplexHeatmap::rowAnnotation(
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_adv,
          labels = marks_lab_adv,
          side = "left",
          labels_gp = grid::gpar(fontsize = 7),
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_l,
          col = pal_marker_row,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(18, 3), c("mm", "mm"))
      )
    }

    # Right (favorable): strip next to heatmap, gene marks on the outer right.
    ha_right <- ComplexHeatmap::rowAnnotation(
      `Phenotype` = ComplexHeatmap::anno_simple(
        strip_r,
        col = pal_marker_row,
        width = grid::unit(3, "mm"),
        na_col = "transparent"
      ),
      marks = ComplexHeatmap::anno_mark(
        at = marks_at_fav,
        labels = marks_lab_fav,
        side = "right",
        labels_gp = grid::gpar(fontsize = 7),
        link_gp = grid::gpar(col = "grey50", lwd = 0.6),
        padding = grid::unit(0.5, "mm")
      ),
      show_annotation_name = FALSE,
      gap = grid::unit(0, "mm"),
      annotation_width = grid::unit(c(3, 18), c("mm", "mm"))
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
      row_title = rep("", nlevels(row_split_g)),
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = ha_top,
      left_annotation = ha_left,
      right_annotation = ha_right,
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

    marks_at_fav <- integer(0)
    marks_lab_fav <- character(0)
    marks_at_adv <- integer(0)
    marks_lab_adv <- character(0)
    start <- 1L
    for (bk in unique(block_key)) {
      ii <- which(block_key == bk)
      n <- length(ii)
      if (n == 0L) next
      nm <- min(n_mark_labels, n)
      at <- start + seq_len(nm) - 1L
      lab <- as.character(gene_info$gene[ii[seq_len(nm)]])
      bin <- as.character(gene_info$phenotype_bin[ii[1L]])
      if (identical(bin, "Most Favorable")) {
        marks_at_fav <- c(marks_at_fav, at)
        marks_lab_fav <- c(marks_lab_fav, lab)
      } else if (identical(bin, "Most Adverse")) {
        marks_at_adv <- c(marks_at_adv, at)
        marks_lab_adv <- c(marks_lab_adv, lab)
      }
      start <- start + n
    }

    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        width = grid::unit(3, "mm")
      ),
      `Cell type` = ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        width = grid::unit(3, "mm")
      ),
      `PhenoMapR score` = ComplexHeatmap::anno_simple(
        score_ann,
        col = score_col_fun
      ),
      annotation_name_side = "right",
      show_annotation_name = TRUE,
      show_legend = TRUE,
      gap = grid::unit(0, "mm")
    )

    row_split <- factor(block_key, levels = unique(block_key))

    n_row_ct <- nrow(mat_plot)
    strip_l_pheno <- rep(NA_character_, n_row_ct)
    strip_l_ct <- rep(NA_character_, n_row_ct)
    strip_r_pheno <- rep(NA_character_, n_row_ct)
    strip_r_ct <- rep(NA_character_, n_row_ct)
    idx_fav <- which(gene_info$phenotype_bin == "Most Favorable")
    idx_adv <- which(gene_info$phenotype_bin == "Most Adverse")
    # Left = adverse strips; right = favorable strips.
    if (length(idx_adv) > 0L) {
      strip_l_pheno[idx_adv] <- "Most Adverse"
      strip_l_ct[idx_adv] <- as.character(gene_info$cell_type[idx_adv])
    }
    if (length(idx_fav) > 0L) {
      strip_r_pheno[idx_fav] <- "Most Favorable"
      strip_r_ct[idx_fav] <- as.character(gene_info$cell_type[idx_fav])
    }

    ha_left <- NULL
    if (length(idx_adv) > 0L) {
      ha_left <- ComplexHeatmap::rowAnnotation(
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_adv,
          labels = marks_lab_adv,
          side = "left",
          labels_gp = grid::gpar(fontsize = 7),
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_l_pheno,
          col = pal_group,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        `Cell type` = ComplexHeatmap::anno_simple(
          strip_l_ct,
          col = pal_celltype,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(18, 3, 3), c("mm", "mm", "mm"))
      )
    }

    ha_right <- NULL
    if (length(idx_fav) > 0L) {
      ha_right <- ComplexHeatmap::rowAnnotation(
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_r_pheno,
          col = pal_group,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        `Cell type` = ComplexHeatmap::anno_simple(
          strip_r_ct,
          col = pal_celltype,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_fav,
          labels = marks_lab_fav,
          side = "right",
          labels_gp = grid::gpar(fontsize = 7),
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(3, 3, 18), c("mm", "mm", "mm"))
      )
    }

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
      row_title = rep("", nlevels(row_split)),
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = ha_top,
      left_annotation = ha_left,
      right_annotation = ha_right,
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
    # Manual legends for top annotations (draw() expects merge_legends, not merge_legend).
    lgd_score <- ComplexHeatmap::Legend(
      title = "PhenoMapR score",
      col_fun = score_col_fun,
      title_position = "leftcenter-rot",
      direction = "vertical",
      legend_height = grid::unit(3, "cm")
    )
    lgd_group <- ComplexHeatmap::Legend(
      title = "Phenotype group",
      at = c("Most Favorable", "Other", "Most Adverse"),
      legend_gp = grid::gpar(
        fill = pal_group[c("Most Favorable", "Other", "Most Adverse")]
      )
    )
    lgd_ct <- ComplexHeatmap::Legend(
      title = "Cell type",
      at = hm_celltype_levels,
      legend_gp = grid::gpar(fill = pal_celltype[hm_celltype_levels])
    )
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      show_annotation_legend = TRUE,
      annotation_legend_list = list(lgd_score, lgd_group, lgd_ct),
      merge_legends = TRUE,
      padding = grid::unit(c(4, 4, 4, 50), "mm")
    )
  }
  invisible(ht)
}


#' ColorBrewer diverging \code{RdGy} palette, 11 classes (standard order:
#' low = red \ldots high = gray). Reversed when mapping so \strong{high} = red,
#' \strong{low} = black.
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
  circlize::colorRamp2(breaks, rev(.rdgy11_brewer))
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


#' Column order for cell-type-agnostic heatmap: increasing PhenoMapR score.
#'
#' @noRd
#' @keywords internal
.global_marker_heatmap_cell_order <- function(meta,
                                               expr_mat,
                                               cell_id_col,
                                               score_col) {
  cn <- colnames(expr_mat)
  if (is.null(cn) || length(cn) == 0L) {
    stop("'expr_mat' must have non-empty colnames (cell IDs).")
  }
  meta_ids <- meta[[cell_id_col]]
  score_vec <- as.numeric(meta[[score_col]])
  idx_match <- match(cn, meta_ids)
  sc <- rep(NA_real_, length(cn))
  ok <- which(!is.na(idx_match))
  if (length(ok) > 0L) {
    sc[ok] <- score_vec[idx_match[ok]]
  }
  o <- order(sc, seq_along(cn), na.last = TRUE)
  cell_order_hm <- cn[o]
  meta_idx_hm <- match(cell_order_hm, meta_ids)
  list(cell_order = cell_order_hm, meta_idx = meta_idx_hm)
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
