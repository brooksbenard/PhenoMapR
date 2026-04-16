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
#' @param rank_by How to rank genes when selecting the top markers for the heatmap.
#'   \code{"lfc"} ranks by decreasing \code{avg_log2FC} among significant genes
#'   (\code{p_adj < p_adj_threshold}). \code{"p_adj"} ranks by increasing
#'   \code{p_adj} but additionally requires \code{avg_log2FC > 1} (and
#'   \code{p_adj < p_adj_threshold}) so that very small effects are not selected
#'   purely due to sample size.
#' @param n_mark_labels Number of row labels to draw per block via
#'   \code{ComplexHeatmap::anno_mark} (top genes by \code{avg_log2FC} within each block).
#' @param p_adj_threshold Only genes with \code{p_adj} below this value (and positive
#'   \code{avg_log2FC}) are candidates before ordering by log fold change.
#' @param scale_clip Length-2 numeric vector \code{c(lo, hi)} applied after row scaling
#'   (values outside are clipped). If \code{NULL}, uses \code{c(-3, 3)} for
#'   \code{heatmap_type = "global"} and \code{c(-5, 5)} for cell-type-specific.
#' @param heatmap_width Optional heatmap width passed to
#'   \code{ComplexHeatmap::Heatmap(width = ...)}. Accepts a \code{grid::unit} or a
#'   single number interpreted as millimeters.
#' @param heatmap_height Optional heatmap height passed to
#'   \code{ComplexHeatmap::Heatmap(height = ...)}. Accepts a \code{grid::unit} or a
#'   single number interpreted as millimeters.
#' @param column_title Optional title above the heatmap.
#' @param draw If \code{TRUE} (default), calls \code{ComplexHeatmap::draw()}. If
#'   \code{FALSE}, returns the \code{Heatmap} object invisibly.
#' @param use_raster Passed to \code{Heatmap()} (default \code{FALSE}).
#' @param outline_marker_blocks If \code{TRUE} (default) and
#'   \code{heatmap_type = "cell_type_specific"}, draws white outline boxes around
#'   each marker-gene block after \code{draw()}. Set \code{FALSE} to omit them.
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
#' PhenoMapR score colors are diverging blue--white--red with \strong{white at
#' 0}: negative scores are blue, positive scores are red. The score column is
#' coerced to a plain numeric vector (so matrix columns, e.g. from
#' \code{scale()}, do not corrupt subsetting). Limits use \strong{heatmap column}
#' scores only; the merged legend uses an explicit \code{Legend(at = ...)} with
#' tick positions rounded to whole numbers (annotation auto-legends are disabled). Row
#' annotations: for \code{left_annotation}, the first track is leftmost
#' (farthest from the matrix). For \code{heatmap_type = "cell_type_specific"},
#' the order is \code{anno_mark} (outer), \strong{cell type}, then
#' \strong{phenotype bin} (adjacent to the heatmap). For \code{right_annotation},
#' strips are first (next to the heatmap): phenotype then cell type, then
#' \code{anno_mark} outermost. Favorable-tail strips and gene marks on the
#' \strong{left}, adverse on the \strong{right}. For
#' \code{heatmap_type = "cell_type_specific"}, optional white outline boxes
#' (\code{outline_marker_blocks}) use \code{decorate_heatmap_body} per row slice
#' with column span in native units so each box covers only the columns whose
#' \code{group_col}/\code{celltype_col} match that block (see
#' \code{outline_marker_blocks}).
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
                                   rank_by = c("lfc", "p_adj"),
                                   n_mark_labels = 5L,
                                   p_adj_threshold = 0.05,
                                   scale_clip = NULL,
                                   heatmap_width = NULL,
                                   heatmap_height = NULL,
                                   column_title = NULL,
                                   draw = TRUE,
                                   use_raster = FALSE,
                                   outline_marker_blocks = TRUE) {
  heatmap_type <- match.arg(heatmap_type)
  rank_by <- match.arg(rank_by)

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

  .as_heatmap_unit_mm <- function(x, arg_name) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "unit")) return(x)
    if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
      return(grid::unit(as.numeric(x), "mm"))
    }
    stop("'", arg_name, "' must be NULL, a grid::unit, or a single finite number (mm).")
  }
  heatmap_width_u <- .as_heatmap_unit_mm(heatmap_width, "heatmap_width")
  heatmap_height_u <- .as_heatmap_unit_mm(heatmap_height, "heatmap_height")

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

  # PhenoMapR score bar: diverging blue–white–red with white at 0. Color limits
  # use scores for heatmap columns only (aligned to cell_order_hm), so the
  # legend matches what is plotted. Coerce meta[[score_col]] to a plain numeric
  # vector first (1-column matrix / scale() columns would otherwise break
  # subsetting and inflate the legend range). Fall back to all meta scores if
  # no finite scores on heatmap columns.
  score_ann <- .score_ann_for_heatmap(meta, score_col, meta_idx_hm)
  sa <- score_ann[is.finite(score_ann)]
  if (length(sa) > 0L) {
    smin <- min(sa)
    smax <- max(sa)
  } else {
    score_all <- .meta_score_numeric_vector(meta, score_col)
    score_all <- score_all[is.finite(score_all)]
    if (length(score_all) == 0L) {
      smin <- -1
      smax <- 1
    } else {
      smin <- min(score_all)
      smax <- max(score_all)
    }
  }
  if (!is.finite(smin) || !is.finite(smax) || smin == smax) {
    smin <- -1
    smax <- 1
  }
  score_col_fun <- .phenomap_score_col_fun(smin, smax)

  decorate_ct_rect <- NULL

  if (heatmap_type == "global") {
    pick_global <- function(df, n_keep) {
      .pick_marker_genes(
        df, n_keep = n_keep,
        p_adj_threshold = p_adj_threshold,
        rank_by = rank_by,
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
      show_legend = FALSE,
      gap = grid::unit(0, "mm")
    )

    pal_marker_row <- c(`Most Favorable` = "#2166AC", `Most Adverse` = "#B2182B")
    # Left = favorable; right = adverse.
    strip_l <- rep(NA_character_, nrow(mat_plot))
    strip_l[seq_len(n_fav)] <- "Most Favorable"
    strip_r <- rep(NA_character_, nrow(mat_plot))
    if (n_adv > 0L) {
      strip_r[n_fav + seq_len(n_adv)] <- "Most Adverse"
    }
    marks_at_fav <- seq_len(min(n_mark_labels, n_fav))
    marks_lab_fav <- fav_genes[seq_len(min(n_mark_labels, n_fav))]
    marks_at_adv <- if (n_adv > 0L) {
      n_fav + seq_len(min(n_mark_labels, n_adv))
    } else {
      integer(0)
    }
    marks_lab_adv <- adv_only[seq_len(min(n_mark_labels, n_adv))]

    # Left (favorable): marks leftmost, strip next to heatmap.
    ha_left <- ComplexHeatmap::rowAnnotation(
      marks = ComplexHeatmap::anno_mark(
        at = marks_at_fav,
        labels = marks_lab_fav,
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

    # Right (adverse): strip next to heatmap, gene marks on the outer right.
    ha_right <- NULL
    if (n_adv > 0L) {
      ha_right <- ComplexHeatmap::rowAnnotation(
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_r,
          col = pal_marker_row,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_adv,
          labels = marks_lab_adv,
          side = "right",
          labels_gp = grid::gpar(fontsize = 7),
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(3, 18), c("mm", "mm"))
      )
    }

    row_split_g <- factor(marker_tail, levels = c("Most Favorable", "Most Adverse"))
    hm_col_fun <- .scaled_expr_col_fun_rdgy11(scale_clip)

    ct <- column_title %||% "Global phenotype marker genes (favorable vs adverse)"

    ht <- ComplexHeatmap::Heatmap(
      mat_plot,
      name = "Scaled expr",
      col = hm_col_fun,
      use_raster = use_raster,
      width = heatmap_width_u,
      height = heatmap_height_u,
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
          rank_by = rank_by,
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
      trimws(as.character(gene_info$phenotype_bin)),
      trimws(as.character(gene_info$cell_type)),
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
      show_legend = FALSE,
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
    # Left = favorable strips; right = adverse strips.
    if (length(idx_fav) > 0L) {
      strip_l_pheno[idx_fav] <- "Most Favorable"
      strip_l_ct[idx_fav] <- as.character(gene_info$cell_type[idx_fav])
    }
    if (length(idx_adv) > 0L) {
      strip_r_pheno[idx_adv] <- "Most Adverse"
      strip_r_ct[idx_adv] <- as.character(gene_info$cell_type[idx_adv])
    }

    ha_left <- NULL
    if (length(idx_fav) > 0L) {
      # Left (outside -> heatmap): marks | cell type | phenotype bin
      ha_left <- ComplexHeatmap::rowAnnotation(
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_fav,
          labels = marks_lab_fav,
          side = "left",
          labels_gp = grid::gpar(fontsize = 7),
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        `Cell type` = ComplexHeatmap::anno_simple(
          strip_l_ct,
          col = pal_celltype,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_l_pheno,
          col = pal_group,
          width = grid::unit(3, "mm"),
          na_col = "transparent"
        ),
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(18, 3, 3), c("mm", "mm", "mm"))
      )
    }

    ha_right <- NULL
    if (length(idx_adv) > 0L) {
      # Right (heatmap -> outside): phenotype | cell type | marks
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
          at = marks_at_adv,
          labels = marks_lab_adv,
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

    # Per-column key must match block_key = paste(phenotype_bin, cell_type, sep = "||")
    # (same convention as gene rows) so white boxes align with column blocks.
    ncol_hm <- ncol(mat_plot)
    col_block_key <- rep(NA_character_, ncol_hm)
    for (j in seq_len(ncol_hm)) {
      mid <- meta_idx_hm[j]
      if (!is.na(mid)) {
        col_block_key[j] <- paste(
          trimws(as.character(meta[[group_col]][mid])),
          trimws(as.character(meta[[celltype_col]][mid])),
          sep = "||"
        )
      }
    }
    decorate_ct_rect <- if (isTRUE(outline_marker_blocks)) {
      list(row_split = row_split, col_block_key = col_block_key)
    } else {
      NULL
    }

    ht <- ComplexHeatmap::Heatmap(
      mat_plot,
      name = "phenomap_ct_markers",
      col = hm_col_fun,
      use_raster = use_raster,
      width = heatmap_width_u,
      height = heatmap_height_u,
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
      border = TRUE,
      heatmap_legend_param = list(
        title = "Scaled expr",
        direction = "vertical",
        title_position = "leftcenter-rot",
        legend_height = grid::unit(3, "cm")
      )
    )
  }

  if (isTRUE(draw)) {
    # Manual legends only (HeatmapAnnotation uses show_legend = FALSE so CH does
    # not auto-build a second score legend with its own break inference).
    lgd_at_raw <- if (smin < 0 && smax > 0) {
      c(smin, 0, smax)
    } else if (smax <= 0) {
      c(smin, 0)
    } else if (smin >= 0) {
      c(0, smax)
    } else {
      c(smin, smax)
    }
    lgd_at <- sort(unique(round(lgd_at_raw)))
    if (length(lgd_at) < 2L) {
      lgd_at <- sort(unique(round(c(smin, smax))))
    }
    if (length(lgd_at) < 2L) {
      lgd_at <- c(-1, 1)
    }
    lgd_labels <- as.character(lgd_at)
    lgd_score <- ComplexHeatmap::Legend(
      title = "PhenoMapR score",
      col_fun = score_col_fun,
      at = lgd_at,
      labels = lgd_labels,
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
    if (isTRUE(outline_marker_blocks) && !is.null(decorate_ct_rect)) {
      rect <- .rect_native_ct_marker_blocks(
        decorate_ct_rect$col_block_key,
        decorate_ct_rect$row_split
      )
      rs <- decorate_ct_rect$row_split
      n_slice <- nlevels(rs)
      for (si in seq_len(n_slice)) {
        bk <- levels(rs)[si]
        kk <- match(bk, rect$block)
        if (is.na(kk)) next
        rows_slice <- which(as.integer(rs) == si)
        n_sr <- length(rows_slice)
        if (!n_sr) next
        j1 <- rect$jmin[kk]
        j2 <- rect$jmax[kk]
        ComplexHeatmap::decorate_heatmap_body(
          "phenomap_ct_markers",
          {
            grid::grid.rect(
              x = grid::unit(j1 - 1L, "native"),
              y = grid::unit(1, "npc"),
              width = grid::unit(j2 - j1 + 1L, "native"),
              height = grid::unit(1, "npc"),
              hjust = 0,
              vjust = 1,
              gp = grid::gpar(col = "white", fill = NA, lty = 1L, lwd = 1)
            )
          },
          slice = si
        )
      }
    }
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

#' Coerce \code{meta[[score_col]]} to a length-\code{nrow(meta)} numeric vector.
#' Handles 1-column matrices (e.g. \code{scale()} bound into a data.frame).
#'
#' @noRd
#' @keywords internal
.meta_score_numeric_vector <- function(meta, score_col) {
  v <- meta[[score_col]]
  n <- nrow(meta)
  if (is.null(v)) {
    return(rep(NA_real_, n))
  }
  if (is.matrix(v)) {
    if (nrow(v) == n && ncol(v) >= 1L) {
      v <- as.numeric(v[, 1L])
    } else if (ncol(v) == n && nrow(v) >= 1L) {
      v <- as.numeric(v[1L, ])
    } else {
      v <- suppressWarnings(as.numeric(v))
      if (length(v) != n) {
        v <- rep(NA_real_, n)
      }
    }
  } else {
    v <- suppressWarnings(as.numeric(v))
  }
  if (length(v) != n) {
    if (length(v) == 1L) {
      v <- rep(v, n)
    } else if (length(v) > n) {
      v <- v[seq_len(n)]
    } else {
      v <- c(v, rep(NA_real_, n - length(v)))
    }
  }
  v
}

#' PhenoMapR scores aligned to heatmap column order (\code{meta_idx_hm}).
#'
#' @noRd
#' @keywords internal
.score_ann_for_heatmap <- function(meta, score_col, meta_idx_hm) {
  scv <- .meta_score_numeric_vector(meta, score_col)
  out <- rep(NA_real_, length(meta_idx_hm))
  ok <- !is.na(meta_idx_hm) & meta_idx_hm >= 1L & meta_idx_hm <= length(scv)
  out[ok] <- scv[as.integer(meta_idx_hm[ok])]
  out
}

#' Diverging PhenoMapR score colors: blue (\eqn{< 0}), white at 0, red (\eqn{> 0}).
#'
#' Breakpoints use \code{smin} and \code{smax} supplied by the caller (typically
#' the range of scores on heatmap columns). When the range crosses zero,
#' breakpoints are \code{c(smin, 0, smax)}.
#'
#' @noRd
#' @keywords internal
.phenomap_score_col_fun <- function(smin, smax) {
  if (!is.finite(smin) || !is.finite(smax) || smin == smax) {
    smin <- -1
    smax <- 1
  }
  if (smin < 0 && smax > 0) {
    circlize::colorRamp2(
      c(smin, 0, smax),
      c("#2166AC", "#FFFFFF", "#B2182B")
    )
  } else if (smax <= 0) {
    circlize::colorRamp2(c(smin, 0), c("#2166AC", "#FFFFFF"))
  } else if (smin >= 0) {
    circlize::colorRamp2(c(0, smax), c("#FFFFFF", "#B2182B"))
  } else {
    circlize::colorRamp2(c(smin, smax), c("#2166AC", "#B2182B"))
  }
}

#' Column and row spans for white boxes around each cell-type-specific marker block.
#'
#' One row per \code{levels(row_split)} with a matching \code{col_block_key};
#' \code{jmin}/\code{jmax} are 1-based column indices (heatmap body native width
#' units). \code{r1}/\code{r2} are global matrix row indices for that block.
#'
#' @noRd
#' @keywords internal
.rect_native_ct_marker_blocks <- function(col_block_key, row_split) {
  empty <- data.frame(
    block = character(),
    jmin = integer(),
    jmax = integer(),
    r1 = integer(),
    r2 = integer(),
    stringsAsFactors = FALSE
  )
  columns <- as.character(col_block_key)
  nc <- length(columns)
  nr <- length(row_split)
  if (nc < 1L || nr < 1L) {
    return(empty)
  }
  rows <- factor(as.character(row_split), levels = levels(row_split))
  lvls <- levels(rows)
  blocks <- character()
  jmin <- jmax <- r1 <- r2 <- integer()
  for (lvl in lvls) {
    idx_rows <- which(rows == lvl)
    if (!length(idx_rows)) next
    jj <- which(!is.na(columns) & columns == lvl)
    if (!length(jj)) next
    blocks <- c(blocks, lvl)
    jmin <- c(jmin, min(jj))
    jmax <- c(jmax, max(jj))
    r1 <- c(r1, min(idx_rows))
    r2 <- c(r2, max(idx_rows))
  }
  if (!length(blocks)) {
    return(empty)
  }
  data.frame(
    block = blocks,
    jmin = jmin,
    jmax = jmax,
    r1 = r1,
    r2 = r2,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
.pick_marker_genes <- function(df,
                               n_keep,
                               p_adj_threshold,
                               rank_by = c("lfc", "p_adj"),
                               valid_genes = NULL) {
  req <- c("gene", "avg_log2FC", "p_adj")
  if (is.null(df) || nrow(df) == 0L || !all(req %in% names(df))) {
    return(character(0))
  }
  # Default: keep significant positive-effect markers; ordering controlled by rank_by.
  df <- df[is.finite(df$avg_log2FC) & is.finite(df$p_adj), , drop = FALSE]
  df <- df[df$p_adj < p_adj_threshold & df$avg_log2FC > 0, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(character(0))
  }
  rank_by <- match.arg(rank_by)
  if (rank_by == "p_adj") {
    # When ranking by p_adj, enforce a stronger effect-size cutoff.
    df <- df[df$avg_log2FC > 1, , drop = FALSE]
    if (nrow(df) == 0L) return(character(0))
    df <- df[order(df$p_adj, -df$avg_log2FC), , drop = FALSE]
  } else {
    # LFC-preferred: pick the largest positive LFC markers among significant genes.
    df <- df[order(-df$avg_log2FC, df$p_adj), , drop = FALSE]
  }
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
  score_vec <- .meta_score_numeric_vector(meta, score_col)
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
  score_vec <- .meta_score_numeric_vector(meta, score_col)
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
