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
#'   cell-type colors come from \code{color_schemes$celltype} (default:
#'   \code{\link{get_celltype_palette}()}).
#' @param color_schemes Optional named list controlling heatmap colors. Each
#'   element may be \code{"default"}, a shorthand string like
#'   \code{"brewer:RdBu"} or \code{"viridis:plasma"}, or a list with
#'   \code{source} (\code{"default"}, \code{"brewer"}, \code{"viridis"}, or
#'   \code{"manual"}) plus \code{name} and/or \code{colors}. Supported names:
#'   \code{phenotype} (discrete groups), \code{score} (PhenoMapR score bar),
#'   \code{expression} (scaled-expression heatmap body), and \code{celltype}.
#'   See \code{\link{list_marker_heatmap_color_palettes}()} for built-in options.
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
#' @param mark_label_fontsize Font size (points) for \code{anno_mark} gene labels.
#'   Default \code{7}. Decrease when many labels are drawn in a compact heatmap.
#' @param color_mark_labels_by_celltype If \code{TRUE}, \code{anno_mark} gene labels
#'   are colored with \code{celltype_palette} according to each gene's cell type.
#'   Applies to \code{heatmap_type = "cell_type_specific"} (and to \code{"global"}
#'   when marker tables include a \code{cell_type} column). Default \code{FALSE}
#'   keeps labels black.
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
#'   \code{heatmap_type = "cell_type_specific"}, draws outline boxes around each
#'   marker-gene block (one per cell-type \eqn{\times} phenotype-bin slice)
#'   after \code{draw()}. Set \code{FALSE} to omit them.
#' @param block_outline_color Outline colour used when
#'   \code{outline_marker_blocks = TRUE} (default \code{"white"}). Set to
#'   \code{"black"} for stronger visual separation between cell-type
#'   \eqn{\times} phenotype-bin blocks.
#' @param block_outline_lwd Outline line width used when
#'   \code{outline_marker_blocks = TRUE} (default \code{1}). Slightly thicker
#'   widths (e.g. \code{1.5} or \code{2}) read better with darker
#'   \code{block_outline_color}.
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
                                   color_schemes = NULL,
                                   heatmap_type = c("global", "cell_type_specific"),
                                   top_n_markers = 20L,
                                   rank_by = c("lfc", "p_adj"),
                                   n_mark_labels = 5L,
                                   mark_label_fontsize = 7,
                                   color_mark_labels_by_celltype = FALSE,
                                   p_adj_threshold = 0.05,
                                   scale_clip = NULL,
                                   heatmap_width = NULL,
                                   heatmap_height = NULL,
                                   column_title = NULL,
                                   draw = TRUE,
                                   use_raster = FALSE,
                                   outline_marker_blocks = TRUE,
                                   block_outline_color = "white",
                                   block_outline_lwd = 1) {
  heatmap_type <- match.arg(heatmap_type)
  rank_by <- match.arg(rank_by)
  if (!is.character(block_outline_color) ||
      length(block_outline_color) != 1L ||
      is.na(block_outline_color) ||
      !nzchar(block_outline_color)) {
    block_outline_color <- "white"
  }
  if (!is.numeric(block_outline_lwd) ||
      length(block_outline_lwd) != 1L ||
      !is.finite(block_outline_lwd) ||
      block_outline_lwd <= 0) {
    block_outline_lwd <- 1
  }

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
  mark_label_fontsize <- as.numeric(mark_label_fontsize)[1L]
  if (!is.finite(mark_label_fontsize) || mark_label_fontsize <= 0) {
    mark_label_fontsize <- 7
  }
  color_mark_labels_by_celltype <- isTRUE(color_mark_labels_by_celltype)

  hm_draw_marks_at <- integer(0)
  hm_draw_n_rows <- 0L
  hm_mark_anno_width_mm <- 18
  decorate_pairs <- list()

  # Strip thickness: row annotations use `width`, column annotations
  # use `height` (ComplexHeatmap convention). Previously the top
  # strips passed `width`, which column anno_simple ignores, so they
  # fell back to ht_opt$simple_anno_size (5 mm) and looked ~2x thicker
  # than the 3 mm row strips.
  strip_mm <- 3
  row_anno_strip_u  <- grid::unit(strip_mm, "mm")
  col_anno_strip_u  <- grid::unit(strip_mm, "mm")
  col_anno_name_gp  <- grid::gpar(fontsize = 7)

  # celltype_col is mandatory for cell-type-specific heatmaps (every block
  # is keyed by a cell type) but only decorative for the global heatmap (used
  # for the top "Cell type" annotation strip). Treat NULL / "" / a column not
  # present in `meta` the same way and gracefully drop the cell-type strip
  # when running in global mode.
  has_celltype <- !is.null(celltype_col) && nzchar(celltype_col) &&
    celltype_col %in% names(meta)
  if (heatmap_type == "cell_type_specific" && !has_celltype) {
    stop("'celltype_col' must name an existing column of 'meta' when heatmap_type = 'cell_type_specific'.")
  }
  req_meta <- c(cell_id_col, group_col, score_col)
  if (has_celltype) req_meta <- c(req_meta, celltype_col)
  if (!all(req_meta %in% names(meta))) {
    stop("meta must contain columns: ", paste(req_meta, collapse = ", "))
  }

  gene_ids <- rownames(expr_mat)
  if (is.null(gene_ids)) {
    stop("'expr_mat' must have row names (gene symbols).")
  }

  hm_group_levels <- c("Most Favorable", "Other", "Most Adverse")
  hm_celltype_levels <- if (has_celltype) {
    # Drop ghost factor levels (levels declared on the factor but with
    # zero rows in the actual data). Without this, a `cell_type` factor
    # whose levels are c("1", "Acinar", "Beta") but whose values are
    # all "Acinar" still reports 3 levels, which (a) wastes a slot in
    # the cell-type legend and (b) re-introduces the very "1" tile
    # this guard is meant to prevent.
    keep <- as.character(unique(stats::na.omit(meta[[celltype_col]])))
    intersect(levels(factor(meta[[celltype_col]])), keep)
  } else {
    character(0)
  }
  # `show_celltype_strip` controls whether the *visual* cell-type strip
  # annotation (and its associated legend tile) is rendered on the
  # heatmap. It is independent of `has_celltype`, which still gates
  # the cell-type-specific *slicing* logic (the row/column splits use
  # block_key = paste(phenotype_bin, cell_type, ...) and need the
  # column to exist even when there is only one level).
  #
  # A uniform colored strip carrying only one value is uninformative
  # AND historically leaked a stray auto-legend tile -- a tiny "1" or
  # single-celltype badge floating next to the main legend column --
  # despite `show_legend = FALSE` on the parent HeatmapAnnotation /
  # rowAnnotation in some ComplexHeatmap versions. Skipping the
  # strip entirely (in BOTH global and cell_type_specific modes) is
  # the most robust way to suppress that artifact: no anno_simple
  # for cell type means no auto-legend to leak.
  show_celltype_strip <- has_celltype && length(hm_celltype_levels) > 1L
  # Preserve the long-standing global-mode behaviour: when the cell
  # type has a single level, the global heatmap also drops the
  # has_celltype gate so legend code paths that take the "no cell
  # type at all" branch keep firing. Cell-type-specific keeps
  # has_celltype = TRUE so block_key splitting still works.
  if (has_celltype && heatmap_type == "global" &&
      length(hm_celltype_levels) <= 1L) {
    has_celltype <- FALSE
    hm_celltype_levels <- character(0)
    show_celltype_strip <- FALSE
  }

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

  schemes <- .normalize_marker_heatmap_color_schemes(color_schemes)
  pal_group <- .resolve_phenotype_palette(schemes$phenotype)
  score_colors <- .resolve_score_colors(schemes$score)
  expr_colors <- .resolve_expression_colors(schemes$expression, n = 11L)

  if (has_celltype) {
    pal_celltype <- .resolve_celltype_palette_from_spec(
      schemes$celltype,
      hm_celltype_levels,
      override = celltype_palette
    )
    pal_celltype[is.na(pal_celltype)] <- "#BBBBBB"
  } else {
    pal_celltype <- character(0)
  }

  # PhenoMapR score bar: diverging blue-white-red with white at 0. Color limits
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
  score_col_fun <- .phenomap_score_col_fun_from_colors(smin, smax, score_colors)

  # Initialise the (row_slice, column_slice) decoration pair list to
  # empty so the post-render block-outline loop is a no-op for the
  # "global" heatmap path (which has no per-cell-type splits to
  # outline). The cell-type-specific branch overwrites this with
  # diagonal pairs.
  decorate_pairs <- list()

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

    # Top stack (CH: first = nearest heatmap): phenotype -> (cell type) -> score.
    # Cell type strip is omitted when the user hasn't mapped a cell-type column.
    top_anno_args <- list(
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        height = col_anno_strip_u
      )
    )
    if (show_celltype_strip) {
      # Same single-level guard as the cell_type_specific path:
      # only emit the cell-type strip when there are at least two
      # represented levels. A uniform-coloured strip is
      # uninformative AND historically leaked a stray "1" tile.
      top_anno_args[["Cell type"]] <- ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        height = col_anno_strip_u
      )
    }
    top_anno_args[["PhenoMapR score"]] <- ComplexHeatmap::anno_simple(
      score_ann,
      col = score_col_fun,
      height = col_anno_strip_u
    )
    ha_top <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(top_anno_args, list(
        annotation_name_side = "right",
        annotation_name_gp = col_anno_name_gp,
        show_annotation_name = TRUE,
        show_legend = FALSE,
        gap = grid::unit(0, "mm")
      ))
    )

    pal_marker_row <- pal_group[c("Most Favorable", "Most Adverse")]
    # Left = favorable; right = adverse.
    strip_l <- rep(NA_character_, nrow(mat_plot))
    if (n_fav > 0L) {
      strip_l[seq_len(n_fav)] <- "Most Favorable"
    }
    strip_r <- rep(NA_character_, nrow(mat_plot))
    if (n_adv > 0L) {
      strip_r[n_fav + seq_len(n_adv)] <- "Most Adverse"
    }
    marks_at_fav <- if (n_fav > 0L) seq_len(min(n_mark_labels, n_fav)) else integer(0)
    marks_lab_fav <- if (n_fav > 0L) fav_genes[marks_at_fav] else character(0)
    marks_at_adv <- if (n_adv > 0L) {
      n_fav + seq_len(min(n_mark_labels, n_adv))
    } else {
      integer(0)
    }
    marks_lab_adv <- if (n_adv > 0L) {
      adv_only[seq_len(min(n_mark_labels, n_adv))]
    } else {
      character(0)
    }
    marks_ct_fav <- if (color_mark_labels_by_celltype && length(marks_lab_fav) > 0L) {
      .lookup_mark_gene_celltypes(marks_lab_fav, favorable_df)
    } else {
      character(0)
    }
    marks_ct_adv <- if (color_mark_labels_by_celltype && length(marks_lab_adv) > 0L) {
      .lookup_mark_gene_celltypes(marks_lab_adv, adverse_df)
    } else {
      character(0)
    }
    mark_labels_gp_fav <- .anno_mark_labels_gp(
      marks_ct_fav, color_mark_labels_by_celltype, pal_celltype,
      fontsize = mark_label_fontsize
    )
    mark_labels_gp_adv <- .anno_mark_labels_gp(
      marks_ct_adv, color_mark_labels_by_celltype, pal_celltype,
      fontsize = mark_label_fontsize
    )
    all_mark_labels <- c(marks_lab_fav, marks_lab_adv)
    mark_anno_width_mm <- .phenomap_mark_anno_width_mm(
      all_mark_labels, mark_label_fontsize
    )
    hm_draw_marks_at <- c(marks_at_fav, marks_at_adv)
    hm_draw_n_rows <- nrow(mat_plot)
    hm_mark_anno_width_mm <- mark_anno_width_mm

    # Left (favorable): marks leftmost, strip next to heatmap. Only built when
    # we actually have favorable rows -- `anno_mark(at = integer(0))` triggers
    # a "select less than one element" error inside ComplexHeatmap.
    ha_left <- NULL
    if (n_fav > 0L) {
      ha_left <- ComplexHeatmap::rowAnnotation(
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_fav,
          labels = marks_lab_fav,
          side = "left",
          labels_gp = mark_labels_gp_fav,
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_l,
          col = pal_marker_row,
          width = row_anno_strip_u,
          na_col = "transparent"
        ),
        # Suppress auto-legends from this rowAnnotation. The manual
        # phenotype-group / cell-type legends below (lgd_group +
        # optional cell-type Legend) are the canonical ones we want
        # rendered. Without this, anno_simple() auto-emits a stray
        # legend per strip -- when a strip has only one level (e.g.
        # the metadata cell-type column has a single value like "1")
        # the auto-legend renders as a lone teal "1" swatch floating
        # next to the Scaled expr legend.
        show_legend = FALSE,
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(mark_anno_width_mm, 3), c("mm", "mm"))
      )
    }

    # Right (adverse): strip next to heatmap, gene marks on the outer right.
    ha_right <- NULL
    if (n_adv > 0L) {
      ha_right <- ComplexHeatmap::rowAnnotation(
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_r,
          col = pal_marker_row,
          width = row_anno_strip_u,
          na_col = "transparent"
        ),
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_adv,
          labels = marks_lab_adv,
          side = "right",
          labels_gp = mark_labels_gp_adv,
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        ),
        # Same auto-legend suppression as ha_left -- see comment there.
        show_legend = FALSE,
        show_annotation_name = FALSE,
        gap = grid::unit(0, "mm"),
        annotation_width = grid::unit(c(3, mark_anno_width_mm), c("mm", "mm"))
      )
    }

    # Drop empty tail levels so ComplexHeatmap doesn't try to slice a row block
    # with zero rows (another way to surface "select less than one element").
    present_tails <- intersect(c("Most Favorable", "Most Adverse"), unique(marker_tail))
    row_split_g <- factor(marker_tail, levels = present_tails)
    hm_col_fun <- .scaled_expr_col_fun_from_colors(scale_clip, expr_colors)

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
    marks_ct_fav <- character(0)
    marks_at_adv <- integer(0)
    marks_lab_adv <- character(0)
    marks_ct_adv <- character(0)
    start <- 1L
    for (bk in unique(block_key)) {
      ii <- which(block_key == bk)
      n <- length(ii)
      if (n == 0L) next
      nm <- min(n_mark_labels, n)
      at <- start + seq_len(nm) - 1L
      lab <- as.character(gene_info$gene[ii[seq_len(nm)]])
      lab_ct <- as.character(gene_info$cell_type[ii[seq_len(nm)]])
      bin <- as.character(gene_info$phenotype_bin[ii[1L]])
      if (identical(bin, "Most Favorable")) {
        marks_at_fav <- c(marks_at_fav, at)
        marks_lab_fav <- c(marks_lab_fav, lab)
        marks_ct_fav <- c(marks_ct_fav, lab_ct)
      } else if (identical(bin, "Most Adverse")) {
        marks_at_adv <- c(marks_at_adv, at)
        marks_lab_adv <- c(marks_lab_adv, lab)
        marks_ct_adv <- c(marks_ct_adv, lab_ct)
      }
      start <- start + n
    }
    mark_labels_gp_fav <- .anno_mark_labels_gp(
      marks_ct_fav, color_mark_labels_by_celltype, pal_celltype,
      fontsize = mark_label_fontsize
    )
    mark_labels_gp_adv <- .anno_mark_labels_gp(
      marks_ct_adv, color_mark_labels_by_celltype, pal_celltype,
      fontsize = mark_label_fontsize
    )
    all_mark_labels <- c(marks_lab_fav, marks_lab_adv)
    mark_anno_width_mm <- .phenomap_mark_anno_width_mm(
      all_mark_labels, mark_label_fontsize
    )
    hm_draw_marks_at <- c(marks_at_fav, marks_at_adv)
    hm_draw_n_rows <- nrow(mat_plot)
    hm_mark_anno_width_mm <- mark_anno_width_mm

    # Build the top annotation as a named list so the cell-type
    # strip can be inserted conditionally. When the cell-type
    # column has just one represented level, dropping the strip
    # entirely (rather than rendering a uniform colored bar) is
    # what prevents the stray "1" auto-legend tile that some
    # ComplexHeatmap versions emit despite `show_legend = FALSE`.
    top_anno_args_ct <- list(
      `Phenotype group` = ComplexHeatmap::anno_simple(
        as.character(meta[[group_col]][meta_idx_hm]),
        col = pal_group,
        height = col_anno_strip_u
      )
    )
    if (show_celltype_strip) {
      top_anno_args_ct[["Cell type"]] <- ComplexHeatmap::anno_simple(
        as.character(meta[[celltype_col]][meta_idx_hm]),
        col = pal_celltype,
        height = col_anno_strip_u
      )
    }
    top_anno_args_ct[["PhenoMapR score"]] <- ComplexHeatmap::anno_simple(
      score_ann,
      col = score_col_fun,
      height = col_anno_strip_u
    )
    ha_top <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(top_anno_args_ct, list(
        annotation_name_side = "right",
        annotation_name_gp = col_anno_name_gp,
        show_annotation_name = TRUE,
        show_legend = FALSE,
        gap = grid::unit(0, "mm")
      ))
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

    # Build left/right rowAnnotations as named lists so the cell-type
    # strip can be omitted whenever `show_celltype_strip` is FALSE
    # (single represented cell-type level). Dropping the strip is
    # the only reliable way to suppress the stray teal "1" tile that
    # ComplexHeatmap can emit despite show_legend = FALSE on the
    # parent rowAnnotation -- see the long comment at the top of the
    # function (`show_celltype_strip` definition) for the full story.
    ha_left <- NULL
    if (length(idx_fav) > 0L) {
      left_args <- list(
        # `which = "row"` must be set explicitly here because we
        # construct the rowAnnotation via do.call() below, which
        # bypasses ComplexHeatmap's normal auto-detection of the
        # annotation orientation from the enclosing call.
        marks = ComplexHeatmap::anno_mark(
          at = marks_at_fav,
          labels = marks_lab_fav,
          which = "row",
          side = "left",
          labels_gp = mark_labels_gp_fav,
          link_gp = grid::gpar(col = "grey50", lwd = 0.6),
          padding = grid::unit(0.5, "mm")
        )
      )
      if (show_celltype_strip) {
        # `which = "row"` is REQUIRED on every anno_simple constructed
        # for a rowAnnotation built via do.call(): the auto-detection
        # of orientation only fires when the annotations are passed
        # to rowAnnotation()'s named-args call (via `match.call()`
        # magic). do.call() bypasses that, so we must mark each
        # annotation as a row annotation up-front.
        left_args[["Cell type"]] <- ComplexHeatmap::anno_simple(
          strip_l_ct,
          col = pal_celltype,
          which = "row",
          width = row_anno_strip_u,
          na_col = "transparent"
        )
      }
      left_args[["Phenotype"]] <- ComplexHeatmap::anno_simple(
        strip_l_pheno,
        col = pal_group,
        which = "row",
        width = row_anno_strip_u,
        na_col = "transparent"
      )
      # Annotation widths track the actual entries in `left_args`:
      # marks (18mm) is always present; cell-type strip is only
      # included when show_celltype_strip is TRUE; phenotype strip
      # is always present.
      left_widths <- if (show_celltype_strip) {
        grid::unit(c(mark_anno_width_mm, 3, 3), c("mm", "mm", "mm"))
      } else {
        grid::unit(c(mark_anno_width_mm, 3), c("mm", "mm"))
      }
      ha_left <- do.call(
        ComplexHeatmap::rowAnnotation,
        c(left_args, list(
          show_legend = FALSE,
          show_annotation_name = FALSE,
          gap = grid::unit(0, "mm"),
          annotation_width = left_widths
        ))
      )
    }

    ha_right <- NULL
    if (length(idx_adv) > 0L) {
      right_args <- list(
        `Phenotype` = ComplexHeatmap::anno_simple(
          strip_r_pheno,
          col = pal_group,
          # See the long note in ha_left's anno_simple block for why
          # `which = "row"` is required on every anno_simple
          # constructed for a do.call(rowAnnotation, ...) build.
          which = "row",
          width = row_anno_strip_u,
          na_col = "transparent"
        )
      )
      if (show_celltype_strip) {
        right_args[["Cell type"]] <- ComplexHeatmap::anno_simple(
          strip_r_ct,
          col = pal_celltype,
          which = "row",
          width = row_anno_strip_u,
          na_col = "transparent"
        )
      }
      right_args[["marks"]] <- ComplexHeatmap::anno_mark(
        at = marks_at_adv,
        labels = marks_lab_adv,
        # See note in ha_left: do.call(rowAnnotation, ...) does
        # not propagate the row-orientation default to anno_mark,
        # so we set it explicitly.
        which = "row",
        side = "right",
        labels_gp = mark_labels_gp_adv,
        link_gp = grid::gpar(col = "grey50", lwd = 0.6),
        padding = grid::unit(0.5, "mm")
      )
      right_widths <- if (show_celltype_strip) {
        grid::unit(c(3, 3, mark_anno_width_mm), c("mm", "mm", "mm"))
      } else {
        grid::unit(c(3, mark_anno_width_mm), c("mm", "mm"))
      }
      ha_right <- do.call(
        ComplexHeatmap::rowAnnotation,
        c(right_args, list(
          show_legend = FALSE,
          show_annotation_name = FALSE,
          gap = grid::unit(0, "mm"),
          annotation_width = right_widths
        ))
      )
    }

    hm_col_fun <- .scaled_expr_col_fun_from_colors(scale_clip, expr_colors)

    ct <- column_title %||% "Cell-type-specific phenotype marker genes"

    # Per-column key must match block_key = paste(phenotype_bin, cell_type, sep = "||")
    # (same convention as gene rows) so block outlines align with column blocks.
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
    # Build a column_split factor that mirrors the row_split block_key
    # scheme. By splitting BOTH axes on the same (phenotype_bin x
    # cell_type) key, ComplexHeatmap renders each block as an
    # independent (row_slice, column_slice) viewport. We can then
    # outline a block by simply calling
    # `decorate_heatmap_body(row_slice = ri, column_slice = cj, ...)`
    # with `grid.rect()` and no explicit coordinates -- following the
    # EcoTyper / jokergoo recommended pattern. This is far more
    # robust than the previous "decorate_heatmap_body + manual native
    # x/width" approach, which failed silently because the inner
    # pushViewport reset the native x-scale to 0..1 and the rectangle
    # ended up clipped outside the visible body.
    #
    # NA columns (cells without a phenotype/celltype label) get
    # bucketed into their own slice so they never trigger an
    # incorrect outline; we filter the decorate loop below to only
    # the slices whose level matches the row_split levels.
    col_split_levels <- unique(col_block_key[!is.na(col_block_key)])
    has_unassigned <- any(is.na(col_block_key))
    col_split_factor <- if (length(col_split_levels) > 0L) {
      vals <- col_block_key
      lvls <- col_split_levels
      if (has_unassigned) {
        vals <- ifelse(is.na(vals), "__phenomap_unassigned__", vals)
        lvls <- c(col_split_levels, "__phenomap_unassigned__")
      }
      factor(vals, levels = lvls)
    } else {
      NULL
    }
    # Match row-slice indices to column-slice indices on the shared
    # block_key. Diagonal pairs (i, j) where row_levels[i] ==
    # column_levels[j] are the only ones we outline.
    row_slice_levels <- levels(row_split)
    decorate_pairs <- if (isTRUE(outline_marker_blocks) &&
                         !is.null(col_split_factor) &&
                         length(row_slice_levels) > 0L) {
      col_lv <- levels(col_split_factor)
      pairs <- list()
      for (ri in seq_along(row_slice_levels)) {
        cj <- match(row_slice_levels[ri], col_lv)
        if (!is.na(cj)) {
          pairs[[length(pairs) + 1L]] <- c(row_slice = ri, column_slice = cj)
        }
      }
      pairs
    } else {
      list()
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
      column_split = col_split_factor,
      cluster_column_slices = FALSE,
      column_gap = grid::unit(0, "mm"),
      # Per-slice column titles are noisy (one block_key per slice
      # repeated above the heatmap), so suppress them. The overall
      # column title is set below via column_title = ct (length 1)
      # which ComplexHeatmap centres across all slices when split.
      column_title = ct,
      column_title_side = "top",
      show_column_names = FALSE,
      show_row_names = FALSE,
      top_annotation = ha_top,
      left_annotation = ha_left,
      right_annotation = ha_right,
      # Avoid black borders around every row-split slice and/or cell.
      border = FALSE,
      rect_gp = grid::gpar(col = NA),
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
    # The underlying data column still holds "Most Favorable" /
    # "Most Adverse" (so every downstream pipeline that pattern-
    # matches on those exact strings keeps working). For the
    # rendered legend we relabel them to the more neutral
    # "Most Phenotype +" / "Most Phenotype -" so the heatmap reads
    # the same regardless of whether higher-z = better or higher-z
    # = worse for the cohort under study. Passing `labels` while
    # leaving `at` (and the colour mapping) untouched is the
    # ComplexHeatmap-supported way to do this without re-encoding
    # the metadata.
    #
    # Convention (must match the rest of the package, including the
    # Shiny app and downstream score-by-cell-type plots):
    #   "Most Adverse"   (red,  #B2182B) -> "Most Phenotype +"
    #   "Most Favorable" (blue, #2166AC) -> "Most Phenotype -"
    # i.e. higher PhenoMapR z-scores = "Most Phenotype +" = red,
    # lower z-scores = "Most Phenotype -" = blue. The labels vector
    # below is paired position-by-position with `at`, so flipping
    # the labels (NOT the `at` order or the `pal_group` fill order)
    # is what gets the +/- assigned to the correct colours.
    lgd_group <- ComplexHeatmap::Legend(
      title = "Phenotype group",
      at = c("Most Favorable", "Other", "Most Adverse"),
      labels = c("Most Phenotype -", "Other", "Most Phenotype +"),
      legend_gp = grid::gpar(
        fill = pal_group[c("Most Favorable", "Other", "Most Adverse")]
      )
    )
    annotation_legend_list <- list(lgd_score, lgd_group)
    # Only include the Cell type legend when there are at least 2 distinct
    # cell types -- a 1-entry legend is just a stray colored swatch with no
    # discriminative value.
    if (has_celltype && length(hm_celltype_levels) > 1L) {
      annotation_legend_list <- c(
        annotation_legend_list,
        list(ComplexHeatmap::Legend(
          title = "Cell type",
          at = hm_celltype_levels,
          legend_gp = grid::gpar(fill = pal_celltype[hm_celltype_levels])
        ))
      )
    }
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      show_annotation_legend = TRUE,
      annotation_legend_list = annotation_legend_list,
      merge_legends = TRUE,
      padding = grid::unit(
        .phenomap_draw_padding_mm(
          n_rows = hm_draw_n_rows,
          marks_at = hm_draw_marks_at,
          fontsize = mark_label_fontsize,
          mark_anno_width_mm = hm_mark_anno_width_mm,
          heatmap_type = heatmap_type
        ),
        "mm"
      )
    )
    if (isTRUE(outline_marker_blocks) && length(decorate_pairs) > 0L) {
      # EcoTyper-style block outlines. Because we column_split on the
      # SAME (phenotype_bin x cell_type) key as the row split,
      # ComplexHeatmap renders each block as its own
      # (row_slice, column_slice) viewport. Calling grid.rect() with
      # no explicit coordinates inside `decorate_heatmap_body(..., row_slice = ri,
      # column_slice = cj, ...)` therefore draws a rectangle that
      # exactly fills that block's viewport -- no fragile native-unit
      # arithmetic, no nested viewport pushes that reset the x-scale.
      # See the EcoTyper marker-gene heatmap and the ComplexHeatmap
      # FAQ ("Draw border on one axis only" / "split on columns?")
      # for the canonical version of this idiom.
      lwd_use <- block_outline_lwd
      col_use <- block_outline_color
      for (pair in decorate_pairs) {
        ri <- pair[["row_slice"]]
        cj <- pair[["column_slice"]]
        ComplexHeatmap::decorate_heatmap_body(
          "phenomap_ct_markers",
          row_slice = ri,
          column_slice = cj,
          {
            grid::grid.rect(
              gp = grid::gpar(
                col = col_use,
                fill = NA,
                lty = 1L,
                lwd = lwd_use
              )
            )
          }
        )
      }
    }
  }
  invisible(ht)
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
.lookup_mark_gene_celltypes <- function(genes, df) {
  if (is.null(genes) || length(genes) == 0L || is.null(df) ||
      !"cell_type" %in% names(df)) {
    return(rep(NA_character_, length(genes)))
  }
  as.character(df$cell_type[match(genes, df$gene)])
}


#' Build `labels_gp` for ComplexHeatmap::anno_mark().
#'
#' When `color_by_celltype` is TRUE and `cell_types` is parallel to the
#' labels, text is colored from `pal_celltype`; otherwise labels use
#' `default_col`.
#' @keywords internal
.phenomap_mark_anno_width_mm <- function(labels,
                                         fontsize = 7,
                                         min_mm = 18,
                                         max_mm = 44) {
  if (is.null(labels) || length(labels) == 0L) {
    return(min_mm)
  }
  nch <- max(nchar(as.character(labels)), na.rm = TRUE)
  if (!is.finite(nch) || nch < 1L) {
    return(min_mm)
  }
  mm_per_char <- 0.28 * (fontsize / 7)
  w <- nch * mm_per_char + 5
  max(min_mm, min(max_mm, w))
}


#' Default heatmap body height (mm) from row count.
#' @keywords internal
.phenomap_hm_height_mm <- function(n_rows, heatmap_type, n_splits = 1L) {
  if (n_rows <= 0L) {
    return(if (heatmap_type == "global") 42 else 58)
  }
  row_mm <- if (heatmap_type == "global") 3.2 else 3.8
  h <- n_rows * row_mm
  min_h <- if (heatmap_type == "global") 42 else 58
  max(h, min_h)
}


#' Top/right/bottom/left draw() padding (mm) for anno_mark label headroom.
#' Kept modest so draw() margins do not blow up the layout; prefer
#' `heatmap_height` / knitr `fig.height` for tall matrices.
#' @keywords internal
.phenomap_draw_padding_mm <- function(n_rows,
                                      marks_at,
                                      fontsize = 7,
                                      mark_anno_width_mm = 18,
                                      heatmap_type = c("global", "cell_type_specific")) {
  heatmap_type <- match.arg(heatmap_type)
  top <- 6
  right <- 50
  bottom <- 6
  left <- 4

  if (identical(heatmap_type, "cell_type_specific")) {
    top <- top + 2
    bottom <- bottom + 2
  }

  if (length(marks_at) > 0L && n_rows > 0L) {
    marks_at <- marks_at[marks_at >= 1L & marks_at <= n_rows]
    if (length(marks_at) > 0L) {
      if (any(marks_at <= 2L)) {
        top <- top + 3
      }
      if (any(marks_at >= n_rows - 1L)) {
        bottom <- bottom + 3
      }
    }
  }

  c(top, right, bottom, left)
}


#' Build `labels_gp` for ComplexHeatmap::anno_mark().
#'
#' When `color_by_celltype` is TRUE and `cell_types` is parallel to the
#' labels, text is colored from `pal_celltype`; otherwise labels use
#' `default_col`.
#' @keywords internal
.anno_mark_labels_gp <- function(cell_types,
                                 color_by_celltype,
                                 pal_celltype,
                                 fontsize = 7,
                                 default_col = "black") {
  if (!isTRUE(color_by_celltype)) {
    return(grid::gpar(fontsize = fontsize, col = default_col))
  }
  if (is.null(cell_types) || length(cell_types) == 0L) {
    return(grid::gpar(fontsize = fontsize, col = default_col))
  }
  ct <- trimws(as.character(cell_types))
  cols <- unname(pal_celltype[ct])
  missing <- is.na(cols) | !nzchar(cols)
  if (any(missing)) {
    cols[missing] <- default_col
  }
  grid::gpar(fontsize = fontsize, col = cols)
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
