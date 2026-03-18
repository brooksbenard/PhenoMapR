#' Plot reference signature z-scores as a 1×N heatmap
#'
#' Takes the gene z-score data.frame returned by \code{\link{derive_reference_from_bulk}}
#' and draws a single-row heatmap with genes ordered by z-score, top/bottom gene
#' labels, and a strip indicating \code{|z| > 2}. Requires the \pkg{ComplexHeatmap}
#' and \pkg{circlize} packages.
#'
#' @param reference Data.frame with genes as rownames and a single numeric column
#'   of z-scores (e.g. the return value of \code{\link{derive_reference_from_bulk}}).
#'   Can also be a numeric vector with names (gene IDs); it will be converted
#'   to the expected format.
#' @param n_label Integer. Number of genes to label at the low and high end of
#'   the z-score axis (default 15). Total labels will be at most \code{2 * n_label}.
#' @param z_cutoff Numeric. Absolute z-score threshold for the significance strip
#'   (default 2). Genes with \code{|z| > z_cutoff} are shown as a colored bar
#'   (blue for z < -z_cutoff, red for z > z_cutoff).
#' @param row_title Character. Title for the heatmap row / legend (e.g. "Survival z-score").
#'   If \code{NULL}, the column name of \code{reference} is used when \code{reference}
#'   is a data.frame; otherwise \code{"z-score"}.
#' @param legend_side Character. Where to draw the heatmap legend: \code{"bottom"}
#'   or \code{"right"} (default \code{"bottom"}).
#' @param ... Optional arguments passed to \code{ComplexHeatmap::Heatmap} (e.g.
#'   \code{heatmap_legend_param}).
#'
#' @return The return value of \code{ComplexHeatmap::draw()} (a \code{HeatmapList}
#'   object), invisibly. The plot is drawn in the current device.
#'
#' @examples
#' \dontrun{
#' ref <- derive_reference_from_bulk(
#'   bulk_expression = bulk, phenotype = pheno,
#'   sample_id_column = "sample_id", phenotype_column = "response",
#'   phenotype_type = "binary")
#' plot_reference_signature(ref, row_title = "Survival z-score")
#' }
#'
#' @export
plot_reference_signature <- function(reference,
                                    n_label = 15L,
                                    z_cutoff = 2,
                                    row_title = NULL,
                                    legend_side = c("bottom", "right"),
                                    ...) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "Package 'ComplexHeatmap' is required for plot_reference_signature(). ",
      "Install with: install.packages('ComplexHeatmap')"
    )
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop(
      "Package 'circlize' is required for plot_reference_signature(). ",
      "Install with: install.packages('circlize')"
    )
  }

  legend_side <- match.arg(legend_side)

  # Normalize input to data.frame with genes as rownames and one z column
  if (is.vector(reference) && is.numeric(reference) && !is.matrix(reference)) {
    ref_df <- data.frame(z = reference, row.names = names(reference))
    if (is.null(row_title)) row_title <- "z-score"
  } else if (is.data.frame(reference) || is.matrix(reference)) {
    ref_df <- as.data.frame(reference)
    if (ncol(ref_df) < 1L) stop("'reference' must have at least one column of z-scores")
    if (is.null(rownames(ref_df))) stop("'reference' must have gene names as rownames")
    if (is.null(row_title)) row_title <- colnames(ref_df)[1L]
  } else {
    stop("'reference' must be a data.frame (e.g. from derive_reference_from_bulk) or a named numeric vector")
  }

  z_col <- ref_df[[1L]]
  df_z <- data.frame(
    gene = rownames(ref_df),
    z = as.numeric(z_col),
    stringsAsFactors = FALSE
  )
  df_z <- df_z[is.finite(df_z$z) & !is.na(df_z$gene) & nzchar(df_z$gene), ]
  if (nrow(df_z) == 0L) stop("No finite, non-missing gene z-scores in 'reference'")

  df_z <- df_z[order(df_z$z), ]

  # 1×N matrix, columns ordered by z
  mat <- matrix(df_z$z, nrow = 1L)
  colnames(mat) <- df_z$gene
  rownames(mat) <- row_title

  lab_n <- min(as.integer(n_label), floor(ncol(mat) / 2L))
  mark_at <- c(seq_len(lab_n), (ncol(mat) - lab_n + 1L):ncol(mat))
  mark_labels <- colnames(mat)[mark_at]

  col_fun <- circlize::colorRamp2(
    c(min(df_z$z), 0, max(df_z$z)),
    c("#2166AC", "white", "#B2182B")
  )

  sig_dir <- ifelse(df_z$z <= -z_cutoff, "neg",
    ifelse(df_z$z >= z_cutoff, "pos", "none")
  )

  ha_top <- ComplexHeatmap::HeatmapAnnotation(
    marks = ComplexHeatmap::anno_mark(
      at = mark_at,
      labels = mark_labels,
      side = "top",
      labels_gp = grid::gpar(fontsize = 8),
      link_gp = grid::gpar(col = "grey50", lwd = 0.6),
      padding = grid::unit(2, "mm")
    ),
    which = "column",
    annotation_height = grid::unit(22, "mm"),
    show_annotation_name = FALSE
  )

  sig_label <- sprintf("|z|>%s", z_cutoff)
  anno_sig <- ComplexHeatmap::anno_simple(
    sig_dir,
    col = c(
      "none" = "transparent",
      "neg"  = "#2166AC",
      "pos"  = "#B2182B"
    ),
    height = grid::unit(3, "mm"),
    border = FALSE
  )
  ha_bottom <- do.call(ComplexHeatmap::HeatmapAnnotation, c(
    stats::setNames(list(anno_sig), sig_label),
    list(
      which = "column",
      show_annotation_name = TRUE,
      annotation_name_side = "right"
    )
  ))

  heatmap_legend_param <- list(
    title = row_title,
    direction = "horizontal",
    title_position = "topcenter"
  )
  extra <- list(...)
  if ("heatmap_legend_param" %in% names(extra)) {
    heatmap_legend_param <- utils::modifyList(heatmap_legend_param, extra$heatmap_legend_param)
    extra <- extra[setdiff(names(extra), "heatmap_legend_param")]
  }

  ht <- ComplexHeatmap::Heatmap(
    mat,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha_top,
    bottom_annotation = ha_bottom,
    heatmap_legend_param = heatmap_legend_param,
    ...
  )

  ComplexHeatmap::draw(ht, heatmap_legend_side = legend_side)
}
