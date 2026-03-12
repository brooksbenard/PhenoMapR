#' Define Prognostic Groups (Top and Bottom Percentile)
#'
#' For each score column (dataset), labels cells as adverse (top percentile,
#' highest scores), favorable (bottom percentile, lowest scores), or middle.
#'
#' @param scores Data.frame of prognostic scores from \code{PhenoMap}.
#'   Rows = cells/samples, columns = score variables (e.g. weighted_sum_score_precog_BRCA).
#' @param percentile Fraction of cells in each tail (default 0.05 for top and bottom 5%).
#' @param score_columns Character vector of score column names to use. If NULL,
#'   all numeric columns in \code{scores} are used.
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item \code{cell_id}: same as row names of \code{scores}
#'     \item For each score column: a \code{prognostic_group_<name>} column with
#'       values \code{"Most Adverse"} (top percentile), \code{"Most Favorable"}
#'       (bottom percentile), or \code{"Other"}
#'   }
#'
#' @examples
#' \dontrun{
#' scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
#' groups <- define_prognostic_groups(scores, percentile = 0.05)
#' table(groups$prognostic_group_weighted_sum_score_precog_BRCA)
#' }
#'
#' @export
define_prognostic_groups <- function(scores,
                                     percentile = 0.05,
                                     score_columns = NULL) {
  if (!is.data.frame(scores)) {
    stop("'scores' must be a data.frame from PhenoMap()")
  }
  if (percentile <= 0 || percentile >= 0.5) {
    stop("'percentile' must be between 0 and 0.5 (e.g. 0.05 for 5%)")
  }

  cell_ids <- rownames(scores)
  if (is.null(cell_ids)) {
    cell_ids <- as.character(seq_len(nrow(scores)))
  }

  if (is.null(score_columns)) {
    score_columns <- names(scores)[vapply(scores, is.numeric, logical(1))]
  }
  if (length(score_columns) == 0) {
    stop("No numeric score columns found in 'scores'")
  }

  out <- data.frame(cell_id = cell_ids, stringsAsFactors = FALSE)
  rownames(out) <- cell_ids

  for (col in score_columns) {
    if (!col %in% names(scores)) {
      warning(glue::glue("Score column '{col}' not found, skipping"))
      next
    }
    v <- scores[[col]]
    if (!is.numeric(v)) next

    n <- sum(!is.na(v))
    if (n == 0) {
      out[[paste0("prognostic_group_", col)]] <- NA_character_
      next
    }

    q_lo <- stats::quantile(v, probs = percentile, na.rm = TRUE)
    q_hi <- stats::quantile(v, probs = 1 - percentile, na.rm = TRUE)

    group <- rep("Other", length(v))
    group[!is.na(v) & v <= q_lo] <- "Most Favorable"
    group[!is.na(v) & v >= q_hi] <- "Most Adverse"
    group[is.na(v)] <- NA_character_

    out[[paste0("prognostic_group_", col)]] <- group
  }

  return(out)
}


#' Find Unique Marker Genes for Adverse and Favorable Prognostic Groups
#'
#' Uses Seurat's \code{FindMarkers} to perform differential expression between
#' the top (adverse) and bottom (favorable) prognostic groups versus the rest.
#' Requires the Seurat package.
#'
#' @param expression Expression matrix (genes x cells), a Seurat object, or a
#'   SingleCellExperiment. Must match cells in \code{group_labels}.
#' @param group_labels Either a character vector of group labels
#'   (\code{"Most Adverse"}, \code{"Most Favorable"}, \code{"Other"}) in the
#'   same order as columns of \code{expression}, or a data.frame from
#'   \code{define_prognostic_groups()} (see \code{group_column}).
#' @param group_column If \code{group_labels} is a data.frame, the name of the
#'   column containing \code{"Most Adverse"} / \code{"Most Favorable"} / \code{"Other"}.
#' @param cell_id_column If \code{group_labels} is a data.frame, the column
#'   name for cell/sample IDs (default \code{"cell_id"}).
#' @param assay Assay name for Seurat/SCE (e.g. \code{"RNA"}).
#' @param slot Layer for Seurat: \code{"data"}, \code{"counts"}, or
#'   \code{"scale.data"} (default \code{"data"}).
#' @param test.use Seurat \code{test.use}: \code{"wilcox"} (default), \code{"t"},
#'   \code{"roc"}, \code{"negbinom"}, \code{"poisson"}, or \code{"LR"}.
#' @param min.pct Minimum fraction of cells expressing the gene in either group
#'   (default 0.1). Passed to \code{FindMarkers}.
#' @param logfc.threshold Minimum absolute log2 fold change (default 0.25).
#'   Passed to \code{FindMarkers}.
#' @param pval_threshold Maximum unadjusted p-value to include (default 0.05).
#' @param verbose Print progress messages (default TRUE).
#' @param max_cells_per_ident When any prognostic group exceeds this many cells,
#'   subsample to this limit before FindMarkers (default 5000). Reduces memory
#'   for large objects. Set to \code{Inf} to disable.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{adverse_markers}: data.frame of genes that are markers of
#'       the adverse (top score) group (vs rest).
#'     \item \code{favorable_markers}: data.frame of genes that are markers of
#'       the favorable (bottom score) group (vs rest).
#'   }
#'   Each data.frame has columns: \code{gene}, \code{avg_log2FC},
#'   \code{pct_in_group}, \code{pct_rest}, \code{p_val}, \code{p_adj}.
#'
#' @examples
#' \dontrun{
#' scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
#' groups <- define_prognostic_groups(scores, percentile = 0.05)
#' markers <- find_prognostic_markers(
#'   seurat_obj,
#'   group_labels = groups,
#'   group_column = "prognostic_group_weighted_sum_score_precog_BRCA",
#'   cell_id_column = "cell_id"
#' )
#' head(markers$adverse_markers)
#' head(markers$favorable_markers)
#' }
#'
#' @export
find_prognostic_markers <- function(expression,
                                    group_labels,
                                    group_column = NULL,
                                    cell_id_column = "cell_id",
                                    assay = NULL,
                                    slot = "data",
                                    test.use = c("wilcox", "t", "roc", "negbinom", "poisson", "LR"),
                                    min.pct = 0.1,
                                    logfc.threshold = 0.25,
                                    pval_threshold = 0.05,
                                    verbose = TRUE,
                                    max_cells_per_ident = 5000L,
                                    ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("find_prognostic_markers() requires the Seurat package. Install with: install.packages('Seurat')")
  }
  test.use <- match.arg(test.use)

  # Resolve expression and cell order
  expr_info <- process_expression_for_markers(
    expression,
    assay = assay,
    slot = slot
  )
  cell_ids <- expr_info$cell_names
  slot_use <- expr_info$slot_used %||% slot

  group_vec <- resolve_group_labels(
    group_labels,
    cell_ids,
    group_column = group_column,
    cell_id_column = cell_id_column
  )

  n_adverse <- sum(group_vec == "Most Adverse", na.rm = TRUE)
  n_favorable <- sum(group_vec == "Most Favorable", na.rm = TRUE)
  if (n_adverse < 5L) {
    warning("Fewer than 5 Most Adverse cells; marker results may be unreliable")
  }
  if (n_favorable < 5L) {
    warning("Fewer than 5 Most Favorable cells; marker results may be unreliable")
  }
  if (verbose) {
    message(glue::glue(
      "Using Seurat FindMarkers: Most Adverse n={n_adverse}, Most Favorable n={n_favorable}"
    ))
  }

  # Get or create Seurat object and add prognostic group
  seurat_obj <- get_or_create_seurat_for_markers(expression, expr_info, group_vec, assay, slot)
  meta_col <- "PhenoMapR_prognostic_group"
  Seurat::Idents(seurat_obj) <- seurat_obj[[meta_col]][, 1]
  assay_use <- assay %||% attr(seurat_obj, "PhenoMapR_assay") %||% "RNA"
  if (slot_use == "counts" && inherits(seurat_obj, "Seurat")) {
    dat_check <- tryCatch(
      Seurat::GetAssayData(seurat_obj, assay = assay_use, layer = "data"),
      error = function(e) NULL
    )
    if (is.null(dat_check) || (nrow(dat_check) == 0 || ncol(dat_check) == 0)) {
      if (verbose) message("Normalizing assay '", assay_use, "' for FindMarkers (data layer was empty)")
      seurat_obj <- Seurat::NormalizeData(seurat_obj, assay = assay_use, verbose = verbose)
      slot_use <- "data"
    }
  }

  # Subsample when any ident exceeds max_cells_per_ident (reduces memory for large objects)
  if (is.finite(max_cells_per_ident)) {
    idents <- as.character(Seurat::Idents(seurat_obj))
    cells_keep <- character(0)
    for (grp in c("Most Adverse", "Most Favorable", "Other")) {
      idx <- which(idents == grp)
      if (length(idx) > max_cells_per_ident) {
        idx <- idx[sample.int(length(idx), max_cells_per_ident)]
        if (verbose) {
          message(glue::glue("Subsampled {grp} from {sum(idents == grp)} to {max_cells_per_ident} cells (memory limit)"))
        }
      }
      cells_keep <- c(cells_keep, colnames(seurat_obj)[idx])
    }
    seurat_obj <- tryCatch(
      seurat_obj[, cells_keep],
      error = function(e) {
        if (verbose) {
          message(glue::glue("Subsetting failed ({e$message}); using full object (may use more memory)"))
        }
        seurat_obj
      }
    )
  }

  # Adverse vs rest
  gc(verbose = FALSE)
  assay_use <- assay %||% attr(seurat_obj, "PhenoMapR_assay") %||% "RNA"
  adverse_markers <- tryCatch(
    Seurat::FindMarkers(
      seurat_obj,
      ident.1 = "Most Adverse",
      ident.2 = NULL,
      assay = assay_use,
      slot = slot_use,
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      verbose = verbose,
      ...
    ),
    error = function(e) {
      warning(glue::glue("FindMarkers (adverse) failed: {e$message}"))
      data.frame(
        p_val = numeric(0), avg_log2FC = numeric(0),
        pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)
      )
    }
  )
  gc(verbose = FALSE)
  # Favorable vs rest
  favorable_markers <- tryCatch(
    Seurat::FindMarkers(
      seurat_obj,
      ident.1 = "Most Favorable",
      ident.2 = NULL,
      assay = assay_use,
      slot = slot_use,
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      verbose = verbose,
      ...
    ),
    error = function(e) {
      warning(glue::glue("FindMarkers (Most Favorable) failed: {e$message}"))
      data.frame(
        p_val = numeric(0), avg_log2FC = numeric(0),
        pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0)
      )
    }
  )

  # Remove temporary metadata if we added it to the user's object
  if (inherits(expression, "Seurat") && meta_col %in% names(seurat_obj[[]])) {
    seurat_obj[[meta_col]] <- NULL
  }

  # Standardize output columns to match documented format
  adverse_markers <- standardize_findmarkers_output(adverse_markers, pval_threshold)
  favorable_markers <- standardize_findmarkers_output(favorable_markers, pval_threshold)

  list(
    adverse_markers = adverse_markers,
    favorable_markers = favorable_markers
  )
}


#' Convert FindMarkers result to standard columns and filter by p-value
#' @keywords internal
standardize_findmarkers_output <- function(df, pval_threshold) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      gene = character(0),
      avg_log2FC = numeric(0),
      pct_in_group = numeric(0),
      pct_rest = numeric(0),
      p_val = numeric(0),
      p_adj = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  df$gene <- rownames(df)
  rownames(df) <- NULL
  # Seurat may return avg_log2FC or avg_logFC
  if (!"avg_log2FC" %in% names(df) && "avg_logFC" %in% names(df)) {
    df$avg_log2FC <- df$avg_logFC
    df$avg_logFC <- NULL
  }
  names(df)[names(df) == "pct.1"] <- "pct_in_group"
  names(df)[names(df) == "pct.2"] <- "pct_rest"
  if ("p_val_adj" %in% names(df)) {
    df$p_adj <- df$p_val_adj
    df$p_val_adj <- NULL
  }
  if (!"p_adj" %in% names(df)) {
    df$p_adj <- df$p_val
  }
  df <- df[df$p_val <= pval_threshold, , drop = FALSE]
  df <- df[order(df$p_val), , drop = FALSE]
  df
}


#' Get existing Seurat object or create one from matrix/SCE; add prognostic group to metadata
#' @keywords internal
get_or_create_seurat_for_markers <- function(expression, expr_info, group_vec, assay, slot) {
  if (inherits(expression, "Seurat")) {
    assay <- assay %||% "RNA"
    obj <- expression
    # Ensure cells in same order as expr_info
    cell_ids <- expr_info$cell_names
    obj <- obj[, cell_ids]
    obj[["PhenoMapR_prognostic_group"]] <- group_vec
    attr(obj, "PhenoMapR_assay") <- assay
    return(obj)
  }
  # Matrix or SCE: create temporary Seurat object
  mat <- expr_info$matrix
  if (inherits(mat, "sparseMatrix")) {
    obj <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA", verbose = FALSE)
  } else {
    obj <- Seurat::CreateSeuratObject(counts = as.matrix(mat), assay = "RNA", verbose = FALSE)
  }
  # For matrix input, CreateSeuratObject only has "counts"; normalize if slot is "data"
  if (slot == "data") {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  obj[["PhenoMapR_prognostic_group"]] <- group_vec
  attr(obj, "PhenoMapR_assay") <- "RNA"
  obj
}


#' Process expression input for marker finding
#' @keywords internal
process_expression_for_markers <- function(expression, assay = NULL, slot = "data") {
  if (is.matrix(expression) || is.data.frame(expression)) {
    if (is.data.frame(expression)) {
      expression <- as.matrix(expression)
    }
    return(list(
      matrix = expression,
      cell_names = colnames(expression),
      gene_names = rownames(expression)
    ))
  }
  if (inherits(expression, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required for Seurat input")
    }
    assay <- assay %||% "RNA"
    if (!assay %in% names(expression@assays)) {
      stop(glue::glue("Assay '{assay}' not found. Available: {paste(names(expression@assays), collapse = ', ')}"))
    }
    layer_map <- c(data = "data", counts = "counts", scale.data = "scale.data")
    layer_name <- if (slot %in% names(layer_map)) layer_map[slot] else "data"
    mat <- tryCatch(
      Seurat::GetAssayData(expression, assay = assay, layer = layer_name),
      error = function(e) NULL
    )
    slot_used <- slot
    if (is.null(mat)) {
      mat <- tryCatch(
        Seurat::GetAssayData(expression, assay = assay, slot = slot),
        error = function(e) NULL
      )
    }
    if (is.null(mat) || (nrow(mat) == 0 || ncol(mat) == 0)) {
      mat <- tryCatch(
        Seurat::GetAssayData(expression, assay = assay, layer = "counts"),
        error = function(e) Seurat::GetAssayData(expression, assay = assay, slot = "counts")
      )
      if (!is.null(mat) && (nrow(mat) > 0 && ncol(mat) > 0)) slot_used <- "counts"
    }
    mat <- as.matrix(mat)
    return(list(
      matrix = mat,
      cell_names = colnames(mat),
      gene_names = rownames(mat),
      slot_used = slot_used
    ))
  }
  if (inherits(expression, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment required for SCE input")
    }
    assay <- assay %||% "logcounts"
    anames <- SummarizedExperiment::assayNames(expression)
    if (!assay %in% anames) assay <- anames[1]
    mat <- as.matrix(SummarizedExperiment::assay(expression, assay))
    return(list(
      matrix = mat,
      cell_names = colnames(mat),
      gene_names = rownames(mat)
    ))
  }
  stop("'expression' must be a matrix, data.frame, Seurat, or SingleCellExperiment")
}


#' Resolve group labels to a vector aligned with cell_ids
#' @keywords internal
resolve_group_labels <- function(group_labels,
                                 cell_ids,
                                 group_column = NULL,
                                 cell_id_column = "cell_id") {
  if (is.vector(group_labels) && !is.data.frame(group_labels)) {
    if (length(group_labels) != length(cell_ids)) {
      stop("Length of 'group_labels' must match number of cells in expression")
    }
    return(group_labels)
  }
  if (!is.data.frame(group_labels)) {
    stop("'group_labels' must be a character vector or a data.frame with group and cell ID columns")
  }

  id_col <- cell_id_column
  if (!id_col %in% names(group_labels)) {
    id_col <- names(group_labels)[1]
  }
  if (is.null(group_column)) {
    group_candidates <- setdiff(names(group_labels), id_col)
    group_candidates <- group_candidates[
      vapply(group_labels[group_candidates], function(x) {
        is.character(x) && all(unique(x) %in% c("Most Adverse", "Most Favorable", "Other", NA_character_), na.rm = TRUE)
      }, logical(1))
    ]
    if (length(group_candidates) == 0) {
      stop("No column in 'group_labels' contains only 'Most Adverse'/'Most Favorable'/'Other'. Specify 'group_column'.")
    }
    if (length(group_candidates) > 1) {
      stop("Multiple group columns found. Specify 'group_column'.")
    }
    group_column <- group_candidates
  }
  if (!group_column %in% names(group_labels)) {
    stop(glue::glue("'group_column' '{group_column}' not found in group_labels"))
  }

  lab_df <- group_labels[, c(id_col, group_column), drop = FALSE]
  lab_df <- lab_df[lab_df[[id_col]] %in% cell_ids, , drop = FALSE]
  lab_df <- lab_df[match(cell_ids, lab_df[[id_col]]), , drop = FALSE]
  vec <- lab_df[[group_column]]
  vec[is.na(lab_df[[id_col]])] <- NA_character_
  vec
}


#' Default value when NULL
#' @name or-null
#' @aliases %||%
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
