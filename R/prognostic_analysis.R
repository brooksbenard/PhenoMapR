#' Define Phenotype Groups (Top and Bottom Percentile)
#'
#' For each score column (dataset), labels cells as adverse (top percentile,
#' highest scores), favorable (bottom percentile, lowest scores), or middle.
#'
#' @param scores Data.frame of phenotype scores from \code{PhenoMap}.
#'   Rows = cells/samples, columns = score variables (e.g. weighted_sum_score_precog_BRCA).
#' @param percentile Fraction of cells in each tail (default 0.05 for top and bottom 5%).
#' @param score_columns Character vector of score column names to use. If NULL,
#'   all numeric columns in \code{scores} are used.
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item \code{cell_id}: same as row names of \code{scores}
#'     \item For each score column: a \code{phenotype_group_<name>} column with
#'       values \code{"Most Adverse"} (top percentile), \code{"Most Favorable"}
#'       (bottom percentile), or \code{"Other"}
#'   }
#'
#' @examples
#' \dontrun{
#' scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
#' groups <- define_phenotype_groups(scores, percentile = 0.05)
#' table(groups$phenotype_group_weighted_sum_score_precog_BRCA)
#' }
#'
#' @export
define_phenotype_groups <- function(scores,
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
      out[[paste0("phenotype_group_", col)]] <- NA_character_
      next
    }

    q_lo <- stats::quantile(v, probs = percentile, na.rm = TRUE)
    q_hi <- stats::quantile(v, probs = 1 - percentile, na.rm = TRUE)

    group <- rep("Other", length(v))
    group[!is.na(v) & v <= q_lo] <- "Most Favorable"
    group[!is.na(v) & v >= q_hi] <- "Most Adverse"
    group[is.na(v)] <- NA_character_

    out[[paste0("phenotype_group_", col)]] <- group
  }

  return(out)
}


#' Find Unique Marker Genes for Adverse and Favorable Phenotype Groups
#'
#' Performs differential expression between the top (adverse) and bottom
#' (favorable) phenotype groups versus the rest. For \code{Seurat} input,
#' uses \code{Seurat::FindMarkers}. For matrix, \code{Matrix}, or
#' \code{SingleCellExperiment} input, uses a Wilcoxon-based path (presto if
#' available, else base R) and does not require Seurat.
#'
#' @param expression Expression matrix (genes x cells), a \code{Matrix}
#'   (e.g. \code{dgCMatrix}), a Seurat object, or a SingleCellExperiment.
#'   Must match cells in \code{group_labels}.
#' @param group_labels Either a character vector of group labels
#'   (\code{"Most Adverse"}, \code{"Most Favorable"}, \code{"Other"}) in the
#'   same order as columns of \code{expression}, or a data.frame from
#'   \code{define_phenotype_groups()} (see \code{group_column}).
#' @param group_column If \code{group_labels} is a data.frame, the name of the
#'   column containing \code{"Most Adverse"} / \code{"Most Favorable"} / \code{"Other"}.
#' @param cell_id_column If \code{group_labels} is a data.frame, the column
#'   name for cell/sample IDs (default \code{"cell_id"}).
#' @param marker_scope Either:
#'   \itemize{
#'     \item \code{"phenotype_groups"}: find markers for the adverse and favorable
#'       phenotype groups globally (cell type agnostic; default).
#'     \item \code{"cell_type_specific"}: find markers separately within each
#'       cell type, comparing the adverse phenotype group vs the rest and the
#'       favorable phenotype group vs the rest, *within the same cell type*.
#'   }
#' @param cell_type_column When \code{marker_scope = "cell_type_specific"}, the
#'   column in \code{group_labels} that contains cell type labels.
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
#' @param max_cells_per_ident When any phenotype group exceeds this many cells,
#'   subsample to this limit before FindMarkers (default 5000). Reduces memory
#'   for large objects. Set to \code{Inf} to disable.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers} when input is a Seurat object (ignored for matrix/SCE/Matrix input).
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
#'   When \code{marker_scope = "cell_type_specific"}, the returned data.frames
#'   additionally include a \code{cell_type} column with the cell type for each
#'   marker result.
#'
#' @examples
#' \dontrun{
#' scores <- PhenoMap(seurat_obj, reference = "precog", cancer_type = "BRCA")
#' groups <- define_phenotype_groups(scores, percentile = 0.05)
#' markers <- find_phenotype_markers(
#'   seurat_obj,
#'   group_labels = groups,
#'   group_column = "phenotype_group_weighted_sum_score_precog_BRCA",
#'   cell_id_column = "cell_id"
#' )
#' head(markers$adverse_markers)
#' head(markers$favorable_markers)
#' }
#'
#' @export
find_phenotype_markers <- function(expression,
                                    group_labels,
                                    group_column = NULL,
                                    cell_id_column = "cell_id",
                                    marker_scope = c("phenotype_groups", "cell_type_specific"),
                                    cell_type_column = NULL,
                                    assay = NULL,
                                    slot = "data",
                                    test.use = c("wilcox", "t", "roc", "negbinom", "poisson", "LR"),
                                    min.pct = 0.1,
                                    logfc.threshold = 0.25,
                                    pval_threshold = 0.05,
                                    verbose = TRUE,
                                    max_cells_per_ident = 5000L,
                                    ...) {
  test.use <- match.arg(test.use)
  marker_scope <- match.arg(marker_scope)

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

  # Optional: resolve cell type labels for cell-type-specific marker detection
  cell_type_vec <- NULL
  if (marker_scope == "cell_type_specific") {
    if (!is.data.frame(group_labels)) {
      stop(
        "marker_scope = 'cell_type_specific' requires 'group_labels' to be a data.frame ",
        "with at least the cell id column and a cell type column."
      )
    }
    cols <- infer_group_labels_columns(
      group_labels = group_labels,
      cell_id_column = cell_id_column,
      group_column = group_column
    )
    cell_type_column_use <- cell_type_column %||%
      infer_cell_type_column(
        group_labels = group_labels,
        id_col = cols$id_col,
        group_column = cols$group_column
      )
    cell_type_vec <- resolve_labels_from_df(
      group_labels = group_labels,
      cell_ids = cell_ids,
      id_col = cols$id_col,
      label_column = cell_type_column_use
    )
  }

  n_adverse <- sum(group_vec == "Most Adverse", na.rm = TRUE)
  n_favorable <- sum(group_vec == "Most Favorable", na.rm = TRUE)
  if (n_adverse < 5L) {
    warning("Fewer than 5 Most Adverse cells; marker results may be unreliable")
  }
  if (n_favorable < 5L) {
    warning("Fewer than 5 Most Favorable cells; marker results may be unreliable")
  }

  # Non-Seurat path: matrix/Matrix/SCE — use internal Wilcoxon (or presto) without requiring Seurat.
  # Also forced for cell-type-specific mode to avoid repeated Seurat metadata/idents reshaping.
  use_matrix_path <- marker_scope == "cell_type_specific" ||
    !inherits(expression, "Seurat") || !requireNamespace("Seurat", quietly = TRUE)
  if (use_matrix_path) {
    if (verbose) {
      message(glue::glue(
        "Using matrix-based marker detection (no Seurat): Most Adverse n={n_adverse}, Most Favorable n={n_favorable}"
      ))
    }
    if (marker_scope == "cell_type_specific") {
      out <- run_markers_on_matrix_by_celltype(
        mat = expr_info$matrix,
        group_vec = group_vec,
        cell_type_vec = cell_type_vec,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        pval_threshold = pval_threshold,
        max_cells_per_ident = max_cells_per_ident,
        verbose = verbose
      )
    } else {
      out <- run_markers_on_matrix(
        mat = expr_info$matrix,
        group_vec = group_vec,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        pval_threshold = pval_threshold,
        max_cells_per_ident = max_cells_per_ident,
        verbose = verbose
      )
    }
    return(out)
  }

  # nocov start - Seurat FindMarkers path (optional; tested in test-find-phenotype-markers when Seurat installed)
  # Seurat path
  if (verbose) {
    message(glue::glue(
      "Using Seurat FindMarkers: Most Adverse n={n_adverse}, Most Favorable n={n_favorable}"
    ))
  }
  seurat_obj <- get_or_create_seurat_for_markers(expression, expr_info, group_vec, assay, slot)
  meta_col <- "PhenoMapR_phenotype_group"
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

  if (inherits(expression, "Seurat") && meta_col %in% names(seurat_obj[[]])) {
    seurat_obj[[meta_col]] <- NULL
  }
  # nocov end

  adverse_markers <- standardize_findmarkers_output(adverse_markers, pval_threshold)
  favorable_markers <- standardize_findmarkers_output(favorable_markers, pval_threshold)

  list(
    adverse_markers = adverse_markers,
    favorable_markers = favorable_markers
  )
}


#' Infer helper columns from a group_labels data.frame
#'
#' @keywords internal
infer_group_labels_columns <- function(group_labels, cell_id_column, group_column) {
  id_col <- cell_id_column
  if (!id_col %in% names(group_labels)) id_col <- names(group_labels)[1]

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

  list(id_col = id_col, group_column = group_column)
}


#' Infer which column contains cell type labels
#'
#' @keywords internal
infer_cell_type_column <- function(group_labels, id_col, group_column) {
  candidates <- setdiff(names(group_labels), c(id_col, group_column))
  if (length(candidates) == 0) {
    stop("marker_scope='cell_type_specific' requires a cell type column in 'group_labels' (other than the cell id and group columns).")
  }

  # Prefer non-group columns with multiple unique values
  group_vals <- c("Most Adverse", "Most Favorable", "Other")
  scored <- sapply(candidates, function(cc) {
    x <- group_labels[[cc]]
    if (!is.character(x) && !is.factor(x)) return(0)
    x_u <- unique(x)
    x_u <- x_u[!is.na(x_u)]
    # Reject if all values are just the phenotype group strings
    if (length(x_u) == 0) return(0)
    if (all(x_u %in% group_vals)) return(0)
    length(x_u)
  })
  if (all(scored == 0)) {
    stop("marker_scope='cell_type_specific' could not infer a cell type column in 'group_labels'. Please pass 'cell_type_column'.")
  }
  candidates[which.max(scored)]
}


#' Resolve a label column from group_labels data.frame to match cell_ids
#'
#' @keywords internal
resolve_labels_from_df <- function(group_labels, cell_ids, id_col, label_column) {
  lab_df <- group_labels[, c(id_col, label_column), drop = FALSE]
  lab_df <- lab_df[lab_df[[id_col]] %in% cell_ids, , drop = FALSE]
  lab_df <- lab_df[match(cell_ids, lab_df[[id_col]]), , drop = FALSE]
  vec <- lab_df[[label_column]]
  vec[is.na(lab_df[[id_col]])] <- NA
  vec
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
  # nocov start - Seurat FindMarkers column names (matrix path uses different names)
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
  # nocov end
  df <- df[df$p_val <= pval_threshold, , drop = FALSE]
  df <- df[order(df$p_val), , drop = FALSE]
  df
}


#' Run marker detection on a matrix (no Seurat). Uses presto::wilcoxauc if available, else base Wilcoxon.
#' @keywords internal
run_markers_on_matrix <- function(mat,
                                 group_vec,
                                 min.pct = 0.1,
                                 logfc.threshold = 0.25,
                                 pval_threshold = 0.05,
                                 max_cells_per_ident = 5000L,
                                 verbose = TRUE) {
  cell_ids <- colnames(mat)
  if (is.null(cell_ids)) cell_ids <- paste0("Cell_", seq_len(ncol(mat)))
  if (length(group_vec) != ncol(mat)) {
    stop("Length of group_vec must match number of columns in expression matrix")
  }
  # Subsample per group
  # nocov start - subsampling only when group size > max_cells_per_ident
  if (is.finite(max_cells_per_ident)) {
    idx_keep <- integer(0)
    for (grp in c("Most Adverse", "Most Favorable", "Other")) {
      idx <- which(group_vec == grp)
      if (length(idx) > max_cells_per_ident) {
        idx <- idx[sample.int(length(idx), max_cells_per_ident)]
        if (verbose) {
          message(glue::glue("Subsampled {grp} from {sum(group_vec == grp)} to {max_cells_per_ident} cells (memory limit)"))
        }
      }
      idx_keep <- c(idx_keep, idx)
    }
    mat <- mat[, idx_keep, drop = FALSE]
    group_vec <- group_vec[idx_keep]
  }
  # nocov end
  mat_dense <- as.matrix(mat)
  genes <- rownames(mat_dense)
  if (is.null(genes)) genes <- paste0("Gene_", seq_len(nrow(mat_dense)))

  # nocov start - presto path (optional dependency)
  # Use presto if available (fast Wilcoxon). X = genes x cells, y = group labels
  if (requireNamespace("presto", quietly = TRUE)) {
    auc_out <- tryCatch(
      presto::wilcoxauc(mat_dense, group_vec),
      error = function(e) NULL
    )
    if (!is.null(auc_out) && nrow(auc_out) > 0) {
      format_one_group <- function(ident) {
        sub <- auc_out[auc_out$group == ident, , drop = FALSE]
        if (nrow(sub) == 0) {
          return(data.frame(
            gene = character(0), avg_log2FC = numeric(0),
            pct_in_group = numeric(0), pct_rest = numeric(0),
            p_val = numeric(0), p_adj = numeric(0),
            stringsAsFactors = FALSE
          ))
        }
        sub$gene <- sub$feature
        sub$avg_log2FC <- log2(exp(1)) * sub$logFC
        sub$pct_in_group <- sub$pct_in
        sub$pct_rest <- sub$pct_out
        sub$p_val <- sub$pval
        sub$p_adj <- sub$pval
        if ("padj" %in% names(sub)) sub$p_adj <- sub$padj
        sub <- sub[sub$p_val <= pval_threshold & abs(sub$avg_log2FC) >= logfc.threshold, , drop = FALSE]
        sub <- sub[order(sub$p_val), , drop = FALSE]
        sub[, c("gene", "avg_log2FC", "pct_in_group", "pct_rest", "p_val", "p_adj")]
      }
      adverse_markers <- format_one_group("Most Adverse")
      favorable_markers <- format_one_group("Most Favorable")
      return(list(adverse_markers = adverse_markers, favorable_markers = favorable_markers))
    }
  }
  # nocov end

  # Fallback: base R Wilcoxon per gene
  adverse_idx <- group_vec == "Most Adverse"
  favorable_idx <- group_vec == "Most Favorable"
  rest_adverse <- !adverse_idx & !is.na(group_vec)
  rest_favorable <- !favorable_idx & !is.na(group_vec)
  n_adv <- sum(adverse_idx)
  n_fav <- sum(favorable_idx)
  if (n_adv < 2 || n_fav < 2) {
    # nocov start - too few cells in either group
    return(list(
      adverse_markers = standardize_findmarkers_output(NULL, pval_threshold),
      favorable_markers = standardize_findmarkers_output(NULL, pval_threshold)
    ))
    # nocov end
  }
  run_wilcox_one_vs_rest <- function(ident_label, is_in_group) {
    in_grp <- mat_dense[, is_in_group, drop = FALSE]
    out_grp <- mat_dense[, !is_in_group & !is.na(group_vec), drop = FALSE]
    n_in <- ncol(in_grp)
    n_out <- ncol(out_grp)
    pct_in <- rowMeans(in_grp > 0, na.rm = TRUE)
    pct_out <- rowMeans(out_grp > 0, na.rm = TRUE)
    keep <- (pct_in >= min.pct | pct_out >= min.pct)
    if (sum(keep) == 0) {
      # nocov start - no genes pass min.pct
      return(data.frame(
        gene = character(0), avg_log2FC = numeric(0),
        pct_in_group = numeric(0), pct_rest = numeric(0),
        p_val = numeric(0), p_adj = numeric(0),
        stringsAsFactors = FALSE
      ))
      # nocov end
    }
    genes_keep <- genes[keep]
    pct_in <- pct_in[keep]
    pct_out <- pct_out[keep]
    in_grp <- in_grp[keep, , drop = FALSE]
    out_grp <- out_grp[keep, , drop = FALSE]
    p_vals <- numeric(length(genes_keep))
    log2fc <- numeric(length(genes_keep))
    for (i in seq_along(genes_keep)) {
      x_in <- as.numeric(in_grp[i, ])
      x_out <- as.numeric(out_grp[i, ])
      wt <- tryCatch(stats::wilcox.test(x_in, x_out, alternative = "two.sided", exact = FALSE), error = function(e) NULL)
      p_vals[i] <- if (!is.null(wt)) wt$p.value else NA_real_
      m_in <- mean(x_in, na.rm = TRUE)
      m_out <- mean(x_out, na.rm = TRUE)
      log2fc[i] <- if (m_out > 0) log2((m_in + 1e-10) / (m_out + 1e-10)) else NA_real_
    }
    df <- data.frame(
      gene = genes_keep,
      avg_log2FC = log2fc,
      pct_in_group = pct_in,
      pct_rest = pct_out,
      p_val = p_vals,
      p_adj = p_vals,
      stringsAsFactors = FALSE
    )
    df <- df[df$p_val <= pval_threshold & abs(df$avg_log2FC) >= logfc.threshold, , drop = FALSE]
    df <- df[order(df$p_val), , drop = FALSE]
    df
  }
  adverse_markers <- run_wilcox_one_vs_rest("Most Adverse", adverse_idx)
  favorable_markers <- run_wilcox_one_vs_rest("Most Favorable", favorable_idx)
  list(adverse_markers = adverse_markers, favorable_markers = favorable_markers)
}


#' Run marker detection on a matrix by cell type
#'
#' For each cell type, compute markers for:
#'   - Most Adverse within that cell type vs all other groups within that cell type
#'   - Most Favorable within that cell type vs all other groups within that cell type
#'
#' @keywords internal
run_markers_on_matrix_by_celltype <- function(mat,
                                               group_vec,
                                               cell_type_vec,
                                               min.pct = 0.1,
                                               logfc.threshold = 0.25,
                                               pval_threshold = 0.05,
                                               max_cells_per_ident = 5000L,
                                               verbose = TRUE) {
  cell_ids <- colnames(mat)
  if (is.null(cell_ids)) cell_ids <- paste0("Cell_", seq_len(ncol(mat)))
  if (length(group_vec) != ncol(mat)) {
    stop("Length of 'group_vec' must match number of columns in expression matrix")
  }
  if (length(cell_type_vec) != ncol(mat)) {
    stop("Length of 'cell_type_vec' must match number of columns in expression matrix")
  }

  cell_type_vec <- as.character(cell_type_vec)
  cell_types <- sort(unique(cell_type_vec[!is.na(cell_type_vec)]))
  if (length(cell_types) == 0) {
    return(list(
      adverse_markers = standardize_findmarkers_output(NULL, pval_threshold),
      favorable_markers = standardize_findmarkers_output(NULL, pval_threshold)
    ))
  }

  adverse_all <- list()
  favorable_all <- list()

  for (ct in cell_types) {
    idx_ct <- which(cell_type_vec == ct)
    if (length(idx_ct) < 4) next

    res_ct <- run_markers_on_matrix(
      mat = mat[, idx_ct, drop = FALSE],
      group_vec = group_vec[idx_ct],
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      pval_threshold = pval_threshold,
      max_cells_per_ident = max_cells_per_ident,
      verbose = verbose
    )

    adverse_ct <- res_ct$adverse_markers
    favorable_ct <- res_ct$favorable_markers
    if (nrow(adverse_ct) > 0) {
      adverse_ct$cell_type <- rep(ct, nrow(adverse_ct))
    } else {
      adverse_ct$cell_type <- character(0)
    }
    if (nrow(favorable_ct) > 0) {
      favorable_ct$cell_type <- rep(ct, nrow(favorable_ct))
    } else {
      favorable_ct$cell_type <- character(0)
    }
    adverse_all[[ct]] <- adverse_ct
    favorable_all[[ct]] <- favorable_ct
  }

  adverse_markers <- if (length(adverse_all) > 0) {
    do.call(rbind, adverse_all)
  } else {
    standardize_findmarkers_output(NULL, pval_threshold)
  }
  favorable_markers <- if (length(favorable_all) > 0) {
    do.call(rbind, favorable_all)
  } else {
    standardize_findmarkers_output(NULL, pval_threshold)
  }

  # Ensure consistent column order (cell_type at end)
  if (!"cell_type" %in% names(adverse_markers)) adverse_markers$cell_type <- NA_character_
  if (!"cell_type" %in% names(favorable_markers)) favorable_markers$cell_type <- NA_character_
  adverse_markers <- adverse_markers[, c(setdiff(names(adverse_markers), "cell_type"), "cell_type"), drop = FALSE]
  favorable_markers <- favorable_markers[, c(setdiff(names(favorable_markers), "cell_type"), "cell_type"), drop = FALSE]

  list(adverse_markers = adverse_markers, favorable_markers = favorable_markers)
}


#' Get existing Seurat object or create one from matrix/SCE; add phenotype group to metadata
#' @keywords internal
# nocov start - only used by Seurat path in find_phenotype_markers (nocov'd)
get_or_create_seurat_for_markers <- function(expression, expr_info, group_vec, assay, slot) {
  if (inherits(expression, "Seurat")) {
    assay <- assay %||% "RNA"
    obj <- expression
    # Ensure cells in same order as expr_info
    cell_ids <- expr_info$cell_names
    obj <- obj[, cell_ids]
    obj[["PhenoMapR_phenotype_group"]] <- group_vec
    attr(obj, "PhenoMapR_assay") <- assay
    return(obj)
  }
  # Matrix or SCE: create temporary Seurat object
  mat <- expr_info$matrix
  if (inherits(mat, "sparseMatrix")) {
    obj <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA")
  } else {
    obj <- Seurat::CreateSeuratObject(counts = as.matrix(mat), assay = "RNA")
  }
  # For matrix input, CreateSeuratObject only has "counts". Use "data" layer;
  # avoid double-normalization when input is already normalized (e.g. TISCH2).
  if (slot == "data") {
    if (is_likely_normalized(mat)) {
      Seurat::SetAssayData(obj, layer = "data", new.data = mat, assay = "RNA")
    } else {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    }
  }
  obj[["PhenoMapR_phenotype_group"]] <- group_vec
  attr(obj, "PhenoMapR_assay") <- "RNA"
  obj
}
# nocov end


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
  if (inherits(expression, "Matrix")) {
    return(list(
      matrix = expression,
      cell_names = colnames(expression),
      gene_names = rownames(expression)
    ))
  }
  # nocov start - Seurat/SCE paths only used when find_phenotype_markers gets Seurat/SCE (Seurat path nocov'd)
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
  # nocov end
  stop("'expression' must be a matrix, Matrix, data.frame, Seurat, or SingleCellExperiment")
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
  # nocov start - data.frame resolution (tested via find_phenotype_markers with group_labels df)
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
  # nocov end
}


#' Default value when NULL
#' @name or-null
#' @aliases %||%
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
