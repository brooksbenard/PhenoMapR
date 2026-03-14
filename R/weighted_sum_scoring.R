#' Calculate Weighted Sum Scores
#'
#' Core function to calculate weighted sum of expression and z-scores
#'
#' @keywords internal
calculate_weighted_scores <- function(expression_matrix,
                                      reference_data,
                                      z_score_cutoff = 2,
                                      pseudobulk = FALSE,
                                      score_name = "weighted_sum_score",
                                      verbose = TRUE) {
  
  # Validate inputs
  validate_expression_matrix(expression_matrix)
  
  if (!is.data.frame(reference_data) && !is.matrix(reference_data)) {
    stop("reference_data must be a data.frame or matrix")
  }
  
  reference_data <- as.data.frame(reference_data)
  
  # Cell IDs: use colnames if present; otherwise generate and assign (validate_expression_matrix
  # may not modify the caller's matrix due to R copy-on-write)
  cell_ids <- colnames(expression_matrix)
  if (is.null(cell_ids) || length(cell_ids) == 0) {
    cell_ids <- paste0("Cell_", seq_len(ncol(expression_matrix)))
    colnames(expression_matrix) <- cell_ids
  }
  
  # Initialize scores dataframe
  scores_all <- data.frame(
    Cell = cell_ids,
    stringsAsFactors = FALSE
  )
  
  # Process each column in reference data
  for (j in seq_len(ncol(reference_data))) {

    prog_data_sub <- as.data.frame(reference_data[, j, drop = FALSE])
    meta_z_label <- colnames(prog_data_sub)

    # Filter by z-score cutoff and remove NAs
    prog_data_sub <- na.omit(prog_data_sub)
    prog_data_sub <- prog_data_sub[
      abs(prog_data_sub[, meta_z_label]) > z_score_cutoff, ,
      drop = FALSE
    ]

    if (nrow(prog_data_sub) == 0) {
      warning(glue::glue("No genes pass z-score cutoff of {z_score_cutoff} for {meta_z_label}"))
      next
    }

    # Match genes between expression and reference
    common_genes <- intersect(rownames(expression_matrix), rownames(prog_data_sub))

    if (length(common_genes) == 0) {
      warning(glue::glue("No common genes found between expression and reference data for {meta_z_label}"))
      next
    }

    expression_data <- expression_matrix[common_genes, , drop = FALSE]
    prog_data_sub <- prog_data_sub[common_genes, , drop = FALSE]

    if (verbose) {
      cat(glue::glue("{length(common_genes)} genes used for scoring against {meta_z_label}\n"))
      cat("\nCalculating scores...\n")
    }

    # Weighted sum: higher score = worse prognosis when reference has positive z = adverse.
    # Compute_scores expects a named numeric vector with gene IDs in names().
    score_vector <- compute_scores(
      expression_data = expression_data,
      prognostic_scores = stats::setNames(prog_data_sub[, 1], rownames(prog_data_sub)),
      pseudobulk = pseudobulk,
      verbose = verbose
    )

    # Add to results
    score_col_name <- paste0("weighted_sum_score_", meta_z_label)
    scores_all[[score_col_name]] <- score_vector
    
    if (verbose) {
      cat(glue::glue("Completed scoring for {meta_z_label}
"))
    }
  }
  
  # Set rownames and remove Cell column
  rownames(scores_all) <- scores_all$Cell
  scores_all$Cell <- NULL
  
  return(scores_all)
}


#' Compute Score Vectors
#'
#' Compute the weighted sum of expression with the reference meta-z vector:
#' score(cell) = sum over genes of \code{expression[gene, cell] * meta_z[gene]}.
#' Implemented via cross product for efficiency: \code{crossprod(meta_z, expression)}
#' yields one value per cell. Higher score = worse prognosis (adverse).
#'
#' This helper explicitly aligns gene IDs between \code{expression_data} and
#' \code{prognostic_scores} so that the cross-product is always computed with
#' a consistent gene order.
#'
#' @keywords internal
compute_scores <- function(expression_data,
                           prognostic_scores,
                           pseudobulk = FALSE,
                           verbose = TRUE) {
  
  # Ensure expression_data has gene IDs as rownames and prognostic_scores is named
  if (is.null(rownames(expression_data))) {
    stop("expression_data must have rownames set to gene IDs")
  }
  if (is.null(names(prognostic_scores))) {
    stop("prognostic_scores must be a named numeric vector with gene IDs in names()")
  }
  
  # Align by common gene IDs
  common_genes <- intersect(rownames(expression_data), names(prognostic_scores))
  
  if (length(common_genes) == 0L) {
    stop("No overlapping genes between expression_data rownames and prognostic_scores names")
  }
  
  # Optionally warn if overlap is very small
  if (length(common_genes) < 50L) {
    warning("Fewer than 50 overlapping genes between expression_data and prognostic_scores")
  }
  
  # Use a consistent, explicit ordering so row i of expression_data and
  # element i of prognostic_scores always refer to the same gene (cross product is correct).
  common_genes <- sort(common_genes)
  expression_data <- expression_data[common_genes, , drop = FALSE]
  prognostic_scores <- prognostic_scores[common_genes]
  stopifnot(
    identical(rownames(expression_data), names(prognostic_scores)),
    identical(rownames(expression_data), common_genes)
  )

  # For non-pseudobulk or when matrix is small, use vectorized crossprod
  if (!pseudobulk || ncol(expression_data) < 100) {

    # crossprod(prognostic_scores, expression_data): (1 x n_genes) %*% (n_genes x n_cells) = (1 x n_cells)
    if (inherits(expression_data, "sparseMatrix")) {
      raw_sum <- as.numeric(
        Matrix::crossprod(prognostic_scores, expression_data)
      )
    } else {
      raw_sum <- as.numeric(
        crossprod(prognostic_scores, expression_data)
      )
    }
    # Higher raw_sum = higher expression of adverse genes (positive z) = worse prognosis
    score_vector <- raw_sum
    
  } else {
    # nocov start - pseudobulk path (large matrices, rarely used in tests)
    # For large pseudobulk data, compute iteratively with progress bar
    if (verbose && requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        total = ncol(expression_data),
        format = "  [:bar] :percent eta: :eta",
        force = TRUE
      )
    }
    
    score_vector <- numeric(ncol(expression_data))
    
    for (j in seq_len(ncol(expression_data))) {
      exp_vector <- expression_data[, j]
      raw_sum_j <- sum(prognostic_scores * exp_vector, na.rm = TRUE)
      score_vector[j] <- raw_sum_j
      
      if (verbose && exists("pb")) {
        pb$tick()
      }
    }
    # nocov end
  }
  
  return(score_vector)
}


#' Calculate Z-Scores
#'
#' Normalize scores to z-scores
#'
#' @param scores Numeric vector of scores
#' @return Numeric vector of z-scores
#'
#' @export
normalize_scores <- function(scores) {
  
  if (!is.numeric(scores)) {
    stop("scores must be numeric")
  }
  
  scores_mean <- mean(scores, na.rm = TRUE)
  scores_sd <- sd(scores, na.rm = TRUE)
  
  if (scores_sd == 0) {
    warning("Standard deviation is 0, returning original scores")
    return(scores)
  }
  
  z_scores <- (scores - scores_mean) / scores_sd
  
  return(z_scores)
}


#' Add Scores to Seurat Object
#'
#' Convenience function to add scores directly to Seurat metadata. Use the
#' **original** (non-subset) Seurat object so that all genes are retained for
#' downstream analyses (e.g. cell type marker gene analysis).
#'
#' @param seurat_obj Seurat object (the same full object passed to \code{PhenoMap})
#' @param scores Data.frame of scores (output from PhenoMap)
#' @param prefix Prefix for metadata column names (default: "")
#'
#' @return Seurat object with scores added to metadata
#'
#' @export
add_scores_to_seurat <- function(seurat_obj, scores, prefix = "") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required")
  }
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  # Match cells
  common_cells <- intersect(rownames(scores), colnames(seurat_obj))
  
  if (length(common_cells) == 0) {
    stop("No common cells found between scores and Seurat object")
  }
  
  if (length(common_cells) < nrow(scores)) {
    warning(glue::glue(
      "{nrow(scores) - length(common_cells)} cells in scores not found in Seurat object"
    ))
  }
  
  # Add each score column to metadata
  for (col in colnames(scores)) {
    col_name <- if (prefix != "") paste0(prefix, col) else col
    seurat_obj@meta.data[common_cells, col_name] <- scores[common_cells, col]
  }
  
  message(glue::glue("Added {ncol(scores)} score column(s) to Seurat metadata"))
  
  return(seurat_obj)
}


#' Add Scores to SingleCellExperiment Object
#'
#' Convenience function to add scores to SCE colData. Use the **original**
#' (non-subset) SCE object so that all genes are retained for downstream
#' analyses (e.g. marker genes).
#'
#' @param sce_obj SingleCellExperiment object (the same full object passed to \code{PhenoMap})
#' @param scores Data.frame of scores (output from PhenoMap)
#' @param prefix Prefix for colData column names (default: "")
#'
#' @return SingleCellExperiment object with scores in colData
#'
#' @export
add_scores_to_sce <- function(sce_obj, scores, prefix = "") {
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required")
  }
  
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required")
  }
  
  if (!inherits(sce_obj, "SingleCellExperiment")) {
    stop("sce_obj must be a SingleCellExperiment object")
  }
  
  # Match cells
  common_cells <- intersect(rownames(scores), colnames(sce_obj))
  
  if (length(common_cells) == 0) {
    stop("No common cells found between scores and SCE object")
  }
  
  # Add to colData
  for (col in colnames(scores)) {
    col_name <- if (prefix != "") paste0(prefix, col) else col
    SummarizedExperiment::colData(sce_obj)[common_cells, col_name] <- scores[common_cells, col]
  }
  
  message(glue::glue("Added {ncol(scores)} score column(s) to SCE colData"))
  
  return(sce_obj)
}
