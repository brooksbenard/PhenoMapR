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
  
  # Initialize scores dataframe
  scores_all <- data.frame(
    Cell = colnames(expression_matrix),
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
      cat(glue::glue("{length(common_genes)} genes used for scoring against {meta_z_label}
"))
      cat("Calculating scores...
")
    }
    
    # Calculate scores
    score_vector <- compute_scores(
      expression_data = expression_data,
      prognostic_scores = prog_data_sub[, 1],
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
#' Efficiently compute weighted sum using matrix operations
#'
#' @keywords internal
compute_scores <- function(expression_data,
                          prognostic_scores,
                          pseudobulk = FALSE,
                          verbose = TRUE) {
  
  # For non-pseudobulk or when matrix is small, use vectorized crossprod
  if (!pseudobulk || ncol(expression_data) < 100) {
    
    # Convert to dense matrix if sparse
    if (inherits(expression_data, "sparseMatrix")) {
      score_vector <- as.numeric(
        Matrix::crossprod(prognostic_scores, expression_data)
      )
    } else {
      score_vector <- as.numeric(
        crossprod(prognostic_scores, expression_data)
      )
    }
    
  } else {
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
      score_vector[j] <- sum(prognostic_scores * exp_vector, na.rm = TRUE)
      
      if (verbose && exists("pb")) {
        pb$tick()
      }
    }
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
#' Convenience function to add scores directly to Seurat metadata
#'
#' @param seurat_obj Seurat object
#' @param scores Data.frame of scores (output from score_expression)
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
#' Convenience function to add scores to SCE colData
#'
#' @param sce_obj SingleCellExperiment object
#' @param scores Data.frame of scores
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
