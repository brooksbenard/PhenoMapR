#' Process Expression Input
#'
#' Convert various expression data formats to matrix
#'
#' @keywords internal
process_expression_input <- function(expression,
                                     pseudobulk = FALSE,
                                     group_by = NULL,
                                     assay = NULL,
                                     slot = "data",
                                     genes_to_extract = NULL,
                                     verbose = TRUE) {

  input_type <- detect_input_type(expression)

  if (verbose) {
    message(glue::glue("Detected input type: {input_type}"))
  }

  result <- switch(input_type,
    "matrix" = process_matrix(expression),
    "seurat" = process_seurat(expression, pseudobulk, group_by, assay, slot, genes_to_extract),
    "seurat_spatial" = process_seurat(expression, pseudobulk, group_by, assay, slot = "counts", genes_to_extract = genes_to_extract),
    "sce" = process_sce(expression, pseudobulk, group_by, assay, genes_to_extract),
    "spatial_experiment" = process_spatial_experiment(expression, pseudobulk, group_by, assay, genes_to_extract),
    "anndata" = process_anndata(expression, pseudobulk, group_by),
    stop("Unsupported input type")
  )

  result$input_type <- input_type
  return(result)
}


#' Detect Input Type
#'
#' @keywords internal
detect_input_type <- function(obj) {

  if (is.matrix(obj) || is.data.frame(obj)) {
    return("matrix")
  }
  # Sparse matrices from Matrix package (e.g. dgCMatrix from Read10X_h5)
  if (inherits(obj, "Matrix")) {
    return("matrix")
  }

  if (inherits(obj, "Seurat")) {
    # Check if it's spatial
    if ("Spatial" %in% names(obj@assays)) {
      return("seurat_spatial")
    }
    return("seurat")
  }

  if (inherits(obj, "SingleCellExperiment")) {
    return("sce")
  }

  # nocov start - optional SpatialExperiment (not in default test deps)
  if (inherits(obj, "SpatialExperiment")) {
    return("spatial_experiment")
  }
  # nocov end

  # nocov start - optional AnnData/reticulate (not in default test deps)
  if (inherits(obj, "python.builtin.object")) {
    if (reticulate::py_has_attr(obj, "X") && reticulate::py_has_attr(obj, "obs")) {
      return("anndata")
    }
  }
  # nocov end

  stop("Unable to detect input type. Supported: matrix, Seurat, SingleCellExperiment, SpatialExperiment, AnnData")
}


#' Process Matrix Input
#'
#' @keywords internal
process_matrix <- function(mat) {

  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }

  # Ensure genes are rows, samples are columns
  if (is.null(rownames(mat))) {
    stop("Matrix must have gene names as rownames")
  }

  list(
    matrix = mat,
    cell_names = colnames(mat),
    gene_names = rownames(mat)
  )
}


#' Process Seurat Object
#'
#' @keywords internal
process_seurat <- function(obj, pseudobulk, group_by, assay, slot, genes_to_extract = NULL) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required but not installed")
  }

  # Determine assay
  if (is.null(assay)) {
    assay <- if ("Spatial" %in% names(obj@assays)) "Spatial" else "RNA"
  }

  # Map slot names to layer names for Seurat v5
  layer_map <- c(
    "data" = "data",
    "counts" = "counts",
    "scale.data" = "scale.data"
  )

  layer_name <- if (!is.null(slot) && slot %in% names(layer_map)) {
    layer_map[slot]
  } else {
    "data"
  }

  # Handle pseudobulk
  if (pseudobulk) {
    if (is.null(group_by)) {
      stop("group_by must be specified for pseudobulk aggregation")
    }

    if (!group_by %in% colnames(obj@meta.data)) {
      stop(glue::glue("'{group_by}' not found in Seurat metadata"))
    }

    # AggregateExpression works with both v4 and v5
    agg_expr <- Seurat::AggregateExpression(
      obj,
      assays = assay,
      group.by = group_by
    )

    expr_matrix <- as.matrix(agg_expr[[assay]])

  } else {
    # Get assay data from the full object (do not subset the Seurat object by
    # genes, to avoid copying spatial images and triggering VisiumV1/etc.
    # validity errors when the object was saved with a different Seurat version).
    assay_data <- tryCatch({
      Seurat::GetAssayData(obj, assay = assay, layer = layer_name)
    }, error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = slot)
    })

    if (length(genes_to_extract) > 0) {
      avail_genes <- rownames(assay_data)
      genes_use <- intersect(genes_to_extract, avail_genes)
      if (length(genes_use) > 0) {
        assay_data <- assay_data[genes_use, , drop = FALSE]
      }
    }

    # Avoid as.matrix on any Matrix class to prevent sparse->dense (9+ GiB for large objects)
    if (inherits(assay_data, "Matrix") || inherits(assay_data, "sparseMatrix")) {
      expr_matrix <- assay_data
    } else if (!is.matrix(assay_data)) {
      expr_matrix <- as.matrix(assay_data)
    } else {
      expr_matrix <- assay_data
    }
  }

  list(
    matrix = expr_matrix,
    cell_names = colnames(expr_matrix),
    gene_names = rownames(expr_matrix),
    assay = assay,
    slot = slot
  )
}


#' Process SingleCellExperiment Object
#'
#' @keywords internal
process_sce <- function(obj, pseudobulk, group_by, assay, genes_to_extract = NULL) {

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required but not installed")
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required but not installed")
  }

  # Determine which assay to use
  if (is.null(assay)) {
    assay <- "logcounts"
    if (!assay %in% SummarizedExperiment::assayNames(obj)) {
      assay <- SummarizedExperiment::assayNames(obj)[1]
    }
  }

  if (pseudobulk) {
    if (is.null(group_by)) {
      stop("group_by must be specified for pseudobulk aggregation")
    }

    if (!group_by %in% colnames(SummarizedExperiment::colData(obj))) {
      stop(glue::glue("'{group_by}' not found in SCE colData"))
    }

    # Aggregate by group
    groups <- SummarizedExperiment::colData(obj)[[group_by]]
    expr_matrix <- SummarizedExperiment::assay(obj, assay)

    # Sum across groups
    agg_matrix <- sapply(unique(groups), function(g) {
      cells <- which(groups == g)
      Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
    })

    expr_matrix <- as.matrix(agg_matrix)

  } else {
    assay_data <- SummarizedExperiment::assay(obj, assay)
    # Subset only for extraction when genes_to_extract set; never return subset.
    # Caller adds scores to original (full) object so all genes remain for analyses.
    if (length(genes_to_extract) > 0) {
      avail_genes <- rownames(assay_data)
      genes_use <- intersect(genes_to_extract, avail_genes)
      if (length(genes_use) > 0) {
        assay_data <- assay_data[genes_use, , drop = FALSE]
      }
    }
    if (inherits(assay_data, "Matrix") || inherits(assay_data, "sparseMatrix")) {
      expr_matrix <- assay_data
    } else {
      expr_matrix <- as.matrix(assay_data)
    }
  }

  list(
    matrix = expr_matrix,
    cell_names = colnames(expr_matrix),
    gene_names = rownames(expr_matrix),
    assay = assay
  )
}


#' Process SpatialExperiment Object
#'
#' @keywords internal
# nocov start - optional SpatialExperiment
process_spatial_experiment <- function(obj, pseudobulk, group_by, assay, genes_to_extract = NULL) {

  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment package required but not installed")
  }

  # SpatialExperiment inherits from SCE, so use same processing
  process_sce(obj, pseudobulk, group_by, assay, genes_to_extract)
}
# nocov end


#' Process AnnData Object
#'
#' @keywords internal
# nocov start - optional AnnData/reticulate
process_anndata <- function(obj, pseudobulk, group_by) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package required for AnnData objects")
  }

  # Extract expression matrix (use .X which is usually normalized)
  expr_matrix <- as.matrix(obj$X)

  # Set rownames and colnames
  if (reticulate::py_has_attr(obj, "var_names")) {
    rownames(expr_matrix) <- reticulate::py_to_r(obj$var_names)
  }
  if (reticulate::py_has_attr(obj, "obs_names")) {
    colnames(expr_matrix) <- reticulate::py_to_r(obj$obs_names)
  }

  # Transpose (AnnData stores as cells x genes)
  expr_matrix <- t(expr_matrix)

  if (pseudobulk) {
    if (is.null(group_by)) {
      stop("group_by must be specified for pseudobulk aggregation")
    }

    obs_df <- reticulate::py_to_r(obj$obs)

    if (!group_by %in% colnames(obs_df)) {
      stop(glue::glue("'{group_by}' not found in AnnData obs"))
    }

    groups <- obs_df[[group_by]]

    # Aggregate
    agg_matrix <- sapply(unique(groups), function(g) {
      cells <- which(groups == g)
      Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
    })

    expr_matrix <- as.matrix(agg_matrix)
  }

  list(
    matrix = expr_matrix,
    cell_names = colnames(expr_matrix),
    gene_names = rownames(expr_matrix)
  )
}
# nocov end


#' Validate Expression Matrix
#'
#' @keywords internal
validate_expression_matrix <- function(mat) {

  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    stop("Expression data must be a matrix or Matrix (sparse)")
  }

  if (is.null(rownames(mat))) {
    stop("Expression matrix must have gene names as rownames")
  }

  if (is.null(colnames(mat))) {
    warning("Expression matrix has no column names, generating generic names")
    colnames(mat) <- paste0("Cell_", seq_len(ncol(mat)))
  }

  # Skip expensive any() on large matrices to avoid memory blow-up (sparse or dense)
  n_elem <- as.numeric(nrow(mat)) * as.numeric(ncol(mat))
  if (n_elem <= 1e6) {
    if (any(is.na(mat))) {
      warning("Expression matrix contains NA values")
    }
    if (any(mat < 0, na.rm = TRUE)) {
      warning("Expression matrix contains negative values")
    }
  }

  invisible(TRUE)
}
