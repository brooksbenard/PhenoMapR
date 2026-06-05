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
    "matrix" = process_matrix(expression, verbose = verbose),
    "seurat" = process_seurat(expression, pseudobulk, group_by, assay, slot, genes_to_extract),
    # nocov start - spatial Seurat (slot=counts) not in default tests
    "seurat_spatial" = process_seurat(expression, pseudobulk, group_by, assay, slot = "counts", genes_to_extract = genes_to_extract),
    # nocov end
    "sce" = process_sce(expression, pseudobulk, group_by, assay, genes_to_extract),
    "spatial_experiment" = process_spatial_experiment(expression, pseudobulk, group_by, assay, genes_to_extract),
    "anndata" = process_anndata(expression, pseudobulk, group_by, genes_to_extract),
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


#' Validate expression matrix axes and gene ID format
#'
#' Ensures genes are rows and samples/cells/spots are columns (heuristic: orders of
#' magnitude more genes than samples). If the matrix appears transposed, it is
#' transposed and a message is printed. Also checks for ENSG-style rownames and
#' warns the user to convert to HUGO symbols.
#'
#' @param mat Matrix (or matrix-coercible object) with rownames and colnames.
#' @param verbose If TRUE, print a message when transposing.
#' @return The matrix, possibly transposed; rownames must be gene IDs, colnames sample IDs.
#' @keywords internal
validate_expression_axes_and_ids <- function(mat, verbose = TRUE) {
  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  # Heuristic: expect orders of magnitude more genes than samples (e.g. 10x or more).
  # If instead samples >> genes, treat as samples x genes and transpose.
  if (n_genes > 0 && n_samples > 0 && !is.null(rownames(mat)) && !is.null(colnames(mat))) {
    if (n_samples > 10 * n_genes) {
      if (isTRUE(verbose)) {
        message(
          "Expression format transposed so rows are gene IDs and columns are samples ",
          "(number of samples is much greater than number of rows)."
        )
      }
      mat <- t(mat)
      # After transpose: rows = genes, columns = cells/samples
      n_genes <- nrow(mat)
      n_samples <- ncol(mat)
    }
    # Sanity check: after correction we should have many more genes than samples
    if (n_genes > 0 && n_samples > 0 && n_genes < 5 * n_samples && isTRUE(verbose)) {
      message(
        "Note: number of rows is not much larger than number of columns. ",
        "Please ensure rows are gene IDs and columns are samples/cells/spots."
      )
    }
  }

  gene_ids <- rownames(mat)
  if (!is.null(gene_ids) && length(gene_ids) > 0) {
    prop_ensg <- mean(grepl("^ENSG\\d", gene_ids))
    if (is.finite(prop_ensg) && prop_ensg > 0.5) {
      warning(
        "Gene names appear to be Ensembl/ENSG IDs (e.g., 'ENSG...'). ",
        "Please convert to HUGO gene symbols before using PhenoMapR (e.g., with biomaRt or AnnotationDbi)."
      )
    }
  }

  mat
}


#' Process Matrix Input
#'
#' @keywords internal
process_matrix <- function(mat, verbose = TRUE) {

  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }

  mat <- validate_expression_axes_and_ids(mat, verbose = verbose)

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

  # Determine assay. Honour the object's DefaultAssay() first so that
  # SCTransform-normalized objects (default assay "SCT") and any
  # other multi-assay Seurat workflows are served by the assay the
  # user actually intends. Falls back to the historical
  # Spatial-vs-RNA priority for objects without a usable default.
  if (is.null(assay)) {
    avail <- names(obj@assays)
    assay <- tryCatch(Seurat::DefaultAssay(obj),
                      error = function(e) NULL)
    if (is.null(assay) || !nzchar(assay) || !assay %in% avail) {
      assay <- if ("Spatial" %in% avail) "Spatial"
               else if ("RNA" %in% avail) "RNA"
               else if ("SCT" %in% avail) "SCT"
               else avail[1L]
    }
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
  # nocov start - pseudobulk with Seurat not exercised in default tests
  if (pseudobulk) {
    if (is.null(group_by)) {
      stop("group_by must be specified for pseudobulk aggregation")
    }

    if (!group_by %in% colnames(obj@meta.data)) {
      stop(glue::glue("'{group_by}' not found in Seurat metadata"))
    }

    grp_vec <- obj@meta.data[[group_by]]
    grp_u <- unique(stats::na.omit(as.character(grp_vec)))
    if (length(grp_u) < 1L) {
      stop(glue::glue("No non-NA values in metadata column '{group_by}' for pseudobulk aggregation"))
    }

    if (length(grp_u) == 1L) {
      # Single pseudobulk group: Seurat::AggregateExpression can fail internally
      # (e.g. colnames<- on a 1D category matrix). Sum across cells = same RNA
      # aggregation as one group in AggregateExpression(counts).
      assay_data <- .get_assay_data_compat(obj, assay = assay, slot = layer_name)
      if (inherits(assay_data, "Matrix") || inherits(assay_data, "sparseMatrix")) {
        summed <- Matrix::rowSums(assay_data)
      } else {
        summed <- rowSums(assay_data)
      }
      expr_matrix <- matrix(
        as.numeric(summed),
        ncol = 1L,
        dimnames = list(rownames(assay_data), grp_u)
      )
    } else {
      agg_expr <- Seurat::AggregateExpression(
        obj,
        assays = assay,
        group.by = group_by
      )
      expr_matrix <- as.matrix(agg_expr[[assay]])
      if (!is.matrix(expr_matrix)) {
        cn <- tryCatch(colnames(agg_expr[[assay]]), error = function(e) NULL)
        if (is.null(cn) || length(cn) < 1L) cn <- grp_u[1]
        rn <- rownames(agg_expr[[assay]])
        if (is.null(rn)) rn <- names(expr_matrix)
        expr_matrix <- matrix(
          as.numeric(expr_matrix),
          ncol = 1L,
          dimnames = list(rn, cn[1])
        )
      }
    }

  } else {
  # nocov end
    # Get assay data from the full object (do not subset the Seurat object by
    # genes, to avoid copying spatial images and triggering VisiumV1/etc.
    # validity errors when the object was saved with a different Seurat version).
    assay_data <- .get_assay_data_compat(obj, assay = assay, slot = layer_name)

    # nocov start - genes_to_extract subset (PhenoMap sets it for Seurat/SCE; tests use small objects)
    if (length(genes_to_extract) > 0) {
      avail_genes <- rownames(assay_data)
      genes_use <- intersect(genes_to_extract, avail_genes)
      if (length(genes_use) > 0) {
        assay_data <- assay_data[genes_use, , drop = FALSE]
      }
    }
    # nocov end

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

  # nocov start - pseudobulk with SCE not exercised in default tests
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
    grp_uniq <- unique(groups)

    # Sum across groups (cbind so a single group stays a matrix, not a vector)
    agg_cols <- lapply(grp_uniq, function(g) {
      cells <- which(groups == g)
      Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
    })
    agg_matrix <- do.call(cbind, agg_cols)
    colnames(agg_matrix) <- as.character(grp_uniq)
    rownames(agg_matrix) <- rownames(expr_matrix)

    expr_matrix <- as.matrix(agg_matrix)

  } else {
  # nocov end
    assay_data <- SummarizedExperiment::assay(obj, assay)
    # nocov start - genes_to_extract subset
    if (length(genes_to_extract) > 0) {
      avail_genes <- rownames(assay_data)
      genes_use <- intersect(genes_to_extract, avail_genes)
      if (length(genes_use) > 0) {
        assay_data <- assay_data[genes_use, , drop = FALSE]
      }
    }
    # nocov end
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
#' Convert a Python \pkg{anndata} object into a genes x cells expression matrix
#' that the rest of PhenoMapR understands. Optimized for very large objects:
#' \itemize{
#'   \item When \code{genes_to_extract} is supplied (e.g. the reference genes
#'         that pass the z-score cutoff in \code{PhenoMap()}), the AnnData is
#'         subset \emph{on the Python side} before any data is copied into R.
#'         This is the single biggest memory win for multi-GB \code{.h5ad}
#'         files: only a few hundred to a few thousand genes typically pass
#'         the cutoff, so we transfer ~1-3\% of the matrix instead of all of
#'         it.
#'   \item scipy-sparse \code{.X} is reinterpreted directly as a
#'         \code{dgCMatrix} in genes x cells orientation by treating the
#'         native CSR-of-(cellsxgenes) storage as CSC-of-(genesxcells); this
#'         avoids \code{Matrix::t()} and the doubling of memory it would
#'         otherwise cause.
#' }
#'
#' @keywords internal
# nocov start - optional AnnData/reticulate
process_anndata <- function(obj, pseudobulk, group_by, genes_to_extract = NULL) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package required for AnnData objects")
  }

  expr_matrix <- .anndata_X_to_genes_cells(obj, gene_subset = genes_to_extract)

  if (pseudobulk) {
    if (is.null(group_by)) {
      stop("group_by must be specified for pseudobulk aggregation")
    }

    obs_df <- .anndata_obs_df(obj)
    if (is.null(obs_df) || !group_by %in% colnames(obs_df)) {
      stop(glue::glue("'{group_by}' not found in AnnData obs"))
    }

    groups <- as.character(obs_df[[group_by]])
    grp_uniq <- unique(stats::na.omit(groups))

    # Sum per group while keeping sparsity (rowSums on dgCMatrix returns dense)
    agg_cols <- lapply(grp_uniq, function(g) {
      cells <- which(groups == g)
      Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
    })
    agg_matrix <- do.call(cbind, agg_cols)
    colnames(agg_matrix) <- as.character(grp_uniq)
    rownames(agg_matrix) <- rownames(expr_matrix)
    expr_matrix <- as.matrix(agg_matrix)
  }

  list(
    matrix = expr_matrix,
    cell_names = colnames(expr_matrix),
    gene_names = rownames(expr_matrix)
  )
}

#' Convert AnnData.X into a genes x cells Matrix
#'
#' Returns a \code{dgCMatrix} (sparse) when \code{adata.X} is a scipy sparse
#' matrix, or a regular dense matrix otherwise. Always genes x cells (i.e.
#' the transpose of AnnData's native cells x genes layout) with
#' \code{rownames = var_names} and \code{colnames = obs_names}.
#'
#' Two memory optimisations versus the naive \code{as.matrix(adata$X)} approach:
#'
#' \enumerate{
#'   \item When \code{gene_subset} is supplied, the AnnData is sliced on the
#'         var axis \emph{in Python} (\code{adata[:, mask].X}) before any data
#'         is copied into R. This is the dominant cost reduction for multi-GB
#'         AnnData objects: only the genes that actually contribute to the
#'         score are transferred.
#'   \item For scipy-sparse \code{.X}, the AnnData native CSR storage of a
#'         (n_obs x n_vars) matrix is the same as the CSC storage of the
#'         transposed (n_vars x n_obs) matrix. We reuse \code{indices},
#'         \code{indptr} and \code{data} arrays directly to build a
#'         \code{dgCMatrix} in genes x cells orientation, with no extra
#'         allocation for a transpose pass.
#' }
#'
#' @keywords internal
.anndata_X_to_genes_cells <- function(obj, gene_subset = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package required for AnnData objects")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required for sparse AnnData conversion")
  }

  # ----- 1. Optionally subset on the Python side --------------------------
  X <- obj$X
  obs_names_use <- .anndata_obs_names(obj)
  var_names_use <- .anndata_var_names(obj)

  # rename_to: when set, we replace rownames of the final matrix with these
  # values (used when we matched via a .var symbol column rather than the
  # raw var_names -- see the Ensembl-ID fallback below).
  rename_to <- NULL

  if (!is.null(gene_subset) && length(gene_subset) > 0L) {
    gene_subset_c <- as.character(gene_subset)
    keep_direct   <- intersect(gene_subset_c, var_names_use)

    if (length(keep_direct) > 0L && length(keep_direct) < length(var_names_use)) {
      # Standard happy path: var_names overlap the reference, subset directly.
      sub <- .anndata_subset_var(obj, keep_direct)
      X             <- sub$X
      var_names_use <- sub$var_names
      obs_names_use <- sub$obs_names

    } else if (length(keep_direct) == 0L) {
      # Zero overlap on var_names. This is overwhelmingly because the file
      # stores Ensembl IDs (ENSG...) as var_names with HUGO symbols hiding in
      # a column of .var. Try to recover automatically.
      sym <- .anndata_find_symbol_column(obj, gene_subset_c)
      if (!is.null(sym) && length(sym$var_names) > 0L) {
        message(sprintf(
          "AnnData.var_names do not overlap the reference (looks like a non-HUGO naming scheme); matched %d genes via AnnData.var$%s and will use those symbols as rownames.",
          length(sym$var_names), sym$column
        ))
        sub <- .anndata_subset_var(obj, sym$var_names)
        X             <- sub$X
        var_names_use <- sub$var_names
        obs_names_use <- sub$obs_names
        # The Python subset preserves the original var_names order, which
        # may differ from the order we found symbols in. Re-look up symbols
        # by the actual var_names returned.
        rename_to <- as.character(setNames(sym$symbols, sym$var_names)[var_names_use])

      } else {
        stop(
          "None of the requested genes are present in AnnData.var_names, ",
          "and no recognisable gene-symbol column was found in AnnData.var. ",
          "First few var_names are: ",
          paste(head(var_names_use, 5L), collapse = ", "), ". ",
          "Either rename AnnData.var_names to HUGO symbols before saving the .h5ad ",
          "(e.g. with `adata.var_names = adata.var['gene_symbols']`), or add a ",
          "column like 'gene_symbol' / 'feature_name' / 'Symbol' to AnnData.var ",
          "with HUGO symbols -- PhenoMapR will pick it up automatically."
        )
      }
    }
    # If keep_direct == length(var_names_use) we don't subset at all
    # (every gene is requested) and fall through to the full conversion.
  }

  # ----- 2. Build the R matrix in genes x cells orientation ---------------
  # Already an R-side object (e.g. user passed convert = TRUE on import)
  if (is.matrix(X) || inherits(X, "Matrix")) {
    m <- if (inherits(X, "Matrix")) Matrix::t(X) else t(X)
    if (length(var_names_use) == nrow(m)) rownames(m) <- var_names_use
    if (length(obs_names_use) == ncol(m)) colnames(m) <- obs_names_use
    return(m)
  }

  scipy_sparse <- tryCatch(
    reticulate::import("scipy.sparse", convert = FALSE),
    error = function(e) NULL
  )
  is_sparse <- !is.null(scipy_sparse) && isTRUE(tryCatch(
    reticulate::py_to_r(scipy_sparse$issparse(X)),
    error = function(e) FALSE
  ))

  if (is_sparse) {
    # Convert to CSR (cheap if already CSR -- AnnData typically stores CSR).
    # Then reuse the CSR(n_obs, n_vars) storage as CSC(n_vars, n_obs):
    # CSR.indices (column indices of X) === CSC.i (row indices of t(X))
    # CSR.indptr  (row pointers of X)    === CSC.p (column pointers of t(X))
    # CSR.data    (values)               === CSC.x (values)
    # Resulting dgCMatrix is genes x cells with no Matrix::t() pass.
    csr <- X$tocsr()
    shape   <- as.integer(reticulate::py_to_r(csr$shape))     # (n_obs, n_vars)
    indices <- as.integer(reticulate::py_to_r(csr$indices))
    indptr  <- as.integer(reticulate::py_to_r(csr$indptr))
    values  <- as.numeric(reticulate::py_to_r(csr$data))
    m <- Matrix::sparseMatrix(
      i      = indices,
      p      = indptr,
      x      = values,
      dims   = c(shape[2L], shape[1L]),     # transpose: n_vars x n_obs
      index1 = FALSE
    )
  } else {
    # Dense path: pull the (n_obs x n_vars) array and transpose to genes x cells.
    dense <- reticulate::py_to_r(X)
    if (!is.matrix(dense) && !inherits(dense, "Matrix")) {
      dense <- as.matrix(dense)
    }
    m <- t(dense)
  }

  # If we matched via a .var symbol column, rownames are the human-readable
  # HUGO symbols rather than the (Ensembl) var_names. Otherwise the matrix
  # is labelled by var_names exactly as in the AnnData.
  effective_rownames <- if (!is.null(rename_to)) rename_to else var_names_use
  if (length(effective_rownames) == nrow(m)) {
    rownames(m) <- effective_rownames
  }
  if (length(obs_names_use) == ncol(m)) colnames(m) <- obs_names_use
  m
}

#' Find a HUGO-symbol-style column in AnnData.var that overlaps the reference
#'
#' Scans \code{AnnData.var} for the most common gene-symbol column conventions
#' (\code{gene_symbol}, \code{gene_symbols}, \code{feature_name}, \code{Symbol},
#' \code{gene_short_name}, \code{hgnc_symbol}, \code{gene_name}) and returns
#' the first column with any overlap against \code{gene_subset}.
#'
#' @return A list with \code{column}, \code{var_names} (the AnnData
#'   \code{var_names} whose symbol matched, suitable for passing to
#'   \code{.anndata_subset_var}) and \code{symbols} (the matched HUGO symbols
#'   in the same order). Returns \code{NULL} when no usable column is found.
#'
#' @keywords internal
.anndata_find_symbol_column <- function(obj, gene_subset) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  if (!reticulate::py_has_attr(obj, "var")) return(NULL)
  var_df <- tryCatch(reticulate::py_to_r(obj$var), error = function(e) NULL)
  if (is.null(var_df)) return(NULL)
  if (!is.data.frame(var_df)) {
    var_df <- tryCatch(as.data.frame(var_df, stringsAsFactors = FALSE),
                       error = function(e) NULL)
    if (is.null(var_df)) return(NULL)
  }
  rn <- rownames(var_df)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(var_df))))) {
    rn <- .anndata_var_names(obj, n = nrow(var_df))
  }

  candidates <- c(
    "gene_symbol", "gene_symbols", "feature_name", "Symbol", "symbol",
    "gene_short_name", "hgnc_symbol", "gene_name", "Gene", "gene"
  )
  candidates <- intersect(candidates, colnames(var_df))
  if (!length(candidates)) return(NULL)

  for (col in candidates) {
    symbols <- as.character(var_df[[col]])
    ok <- which(!is.na(symbols) & nzchar(symbols) & symbols %in% gene_subset)
    if (length(ok) == 0L) next
    # If duplicate symbols exist, keep the first occurrence so the resulting
    # matrix has unique rownames downstream.
    dedup <- !duplicated(symbols[ok])
    ok <- ok[dedup]
    return(list(
      column    = col,
      var_names = rn[ok],
      symbols   = symbols[ok]
    ))
  }
  NULL
}

#' Subset an AnnData on the var (gene) axis in Python
#'
#' Uses a tiny Python helper defined in \code{__main__} to mask AnnData on
#' \code{var_names} and return the subset's \code{.X} plus the new
#' \code{var_names} / \code{obs_names}. All return values are kept as Python
#' objects so the caller can apply the most memory-efficient conversion.
#'
#' @keywords internal
.anndata_subset_var <- function(obj, gene_keep) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate required for AnnData subsetting")
  }
  # py_run_string is idempotent in __main__ -- re-defining is cheap so we
  # don't bother caching the module reference across calls.
  reticulate::py_run_string(
    paste(
      "def _phenomapr_subset_anndata_var(adata, gene_names):",
      "    keep_list = list(gene_names)",
      "    mask = adata.var_names.isin(keep_list)",
      "    sub = adata[:, mask]",
      "    return {",
      "        'X': sub.X,",
      "        'var_names': sub.var_names.tolist(),",
      "        'obs_names': sub.obs_names.tolist(),",
      "    }",
      sep = "\n"
    ),
    convert = FALSE
  )
  main <- reticulate::import_main(convert = FALSE)
  res <- main$`_phenomapr_subset_anndata_var`(obj, as.character(gene_keep))
  list(
    X         = reticulate::py_get_item(res, "X"),
    var_names = as.character(reticulate::py_to_r(reticulate::py_get_item(res, "var_names"))),
    obs_names = as.character(reticulate::py_to_r(reticulate::py_get_item(res, "obs_names")))
  )
}

#' Extract AnnData.obs as an R data.frame with cell IDs in rownames
#'
#' @keywords internal
.anndata_obs_df <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  if (!reticulate::py_has_attr(obj, "obs")) return(NULL)
  df <- tryCatch(reticulate::py_to_r(obj$obs), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(df))))) {
    ids <- .anndata_obs_names(obj, n = nrow(df))
    if (length(ids) == nrow(df)) rownames(df) <- ids
  }
  df
}

#' List AnnData.obsm keys whose value is at least 2 columns wide
#' @keywords internal
.anndata_obsm_keys <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(character(0))
  if (!reticulate::py_has_attr(obj, "obsm")) return(character(0))
  # Importantly, when the AnnData was imported with convert = FALSE,
  # iterate(obj$obsm$keys()) returns a list of Python str objects whose
  # R-side print is "<environment>". Converting via builtins.list() first
  # gives us a Python list of strs that py_to_r() turns into an R character
  # vector cleanly. We also fall back to iterate() with per-element py_to_r
  # in case builtins is unavailable.
  builtins <- tryCatch(
    reticulate::import_builtins(convert = FALSE),
    error = function(e) NULL
  )
  keys <- tryCatch({
    if (!is.null(builtins)) {
      as.character(reticulate::py_to_r(builtins$list(obj$obsm$keys())))
    } else {
      vapply(
        reticulate::iterate(obj$obsm$keys()),
        function(k) as.character(reticulate::py_to_r(k)),
        character(1L)
      )
    }
  }, error = function(e) character(0))
  if (!length(keys)) return(character(0))
  out <- character(0)
  for (k in keys) {
    arr <- tryCatch(reticulate::py_get_item(obj$obsm, k),
                    error = function(e) NULL)
    if (is.null(arr)) next
    shape <- tryCatch(reticulate::py_to_r(arr$shape),
                      error = function(e) NULL)
    if (!is.null(shape) && length(shape) >= 2L && shape[[2L]] >= 2L) {
      out <- c(out, k)
    }
  }
  out
}

#' Pull a 2D embedding from AnnData.obsm
#' @keywords internal
.anndata_obsm_array <- function(obj, key) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  if (!reticulate::py_has_attr(obj, "obsm")) return(NULL)
  arr <- tryCatch(reticulate::py_get_item(obj$obsm, key),
                  error = function(e) NULL)
  if (is.null(arr)) return(NULL)
  emb <- tryCatch(reticulate::py_to_r(arr), error = function(e) NULL)
  if (is.null(emb)) return(NULL)
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  ids <- .anndata_obs_names(obj, n = nrow(emb))
  if (length(ids) == nrow(emb)) rownames(emb) <- ids
  emb
}

#' AnnData obs_names / var_names with safe fallback
#' @keywords internal
.anndata_obs_names <- function(obj, n = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(if (is.null(n)) character(0) else paste0("Cell_", seq_len(n)))
  }
  ids <- tryCatch(
    as.character(reticulate::py_to_r(obj$obs_names$tolist())),
    error = function(e) NULL
  )
  if (is.null(ids) && !is.null(n)) ids <- paste0("Cell_", seq_len(n))
  if (is.null(ids)) character(0) else ids
}
#' @keywords internal
.anndata_var_names <- function(obj, n = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(if (is.null(n)) character(0) else paste0("Gene_", seq_len(n)))
  }
  ids <- tryCatch(
    as.character(reticulate::py_to_r(obj$var_names$tolist())),
    error = function(e) NULL
  )
  if (is.null(ids) && !is.null(n)) ids <- paste0("Gene_", seq_len(n))
  if (is.null(ids)) character(0) else ids
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
