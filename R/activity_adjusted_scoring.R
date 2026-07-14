#' Regress technical covariates and return scaled expression for scoring
#'
#' Uses Seurat normalization, optional cell-cycle scoring, and \code{ScaleData}
#' with \code{vars.to.regress}. Intended for activity-adjusted PhenoMap scoring.
#'
#' @param expression_matrix Genes-by-cells count or normalized matrix.
#' @param signature_genes Gene symbols to scale (typically reference signature genes).
#' @param cell_metadata Optional data.frame with cell covariates; must include a
#'   column matching \code{colnames(expression_matrix)} via \code{cell_id_column}.
#' @param cell_id_column Metadata column with cell barcodes (default \code{"cell_id"}).
#' @param vars_to_regress Character vector of metadata columns to regress out.
#' @param use_counts Logical. If \code{TRUE} (default), treat \code{expression_matrix}
#'   as counts and run \code{NormalizeData}; if \code{FALSE}, use the matrix as
#'   normalized input without re-normalizing.
#' @param verbose Print progress messages.
#'
#' @return A genes-by-cells matrix from the Seurat \code{scale.data} layer.
#' @keywords internal
regress_expression_for_scoring <- function(expression_matrix,
                                         signature_genes,
                                         cell_metadata = NULL,
                                         cell_id_column = "cell_id",
                                         vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
                                         use_counts = TRUE,
                                         verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for activity-adjusted scoring")
  }
  validate_expression_matrix(expression_matrix)
  cells <- colnames(expression_matrix)
  if (is.null(cells) || !length(cells)) {
    stop("expression_matrix must have colnames for activity-adjusted scoring")
  }

  obj <- if (isTRUE(use_counts)) {
    Seurat::CreateSeuratObject(counts = expression_matrix, assay = "RNA")
  } else {
    Seurat::CreateSeuratObject(counts = expression_matrix, assay = "RNA")
  }
  if (isTRUE(use_counts)) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  } else {
    obj <- .set_assay_data_compat(
      obj, assay = "RNA", slot = "data", new.data = expression_matrix
    )
  }

  nfeat <- min(2000L, nrow(obj))
  if (nfeat > 0L) {
    obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE, nfeatures = nfeat)
  }

  if (all(c("S.Score", "G2M.Score") %in% vars_to_regress)) {
    s_genes <- intersect(Seurat::cc.genes.updated.2019$s.genes, rownames(obj))
    g2m_genes <- intersect(Seurat::cc.genes.updated.2019$g2m.genes, rownames(obj))
    if (length(s_genes) >= 5L && length(g2m_genes) >= 5L) {
      obj <- Seurat::CellCycleScoring(
        obj,
        s.features = s_genes,
        g2m.features = g2m_genes,
        set.ident = FALSE
      )
    }
  }

  if (!is.null(cell_metadata)) {
    md <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
    if (!cell_id_column %in% names(md)) {
      stop(glue::glue("cell_metadata is missing column '{cell_id_column}'"))
    }
    rownames(md) <- md[[cell_id_column]]
    for (col in unique(vars_to_regress)) {
      if (col %in% names(md)) {
        obj[[col]] <- md[[col]][match(colnames(obj), md[[cell_id_column]])]
      }
    }
  }

  if (!"nCount_RNA" %in% colnames(obj@meta.data) &&
      "nCount_RNA" %in% vars_to_regress) {
    obj$nCount_RNA <- Matrix::colSums(expression_matrix)
  }

  regressors <- unique(intersect(vars_to_regress, colnames(obj@meta.data)))
  if (!length(regressors) && "nCount_RNA" %in% colnames(obj@meta.data)) {
    regressors <- "nCount_RNA"
  }
  if (!length(regressors)) {
    stop("No regression covariates found for activity-adjusted scoring")
  }

  sig_genes <- intersect(signature_genes, rownames(obj))
  if (length(sig_genes) < 10L) {
    stop("Fewer than 10 signature genes overlap expression after filtering")
  }

  if (verbose) {
    message(glue::glue(
      "Activity-adjusted scoring: regressing {paste(regressors, collapse = ', ')} on {length(sig_genes)} genes"
    ))
  }

  obj <- Seurat::ScaleData(
    obj,
    features = sig_genes,
    vars.to.regress = regressors,
    verbose = FALSE
  )
  reg_mat <- .get_assay_data_compat(obj, assay = "RNA", slot = "scale.data")
  reg_mat
}


#' Shuffle reference z-scores among signature genes (permutation null).
#' @keywords internal
permute_reference_z <- function(meta_z, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  vals <- meta_z
  shuffled <- sample(vals)
  names(shuffled) <- names(meta_z)
  shuffled
}


#' Gene-shuffle permutation p-values for weighted-sum scores.
#'
#' @param expression_matrix Genes-by-cells matrix aligned to \code{meta_z}.
#' @param meta_z Named numeric vector of reference z-scores.
#' @param n_perm Number of label shuffles (default 100).
#' @param seed Random seed.
#' @param pseudobulk Passed to \code{\link{compute_scores}}.
#' @param verbose Print progress.
#'
#' @return List with \code{observed} scores and \code{empirical_p} per cell.
#' @keywords internal
permutation_score_pvalues <- function(expression_matrix,
                                    meta_z,
                                    n_perm = 100L,
                                    seed = 42L,
                                    pseudobulk = FALSE,
                                    verbose = TRUE) {
  n_perm <- as.integer(n_perm)[1L]
  if (n_perm < 1L) stop("'n_perm' must be >= 1")

  common_genes <- intersect(rownames(expression_matrix), names(meta_z))
  common_genes <- sort(common_genes)
  expression_data <- expression_matrix[common_genes, , drop = FALSE]
  meta_z <- meta_z[common_genes]

  observed <- compute_scores(
    expression_data = expression_data,
    prognostic_scores = meta_z,
    pseudobulk = pseudobulk,
    verbose = FALSE
  )
  names(observed) <- colnames(expression_data)

  null_mat <- matrix(NA_real_, nrow = length(observed), ncol = n_perm)
  if (verbose) {
    message(glue::glue("Computing {n_perm} permutation null score(s)..."))
  }
  set.seed(seed)
  for (i in seq_len(n_perm)) {
    perm_z <- permute_reference_z(meta_z, seed = seed + i)
    null_mat[, i] <- compute_scores(
      expression_data = expression_data,
      prognostic_scores = perm_z,
      pseudobulk = pseudobulk,
      verbose = FALSE
    )
  }

  emp_p <- vapply(seq_along(observed), function(j) {
    (sum(null_mat[j, ] >= observed[j], na.rm = TRUE) + 1L) / (n_perm + 1L)
  }, numeric(1))
  names(emp_p) <- names(observed)

  list(observed = observed, empirical_p = emp_p)
}
