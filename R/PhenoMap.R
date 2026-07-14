#' Score Expression Data with Prognostic Gene Signatures
#'
#' Calculate weighted sum scores for expression data using prognostic z-scores
#' from various reference datasets (PRECOG, TCGA, Pediatric, ICI).
#'
#' @param expression Expression data. Can be:
#'   \itemize{
#'     \item Matrix or data.frame (genes x samples/cells)
#'     \item Seurat object
#'     \item SingleCellExperiment object
#'     \item SpatialExperiment object
#'     \item AnnData object (via reticulate)
#'   }
#' @param reference Reference dataset name ("precog", "tcga", "pediatric_precog", 
#'   "ici_precog") or a custom data.frame with genes as rownames and z-scores 
#'   as values
#' @param cancer_type Cancer type label. Required if using built-in reference 
#'   datasets. Should match column names in reference data. For ICI PRECOG, 
#'   use format "ABBREV" or "ABBREV_Metastatic" (e.g., "MELANOMA", "MELANOMA_Metastatic")
#' @param z_score_cutoff Absolute z-score threshold for filtering genes 
#'   (default: 2)
#' @param pseudobulk Logical. If TRUE, aggregate expression before scoring 
#'   (default: FALSE)
#' @param group_by Column name for pseudobulk grouping. Required if 
#'   pseudobulk = TRUE. For Seurat/SCE objects, should be a metadata column.
#' @param assay Assay name for Seurat/SCE objects (default: "RNA" for sc, 
#'   "Spatial" for spatial)
#' @param layer Seurat matrix to score (\code{"data"} normalized,
#'   \code{"counts"} raw, or \code{"scale.data"}; default \code{"data"}).
#'   Seurat v5 terminology. PhenoMapR maps this to
#'   \code{GetAssayData(layer=...)} or \code{GetAssayData(slot=...)} based on
#'   the installed SeuratObject version.
#' @param slot Alias for \code{layer} (Seurat v4 name for the same matrix).
#'   Prefer \code{layer} in new code. If both are set to different values, an
#'   error is thrown.
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param reference_sign Multiply all reference gene z-scores by this value before
#'   the weighted sum. Use \code{1} (default) or \code{-1} to flip the sign of
#'   resulting PhenoMap scores without changing the reference matrix (e.g. if you
#'   want positive scores to mean the opposite association). For custom references
#'   from \code{\link{derive_reference_from_bulk}}, prefer setting
#'   \code{binary_positive_reference} so the reference itself matches your factor
#'   levels (first level = positive association when \code{"first"}).
#' @param score_mode Scoring mode: \code{"weighted_sum"} (default) uses the
#'   standard dot product of expression and reference z-scores;
#'   \code{"activity_adjusted"} regresses technical covariates per gene (via
#'   Seurat \code{ScaleData}) before computing the weighted sum on the regressed
#'   layer. Requires Seurat for matrix inputs.
#' @param vars_to_regress Metadata columns to regress when
#'   \code{score_mode = "activity_adjusted"}. Default:
#'   \code{c("S.Score", "G2M.Score", "nCount_RNA")}. For matrix input, pass
#'   covariates via \code{sample_metadata} (see below). \code{nCount_RNA} is
#'   computed from column sums when absent.
#' @param sample_metadata Optional data.frame of per-cell/sample covariates for
#'   \code{score_mode = "activity_adjusted"} when \code{expression} is a matrix.
#'   Must include a cell identifier column (default \code{"cell_id"}).
#' @param cell_id_column Column in \code{sample_metadata} matching matrix colnames
#'   (default \code{"cell_id"}).
#' @param permutation_n If a positive integer, compute gene-shuffle permutation
#'   null scores and add \code{empirical_p_*} columns alongside raw score columns.
#' @param permutation_seed Random seed for permutation null (default \code{42}).
#'
#' @return A data.frame with samples/cells as rows and score columns. Column 
#'   names follow pattern: \code{PhenoMapR_\{reference\}_\{cancer_type\}}.
#'   **Directionality (built-in PRECOG/TCGA/ICI references)**: higher score =
#'   worse prognosis (adverse); lower score = better prognosis (favorable),
#'   matching positive reference z = worse survival. **Custom references** from
#'   \code{\link{derive_reference_from_bulk}} follow the sign convention you chose
#'   there (see \code{binary_positive_reference}) times \code{reference_sign}.
#'   For Seurat/SCE objects, scoring may use only reference genes internally for 
#'   memory efficiency. Always add scores to the **same (full) object** using 
#'   \code{add_scores_to_seurat} or \code{add_scores_to_sce} so that all genes 
#'   are retained for downstream analyses (e.g. cell type marker genes).
#'
#' @examples
#' \dontrun{
#' # Bulk expression matrix
#' scores <- PhenoMap(
#'   expression = bulk_matrix,
#'   reference = "precog",
#'   cancer_type = "BRCA"
#' )
#'
#' # Single-cell Seurat object (prefer layer=; Seurat v5)
#' scores <- PhenoMap(
#'   expression = seurat_obj,
#'   reference = "tcga",
#'   cancer_type = "LUAD",
#'   assay = "RNA",
#'   layer = "data"
#' )
#'
#' # Spatial with pseudobulk
#' scores <- PhenoMap(
#'   expression = spatial_seurat,
#'   reference = "ici_precog",
#'   cancer_type = "MELANOMA_Metastatic",
#'   assay = "Spatial",
#'   layer = "data",
#'   pseudobulk = TRUE,
#'   group_by = "sample_id"
#' )
#'
#' # Custom reference data
#' custom_ref <- data.frame(
#'   row.names = c("TP53", "MYC", "EGFR"),
#'   my_signature = c(3.2, -2.5, 2.8)
#' )
#' scores <- PhenoMap(
#'   expression = my_data,
#'   reference = custom_ref
#' )
#' }
#'
#' @export
PhenoMap <- function(expression,
                    reference,
                    cancer_type = NULL,
                    z_score_cutoff = 2,
                    pseudobulk = FALSE,
                    group_by = NULL,
                    assay = NULL,
                    layer = "data",
                    slot = NULL,
                    verbose = TRUE,
                    reference_sign = 1L,
                    score_mode = c("weighted_sum", "activity_adjusted"),
                    vars_to_regress = c("S.Score", "G2M.Score", "nCount_RNA"),
                    sample_metadata = NULL,
                    cell_id_column = "cell_id",
                    permutation_n = NULL,
                    permutation_seed = 42L) {

  score_mode <- match.arg(score_mode)
  slot <- .resolve_seurat_layer_or_slot(layer = layer, slot = slot)
  reference_sign <- as.integer(reference_sign)[1L]
  if (!reference_sign %in% c(-1L, 1L)) {
    stop("'reference_sign' must be 1 or -1")
  }
  if (pseudobulk && is.null(group_by)) {
    stop("'group_by' must be specified when pseudobulk = TRUE")
  }
  if (!is.null(permutation_n)) {
    permutation_n <- as.integer(permutation_n)[1L]
    if (is.na(permutation_n) || permutation_n < 1L) {
      stop("'permutation_n' must be a positive integer")
    }
  }
  if (score_mode == "activity_adjusted" && !requireNamespace("Seurat", quietly = TRUE)) {
    stop("score_mode = 'activity_adjusted' requires the Seurat package")
  }

  # Handle reference data
  reference_data <- get_reference_data(reference, cancer_type)

  # For large objects (Seurat/SCE), extract only reference genes to avoid memory blow-up
  genes_to_extract <- get_reference_genes_for_extraction(reference_data, z_score_cutoff)

  # Convert expression to matrix format
  expr_info <- process_expression_input(
    expression = expression,
    pseudobulk = pseudobulk,
    group_by = group_by,
    assay = assay,
    slot = slot,
    genes_to_extract = genes_to_extract,
    verbose = verbose
  )

  expr_matrix <- expr_info$matrix
  if (score_mode == "activity_adjusted") {
    cell_md <- sample_metadata
    if (is.null(cell_md) && expr_info$input_type == "seurat" && inherits(expression, "Seurat")) {
      cell_md <- expression@meta.data
      if (!cell_id_column %in% names(cell_md)) {
        cell_md[[cell_id_column]] <- rownames(cell_md)
      }
    }
    use_counts <- identical(slot, "counts") || expr_info$input_type == "matrix"
    expr_matrix <- regress_expression_for_scoring(
      expression_matrix = expr_matrix,
      signature_genes = genes_to_extract,
      cell_metadata = cell_md,
      cell_id_column = cell_id_column,
      vars_to_regress = vars_to_regress,
      use_counts = use_counts,
      verbose = verbose
    )
  }

  # Calculate scores
  scores <- calculate_weighted_scores(
    expression_matrix = expr_matrix,
    reference_data = reference_data,
    z_score_cutoff = z_score_cutoff,
    pseudobulk = pseudobulk,
    score_name = attr(reference_data, "score_name"),
    reference_sign = reference_sign,
    verbose = verbose
  )

  if (!is.null(permutation_n)) {
    ref_col <- colnames(reference_data)[1]
    meta_z <- reference_data[[ref_col]]
    names(meta_z) <- rownames(reference_data)
    meta_z <- meta_z[!is.na(meta_z) & abs(meta_z) > z_score_cutoff]
    if (reference_sign == -1L) meta_z <- -meta_z
    perm <- permutation_score_pvalues(
      expression_matrix = expr_matrix,
      meta_z = meta_z,
      n_perm = permutation_n,
      seed = permutation_seed,
      pseudobulk = pseudobulk,
      verbose = verbose
    )
    score_col <- colnames(scores)[1]
    p_col <- sub("^PhenoMapR_", "empirical_p_", score_col)
    scores[[p_col]] <- perm$empirical_p[rownames(scores)]
  }

  return(scores)
}


#' Get Reference Prognostic Data
#'
#' @keywords internal
get_reference_data <- function(reference, cancer_type) {

  # If reference is a data.frame (e.g. custom from Cox), use as-is
  if (is.data.frame(reference) || is.matrix(reference)) {
    ref_data <- as.data.frame(reference)
    # nocov start - as.data.frame() always gives colnames in R
    score_name <- if (!is.null(colnames(ref_data))) {
      colnames(ref_data)[1]
    } else {
      "custom"
    }
    # nocov end
    attr(ref_data, "score_name") <- score_name
    return(ref_data)
  }

  # Otherwise, reference should be a character string
  if (!is.character(reference) || length(reference) != 1) {
    stop("'reference' must be one of: 'precog', 'tcga', 'pediatric_precog', 'ici_precog', or a data.frame")
  }

  reference <- tolower(reference)

  # Load appropriate reference dataset
  ref_obj <- switch(reference,
    "precog" = get_data("precog"),
    "tcga" = get_data("tcga"),
    "pediatric_precog" = get_data("pediatric"),
    "ici_precog" = get_data("ici"),
    stop("Unknown reference: ", reference)
  )

  if (is.null(cancer_type)) {
    stop("'cancer_type' must be specified when using a built-in reference")
  }

  # Extract appropriate column(s)
  if (reference == "ici_precog") {
    ref_data <- extract_ici_column(ref_obj, cancer_type)
  } else {
    if (!cancer_type %in% colnames(ref_obj)) {
      stop(glue::glue("'{cancer_type}' not found in {reference}. Available: {paste(colnames(ref_obj), collapse = ', ')}"))
    }
    ref_data <- as.data.frame(ref_obj[, cancer_type, drop = FALSE])
  }

  score_name <- paste0(reference, "_", colnames(ref_data)[1])
  attr(ref_data, "score_name") <- score_name

  return(ref_data)
}


#' Get genes to extract for memory-efficient scoring (Seurat/SCE)
#' @keywords internal
get_reference_genes_for_extraction <- function(reference_data, z_score_cutoff) {
  ref_df <- as.data.frame(reference_data)
  genes <- character(0)
  for (j in seq_len(ncol(ref_df))) {
    col_name <- colnames(ref_df)[j]
    vec <- ref_df[[col_name]]
    keep <- !is.na(vec) & abs(vec) > z_score_cutoff
    genes <- c(genes, rownames(ref_df)[keep])
  }
  unique(genes)
}


#' Extract ICI PRECOG Column
#'
#' Resolve a user-supplied \code{cancer_type} to a single ICI-PRECOG column.
#' The ICI reference's column names are full cohort labels (e.g.
#' \code{"LUAD_PD-L1_Primary_Naive"}, \code{"KIRC_PD1_Metastatic_Naive"}), and
#' \code{\link{list_cancer_types}("ici_precog")} returns those labels
#' verbatim, so most callers pass an exact column name. We accept that
#' directly. A legacy short-form (e.g. \code{"LUAD"} or
#' \code{"LUAD_Metastatic"}) is also still supported via a regex fallback,
#' which picks the first matching column and warns when there are ties.
#'
#' @keywords internal
extract_ici_column <- function(ici_data, cancer_type) {

  if (is.null(cancer_type) || (length(cancer_type) == 1L && is.na(cancer_type))) {
    stop("ICI PRECOG label is NA")
  }
  if (!is.character(cancer_type) || length(cancer_type) != 1L || !nzchar(cancer_type)) {
    stop("`cancer_type` must be a single non-empty character string for ICI PRECOG.")
  }

  # Fast path: exact match against an existing column name. This is what the
  # Shiny app and `list_cancer_types("ici_precog")` always produce.
  if (cancer_type %in% colnames(ici_data)) {
    return(ici_data[, cancer_type, drop = FALSE])
  }

  # Legacy short-form fallback: cancer abbreviation (optionally with
  # "_Metastatic" suffix), e.g. "LUAD" or "LUAD_Metastatic". Build a regex
  # against the full column names as before.
  is_metastatic <- grepl("_Metastatic$", cancer_type)
  cancer_abbrev <- sub("_Metastatic$", "", cancer_type)
  pattern <- if (is_metastatic) {
    paste0("^", cancer_abbrev, "_.+_Metastatic_")
  } else {
    paste0("^", cancer_abbrev, "_.+_Primary_")
  }
  matched_cols <- grep(pattern, colnames(ici_data), value = TRUE)

  if (length(matched_cols) == 0) {
    stop(glue::glue(
      "No ICI columns matched cancer_type = '{cancer_type}'. ",
      "Pick one of: {paste(colnames(ici_data), collapse = ', ')}"
    ))
  }
  if (length(matched_cols) > 1) {
    warning(glue::glue(
      "Multiple ICI columns matched short-form '{cancer_type}' ",
      "({length(matched_cols)}); using first: {matched_cols[1]}"
    ))
  }
  ici_data[, matched_cols[1], drop = FALSE]
}


#' Access Package Data
#'
#' @keywords internal
get_data <- function(name) {
  # This will load data from the package data/ directory
  data_env <- new.env()
  data(list = name, package = "PhenoMapR", envir = data_env)
  return(data_env[[name]])
}
