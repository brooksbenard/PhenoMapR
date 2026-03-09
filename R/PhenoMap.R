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
#' @param slot Data layer for Seurat objects: "data" for normalized, 
#'   "counts" for raw, "scale.data" for scaled (default: "data"). 
#'   In Seurat v5+ this maps to the layer parameter; in Seurat v4, slot. 
#'   The function handles both automatically.
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A data.frame with samples/cells as rows and score columns. Column 
#'   names follow pattern: \code{weighted_sum_score_\{reference\}_\{cancer_type\}}.
#'   **Directionality**: higher score = worse prognosis (adverse); lower score = 
#'   better prognosis (favorable), matching positive reference z = worse survival.
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
#' # Single cell Seurat object
#' scores <- PhenoMap(
#'   expression = seurat_obj,
#'   reference = "tcga",
#'   cancer_type = "LUAD",
#'   assay = "RNA",
#'   slot = "data"
#' )
#'
#' # Spatial with pseudobulk
#' scores <- PhenoMap(
#'   expression = spatial_seurat,
#'   reference = "ici_precog",
#'   cancer_type = "MELANOMA_Metastatic",
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
                    slot = "data",
                    verbose = TRUE) {

  # Validate inputs
  if (pseudobulk && is.null(group_by)) {
    stop("'group_by' must be specified when pseudobulk = TRUE")
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

  # Calculate scores
  scores <- calculate_weighted_scores(
    expression_matrix = expr_info$matrix,
    reference_data = reference_data,
    z_score_cutoff = z_score_cutoff,
    pseudobulk = pseudobulk,
    score_name = attr(reference_data, "score_name"),
    verbose = verbose
  )

  return(scores)
}


#' Get Reference Prognostic Data
#'
#' @keywords internal
get_reference_data <- function(reference, cancer_type) {

  # If reference is a data.frame (e.g. custom from Cox), use as-is
  if (is.data.frame(reference) || is.matrix(reference)) {
    ref_data <- as.data.frame(reference)
    score_name <- if (!is.null(colnames(ref_data))) {
      colnames(ref_data)[1]
    } else {
      "custom"
    }
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
#' @keywords internal
extract_ici_column <- function(ici_data, cancer_type) {

  if (is.na(cancer_type)) {
    stop("ICI PRECOG label is NA")
  }

  is_metastatic <- grepl("_Metastatic$", cancer_type)
  cancer_abbrev <- sub("_Metastatic$", "", cancer_type)

  pattern <- if (is_metastatic) {
    paste0("^", cancer_abbrev, "_.+_Metastatic_")
  } else {
    paste0("^", cancer_abbrev, "_.+_Primary_")
  }

  matched_cols <- grep(pattern, colnames(ici_data), value = TRUE)

  if (length(matched_cols) == 0) {
    stop(glue::glue("No ICI columns found for pattern: {pattern}"))
  }

  if (length(matched_cols) > 1) {
    warning(glue::glue("Multiple ICI columns found ({length(matched_cols)}), using first: {matched_cols[1]}"))
  }

  return(ici_data[, matched_cols[1], drop = FALSE])
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
