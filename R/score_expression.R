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
#' @param slot/layer Data layer for Seurat objects: "data" for normalized, 
#'   "counts" for raw, "scale.data" for scaled (default: "data"). 
#'   Note: In Seurat v5+, this uses the 'layer' parameter; in Seurat v4, 
#'   this uses the 'slot' parameter. The function handles both automatically.
#' @param use_dataset_info Logical. If TRUE, use built-in dataset mapping 
#'   (default: FALSE). Requires 'dataset' parameter.
#' @param dataset Dataset name for automatic cancer type mapping (optional)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A data.frame with samples/cells as rows and score columns. Column 
#'   names follow pattern: "weighted_sum_score_{reference}_{cancer_type}"
#'
#' @examples
#' \dontrun{
#' # Bulk expression matrix
#' scores <- score_expression(
#'   expression = bulk_matrix,
#'   reference = "precog",
#'   cancer_type = "BRCA"
#' )
#'
#' # Single cell Seurat object
#' scores <- score_expression(
#'   expression = seurat_obj,
#'   reference = "tcga",
#'   cancer_type = "LUAD",
#'   assay = "RNA",
#'   slot = "data"
#' )
#'
#' # Spatial with pseudobulk
#' scores <- score_expression(
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
#' scores <- score_expression(
#'   expression = my_data,
#'   reference = custom_ref
#' )
#' }
#'
#' @export
score_expression <- function(expression,
                             reference,
                             cancer_type = NULL,
                             z_score_cutoff = 2,
                             pseudobulk = FALSE,
                             group_by = NULL,
                             assay = NULL,
                             slot = "data",
                             use_dataset_info = FALSE,
                             dataset = NULL,
                             verbose = TRUE) {
  
  # Validate inputs
  if (pseudobulk && is.null(group_by)) {
    stop("'group_by' must be specified when pseudobulk = TRUE")
  }
  
  # Handle reference data
  reference_data <- get_reference_data(reference, cancer_type, use_dataset_info, dataset)
  
  # Convert expression to matrix format
  expr_info <- process_expression_input(
    expression = expression,
    pseudobulk = pseudobulk,
    group_by = group_by,
    assay = assay,
    slot = slot,
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
get_reference_data <- function(reference, cancer_type, use_dataset_info, dataset) {
  
  # If reference is a data.frame, use it directly
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
  
  # Get cancer type
  if (use_dataset_info && !is.null(dataset)) {
    cancer_type <- get_cancer_type_from_dataset(dataset, reference)
  }
  
  if (is.null(cancer_type)) {
    stop("'cancer_type' must be specified or set use_dataset_info = TRUE with 'dataset'")
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


#' Get Cancer Type from Dataset Info
#'
#' @keywords internal
get_cancer_type_from_dataset <- function(dataset, reference) {
  
  datasets_info <- get_data("datasets_info")
  
  if (!dataset %in% datasets_info$dataset_name) {
    stop(glue::glue("Dataset '{dataset}' not found in datasets_info"))
  }
  
  label_col <- switch(reference,
    "precog" = "precog_label",
    "tcga" = "tcga_label",
    "pediatric_precog" = "pediatric_precog_label",
    "ici_precog" = "ici_precog_label",
    stop("Unknown reference for dataset mapping")
  )
  
  cancer_type <- datasets_info %>%
    dplyr::filter(dataset_name == dataset) %>%
    dplyr::pull(!!label_col)
  
  if (length(cancer_type) == 0 || is.na(cancer_type)) {
    stop(glue::glue("No {reference} label found for dataset '{dataset}'"))
  }
  
  return(cancer_type)
}


#' Access Package Data
#'
#' @keywords internal
get_data <- function(name) {
  # This will load data from the package data/ directory
  data_env <- new.env()
  data(list = name, package = "PrognosticScorer", envir = data_env)
  return(data_env[[name]])
}
