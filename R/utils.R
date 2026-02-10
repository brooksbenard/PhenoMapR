#' List Available Cancer Types
#'
#' Show available cancer types for each reference dataset
#'
#' @param reference Reference dataset name ("precog", "tcga", "pediatric_precog", 
#'   "ici_precog"). If NULL, shows all.
#'
#' @return Named list of cancer types for each reference
#'
#' @examples
#' \dontrun{
#' # List all cancer types
#' list_cancer_types()
#'
#' # List only TCGA cancer types
#' list_cancer_types("tcga")
#' }
#'
#' @export
list_cancer_types <- function(reference = NULL) {
  
  refs <- c("precog", "tcga", "pediatric_precog", "ici_precog")
  
  if (!is.null(reference)) {
    if (!reference %in% c(refs, "pediatric", "ici")) {
      stop("reference must be one of: ", paste(refs, collapse = ", "))
    }
    
    # Handle aliases
    if (reference == "pediatric") reference <- "pediatric_precog"
    if (reference == "ici") reference <- "ici_precog"
    
    refs <- reference
  }
  
  result <- list()
  
  for (ref in refs) {
    data_obj <- switch(ref,
      "precog" = get_data("precog"),
      "tcga" = get_data("tcga"),
      "pediatric_precog" = get_data("pediatric"),
      "ici_precog" = get_data("ici")
    )
    
    result[[ref]] <- colnames(data_obj)
  }
  
  if (length(result) == 1) {
    return(result[[1]])
  } else {
    return(result)
  }
}


#' Get Gene Coverage
#'
#' Check which genes from your data are covered in reference datasets
#'
#' @param genes Character vector of gene names
#' @param reference Reference dataset(s) to check. If NULL, checks all.
#'
#' @return Data.frame with gene coverage statistics
#'
#' @examples
#' \dontrun{
#' # Check gene coverage
#' my_genes <- c("TP53", "MYC", "EGFR", "BRCA1")
#' get_gene_coverage(my_genes)
#' }
#'
#' @export
get_gene_coverage <- function(genes, reference = NULL) {
  
  refs <- c("precog", "tcga", "pediatric_precog", "ici_precog")
  
  if (!is.null(reference)) {
    if (reference == "pediatric") reference <- "pediatric_precog"
    if (reference == "ici") reference <- "ici_precog"
    refs <- reference
  }
  
  coverage <- data.frame(
    reference = character(),
    total_genes = integer(),
    covered_genes = integer(),
    coverage_pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (ref in refs) {
    data_obj <- switch(ref,
      "precog" = get_data("precog"),
      "tcga" = get_data("tcga"),
      "pediatric_precog" = get_data("pediatric"),
      "ici_precog" = get_data("ici")
    )
    
    ref_genes <- rownames(data_obj)
    common <- intersect(genes, ref_genes)
    
    coverage <- rbind(coverage, data.frame(
      reference = ref,
      total_genes = length(genes),
      covered_genes = length(common),
      coverage_pct = round(100 * length(common) / length(genes), 2),
      stringsAsFactors = FALSE
    ))
  }
  
  return(coverage)
}


#' Get Top Prognostic Genes
#'
#' Extract top prognostic genes for a specific cancer type
#'
#' @param reference Reference dataset name
#' @param cancer_type Cancer type label
#' @param n Number of top genes to return (default: 50)
#' @param direction "positive", "negative", or "both" (default: "both")
#'
#' @return Data.frame with genes and their z-scores
#'
#' @examples
#' \dontrun{
#' # Get top 100 genes for breast cancer
#' top_genes <- get_top_prognostic_genes("precog", "BRCA", n = 100)
#'
#' # Get top negative prognostic genes
#' bad_genes <- get_top_prognostic_genes("tcga", "LUAD", n = 50, direction = "negative")
#' }
#'
#' @export
get_top_prognostic_genes <- function(reference,
                                      cancer_type,
                                      n = 50,
                                      direction = c("both", "positive", "negative")) {
  
  direction <- match.arg(direction)
  
  # Get reference data
  data_obj <- switch(reference,
    "precog" = get_data("precog"),
    "tcga" = get_data("tcga"),
    "pediatric_precog" = get_data("pediatric"),
    "pediatric" = get_data("pediatric"),
    "ici_precog" = get_data("ici"),
    "ici" = get_data("ici"),
    stop("Unknown reference: ", reference)
  )
  
  if (!cancer_type %in% colnames(data_obj)) {
    stop(glue::glue("'{cancer_type}' not found. Available: {paste(colnames(data_obj), collapse = ', ')}"))
  }
  
  scores <- data.frame(
    gene = rownames(data_obj),
    z_score = data_obj[, cancer_type],
    stringsAsFactors = FALSE
  )
  
  # Remove NAs
  scores <- na.omit(scores)
  
  # Filter by direction
  if (direction == "positive") {
    scores <- scores[scores$z_score > 0, ]
    scores <- scores[order(-scores$z_score), ]
  } else if (direction == "negative") {
    scores <- scores[scores$z_score < 0, ]
    scores <- scores[order(scores$z_score), ]
  } else {
    scores$abs_z <- abs(scores$z_score)
    scores <- scores[order(-scores$abs_z), ]
    scores$abs_z <- NULL
  }
  
  # Return top n
  scores <- head(scores, n)
  rownames(scores) <- NULL
  
  return(scores)
}


#' Plot Score Distribution
#'
#' Simple histogram of score distribution
#'
#' @param scores Numeric vector or data.frame of scores
#' @param score_column If scores is data.frame, which column to plot
#' @param main Plot title
#'
#' @export
plot_score_distribution <- function(scores, score_column = NULL, main = "Score Distribution") {
  
  if (is.data.frame(scores)) {
    if (is.null(score_column)) {
      score_column <- colnames(scores)[1]
    }
    scores <- scores[[score_column]]
  }
  
  hist(scores,
       breaks = 50,
       col = "steelblue",
       border = "white",
       main = main,
       xlab = "Score",
       ylab = "Frequency")
  
  abline(v = median(scores, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = "Median", col = "red", lty = 2, lwd = 2)
}


#' Validate Package Data
#'
#' Check that all required data objects are present and valid
#'
#' @keywords internal
validate_package_data <- function() {
  
  required_data <- c("precog", "tcga", "pediatric", "ici", "datasets_info")
  
  for (data_name in required_data) {
    tryCatch({
      obj <- get_data(data_name)
      
      if (data_name != "datasets_info") {
        if (!is.data.frame(obj) && !is.matrix(obj)) {
          warning(glue::glue("{data_name} is not a data.frame or matrix"))
        }
        if (is.null(rownames(obj))) {
          warning(glue::glue("{data_name} has no rownames (gene names)"))
        }
      }
      
    }, error = function(e) {
      warning(glue::glue("Could not load {data_name}: {e$message}"))
    })
  }
  
  invisible(TRUE)
}
