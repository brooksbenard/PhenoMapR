#' Default Cell Type Color Palette
#'
#' Named vector of hex colors for common cell types. Used in vignettes for
#' consistent coloring across plots. Cell types not in this palette receive
#' fallback colors from \code{get_celltype_palette}.
#'
#' @export
cell_type_colors <- c(
  "Acinar"     = "#66c2a5",
  "Alpha"      = "#a6cee3",
  "B_cell"     = "#1f78b4",
  "Beta"       = "#8c510a",
  "CD4T"       = "#b2df8a",
  "CD8T"       = "#33a02c",
  "Delta"      = "#c51b7d",
  "Dendritic"  = "#fb9a99",
  "Ductal"     = "#e31a1c",
  "Endothelial" = "#8073ac",
  "Fibroblast" = "#fdbf6f",
  "Macrophage" = "#ff7f00",
  "Mast"       = "#cab2d6",
  "NK"         = "#6a3d9a",
  "Plasma"     = "#e5e572",
  "Schwann"    = "#b15928"
)


#' Get Cell Type Color Palette
#'
#' Returns a named vector of colors for the given cell types. Uses
#' \code{cell_type_colors} for known types; assigns fallback colors for others.
#'
#' @param types Character vector of cell type names.
#' @return Named character vector of hex colors.
#' @export
get_celltype_palette <- function(types) {
  types <- unique(as.character(na.omit(types)))
  known <- intersect(types, names(cell_type_colors))
  unknown <- setdiff(types, names(cell_type_colors))
  pal <- character(length(types))
  names(pal) <- types
  if (length(known) > 0) pal[known] <- cell_type_colors[known]
  if (length(unknown) > 0) {
    fallback <- colorRampPalette(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", "#EE8866"))(length(unknown))
    pal[unknown] <- fallback
  }
  pal
}


#' Load RDS File with Fast Parallel Decompression
#'
#' Loads RDS files using \code{fastSave::readRDS.pigz} when available for
#' faster parallel decompression. Installs fastSave from GitHub if not present.
#' Falls back to base \code{readRDS} if fastSave is unavailable or fails.
#'
#' @param file Path to the RDS file.
#' @param install_if_missing If TRUE (default), attempts to install fastSave
#'   via \code{devtools::install_github('barkasn/fastSave')} when not found.
#' @return The R object loaded from the RDS file.
#'
#' @examples
#' \dontrun{
#' obj <- load_rds_fast("path/to/object.rds")
#' }
#' @export
load_rds_fast <- function(file, install_if_missing = TRUE) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  use_fast <- FALSE
  if (requireNamespace("fastSave", quietly = TRUE)) {
    use_fast <- TRUE
  } else if (install_if_missing) {
    tryCatch({
      if (requireNamespace("devtools", quietly = TRUE)) {
        suppressMessages(devtools::install_github("barkasn/fastSave", upgrade = "never", quiet = TRUE))
      }
      if (requireNamespace("fastSave", quietly = TRUE)) {
        use_fast <- TRUE
      }
    }, error = function(e) {
      message("Could not install fastSave: ", conditionMessage(e), "; using base readRDS")
    })
  }
  if (use_fast) {
    tryCatch({
      return(fastSave::readRDS.pigz(file))
    }, error = function(e) {
      message("fastSave::readRDS.pigz failed: ", conditionMessage(e), "; falling back to readRDS")
    })
  }
  readRDS(file)
}


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
#' Simple histogram of score distribution using ggplot2
#'
#' @param scores Numeric vector or data.frame of scores
#' @param score_column If scores is data.frame, which column to plot
#' @param main Plot title
#' @param base_size Base font size for the plot (default 14)
#'
#' @export
plot_score_distribution <- function(scores, score_column = NULL, main = "Score Distribution", base_size = 14) {

  if (is.data.frame(scores)) {
    if (is.null(score_column)) {
      score_column <- colnames(scores)[1]
    }
    scores <- scores[[score_column]]
  }

  df <- data.frame(score = as.numeric(stats::na.omit(scores)))
  med <- stats::median(df$score)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$score)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::geom_vline(xintercept = med, color = "red", linewidth = 1, linetype = 2) +
    ggplot2::labs(title = main, x = "Score", y = "Frequency") +
    ggplot2::theme_minimal(base_size = base_size)

  p
}


#' Validate Package Data
#'
#' Check that all required data objects are present and valid
#'
#' @keywords internal
validate_package_data <- function() {
  
  required_data <- c("precog", "tcga", "pediatric", "ici")
  
  for (data_name in required_data) {
    tryCatch({
      obj <- get_data(data_name)
      if (!is.data.frame(obj) && !is.matrix(obj)) {
        warning(glue::glue("{data_name} is not a data.frame or matrix"))
      }
      if (is.null(rownames(obj))) {
        warning(glue::glue("{data_name} has no rownames (gene names)"))
      }
    }, error = function(e) {
      warning(glue::glue("Could not load {data_name}: {e$message}"))
    })
  }
  
  invisible(TRUE)
}
