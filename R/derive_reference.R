#' Derive Reference Z-Scores from Bulk Expression and Phenotype
#'
#' When you have bulk expression (**genes as rows, samples as columns**) and a
#' phenotype (binary, continuous, or survival), this function computes gene-level
#' association z-scores that can be used as a custom reference in \code{\link{PhenoMap}}.
#'
#' Steps: (1) ensure genes × samples format (transpose with message if heuristic
#' suggests the matrix was provided as samples × genes), (2) clean gene names to
#' approved HUGO symbols, (3) check if expression is normalized and normalize
#' if needed, (4) compute phenotype association z-scores per gene (Cox for
#' survival, logistic regression for binary, correlation for continuous).
#'
#' @param bulk_expression Matrix or data.frame. Bulk expression with **genes as
#'   rows** and **samples as columns**. If the matrix appears to be samples ×
#'   genes (e.g. fewer rows than columns), the function will transpose and
#'   message the user.
#' @param phenotype Data.frame with sample identifiers and phenotype column(s).
#'   Must align with \code{bulk_expression} by sample ID (see \code{sample_id_column});
#'   sample IDs must match the **column names** of \code{bulk_expression}.
#' @param sample_id_column Character. Column name in \code{phenotype} that
#'   matches the column names of \code{bulk_expression} (sample IDs).
#'   If \code{NULL}, the first column of \code{phenotype} is used.
#' @param phenotype_column Character. Column name in \code{phenotype} for the
#'   outcome. For \code{phenotype_type = "survival"} this is ignored; use
#'   \code{survival_time} and \code{survival_event} instead.
#' @param phenotype_type One of \code{"auto"}, \code{"survival"}, \code{"binary"},
#'   \code{"continuous"}. If \code{"auto"}, inferred from the phenotype column
#'   (numeric with >2 unique → continuous; 2 unique → binary; or use survival
#'   if \code{survival_time} and \code{survival_event} are provided).
#' @param survival_time Character. Column name in \code{phenotype} for
#'   time-to-event (e.g. overall survival time). Required when
#'   \code{phenotype_type = "survival"}.
#' @param survival_event Character. Column name in \code{phenotype} for event
#'   indicator (0/1 or FALSE/TRUE). Required when
#'   \code{phenotype_type = "survival"}.
#' @param normalize Logical. If \code{TRUE}, run normalization when expression
#'   looks like counts (default \code{TRUE}). Set \code{FALSE} to skip.
#' @param hugo_species Character. Species for HUGO symbol cleaning:
#'   \code{"human"} or \code{"mouse"} (default \code{"human"}).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return A data.frame with genes as rownames and a single column of
#'   phenotype-association z-scores, suitable for \code{reference} in
#'   \code{\link{PhenoMap}}.
#'
#' @examples
#' \dontrun{
#' # Simulated bulk: genes (rows) x samples (columns)
#' set.seed(1)
#' bulk <- matrix(rnorm(100 * 50), 100, 50,
#'   dimnames = list(paste0("G", 1:100), paste0("S", 1:50)))
#' pheno <- data.frame(
#'   sample_id = paste0("S", 1:50),
#'   response = sample(c("R", "NR"), 50, replace = TRUE))
#'
#' ref <- derive_reference_from_bulk(
#'   bulk_expression = bulk,
#'   phenotype = pheno,
#'   sample_id_column = "sample_id",
#'   phenotype_column = "response",
#'   phenotype_type = "binary")
#'
#' # Use in scoring
#' scores <- PhenoMap(expression = my_single_cell_data, reference = ref)
#' }
#'
#' @export
derive_reference_from_bulk <- function(bulk_expression,
                                       phenotype,
                                       sample_id_column = NULL,
                                       phenotype_column = NULL,
                                       phenotype_type = c("auto", "survival", "binary", "continuous"),
                                       survival_time = NULL,
                                       survival_event = NULL,
                                       normalize = TRUE,
                                       hugo_species = c("human", "mouse"),
                                       verbose = TRUE) {

  phenotype_type <- match.arg(phenotype_type)
  hugo_species <- match.arg(hugo_species)

  if (is.data.frame(bulk_expression)) {
    bulk_expression <- as.matrix(bulk_expression)
  }
  if (!is.matrix(bulk_expression)) {
    stop("'bulk_expression' must be a matrix or data.frame")
  }

  # Expected input: genes (rows) x samples (columns). Heuristic: if there are not
  # far more genes than samples, assume user provided samples x genes and transpose.
  n_row <- nrow(bulk_expression)
  n_col <- ncol(bulk_expression)
  if (n_row > 0 && n_col > 0 && !is.null(rownames(bulk_expression)) && !is.null(colnames(bulk_expression))) {
    if (n_row < 10 * n_col) {
      if (verbose) {
        message(
          "Expression format transposed: expected genes as rows and samples as columns ",
          "(many more rows than columns). Matrix was transposed."
        )
      }
      bulk_expression <- t(bulk_expression)
      # After transpose: rows = samples, cols = genes
      sample_ids <- rownames(bulk_expression)
      gene_names <- colnames(bulk_expression)
      if (is.null(sample_ids)) sample_ids <- paste0("S", seq_len(nrow(bulk_expression)))
      if (is.null(gene_names)) gene_names <- paste0("G", seq_len(ncol(bulk_expression)))
      rownames(bulk_expression) <- sample_ids
      colnames(bulk_expression) <- gene_names
    } else {
      # Input is genes x samples. Transpose to samples x genes for downstream.
      gene_names <- rownames(bulk_expression)
      sample_ids <- colnames(bulk_expression)
      if (is.null(gene_names)) gene_names <- paste0("G", seq_len(nrow(bulk_expression)))
      if (is.null(sample_ids)) sample_ids <- paste0("S", seq_len(ncol(bulk_expression)))
      bulk_expression <- t(bulk_expression)
      rownames(bulk_expression) <- sample_ids
      colnames(bulk_expression) <- gene_names
    }
  } else {
    # No dimnames: assume genes x samples, transpose to samples x genes
    gene_names <- rownames(bulk_expression)
    sample_ids <- colnames(bulk_expression)
    if (is.null(gene_names)) gene_names <- paste0("G", seq_len(nrow(bulk_expression)))
    if (is.null(sample_ids)) sample_ids <- paste0("S", seq_len(ncol(bulk_expression)))
    bulk_expression <- t(bulk_expression)
    rownames(bulk_expression) <- sample_ids
    colnames(bulk_expression) <- gene_names
  }
  # Now: rows = samples, cols = genes

  phenotype <- as.data.frame(phenotype)
  if (nrow(phenotype) == 0) stop("'phenotype' has no rows")

  id_col <- if (is.null(sample_id_column)) names(phenotype)[1] else sample_id_column
  if (!id_col %in% names(phenotype)) {
    stop("Sample ID column '", id_col, "' not found in phenotype")
  }
  pheno_ids <- as.character(phenotype[[id_col]])

  common <- intersect(sample_ids, pheno_ids)
  if (length(common) == 0) {
    stop("No overlapping sample IDs between bulk_expression and phenotype")
  }
  if (verbose) {
    message(glue::glue("Using {length(common)} samples common between expression and phenotype"))
  }

  bulk_expression <- bulk_expression[common, , drop = FALSE]
  phenotype <- phenotype[match(common, pheno_ids), , drop = FALSE]

  # 1) Clean gene names to HUGO
  if (requireNamespace("HGNChelper", quietly = TRUE)) {
    if (verbose) message("Cleaning gene symbols to approved HUGO IDs...")
    gene_names <- colnames(bulk_expression)
    checked <- HGNChelper::checkGeneSymbols(gene_names, species = hugo_species, unmapped.as.na = FALSE)
    new_symbols <- if ("Suggested.Symbol" %in% names(checked)) {
      checked$Suggested.Symbol
    } else {
      checked[[ncol(checked)]]
    }
    na_idx <- is.na(new_symbols) | new_symbols == "" | new_symbols == "NA"
    new_symbols[na_idx] <- gene_names[na_idx]
    colnames(bulk_expression) <- new_symbols
    # Collapse duplicates by mean
    # nocov start - requires genes that collapse to same HUGO symbol
    ugenes <- unique(new_symbols)
    if (length(ugenes) < length(new_symbols)) {
      bulk_expression <- apply(bulk_expression, 1, function(x) {
        tapply(x, new_symbols, mean, na.rm = TRUE)
      })
      bulk_expression <- t(bulk_expression)
      if (verbose) message(glue::glue("Collapsed to {ncol(bulk_expression)} unique genes"))
    }
    # nocov end
  } else if (verbose) {
    message("Package 'HGNChelper' not installed; skipping HUGO symbol cleaning. Install with: install.packages('HGNChelper')")
  }

  # 2) Check normalization and optionally normalize
  if (normalize) {
    x <- as.numeric(bulk_expression)
    looks_like_counts <- (all(x == floor(x)) && max(x, na.rm = TRUE) > 100) ||
      (max(x, na.rm = TRUE) > 1e4)
    if (looks_like_counts) {
      if (verbose) message("Expression looks like counts; applying log2(CPM+1)...")
      bulk_expression <- log2(edgeR_cpm_safe(bulk_expression) + 1)
    }
  }

  # 3) Resolve phenotype type and compute gene z-scores
  if (phenotype_type == "survival") {
    if (is.null(survival_time) || is.null(survival_event)) {
      # nocov start - error path
      stop("For phenotype_type = 'survival', provide 'survival_time' and 'survival_event' column names")
      # nocov end
    }
    time_vec <- phenotype[[survival_time]]
    event_vec <- phenotype[[survival_event]]
    if (is.logical(event_vec)) event_vec <- as.integer(event_vec)
    keep <- !is.na(time_vec) & !is.na(event_vec)
    bulk_expression <- bulk_expression[keep, , drop = FALSE]
    time_vec <- time_vec[keep]
    event_vec <- event_vec[keep]
    z_scores <- z_scores_survival(bulk_expression, time_vec, event_vec, verbose = verbose)
    score_label <- "survival_z"
  } else {
    if (is.null(phenotype_column)) phenotype_column <- names(phenotype)[2]
    if (!phenotype_column %in% names(phenotype)) {
      stop("Phenotype column '", phenotype_column, "' not found")
    }
    y <- phenotype[[phenotype_column]]
    if (phenotype_type == "auto") {
      u <- unique(na.omit(y))
      # nocov start - auto infer survival (needs phenotype with time/event columns)
      if (!is.null(survival_time) && survival_time %in% names(phenotype) &&
          !is.null(survival_event) && survival_event %in% names(phenotype)) {
        phenotype_type <- "survival"
        time_vec <- phenotype[[survival_time]]
        event_vec <- phenotype[[survival_event]]
        if (is.logical(event_vec)) event_vec <- as.integer(event_vec)
        keep <- !is.na(time_vec) & !is.na(event_vec)
        bulk_expression <- bulk_expression[keep, , drop = FALSE]
        time_vec <- time_vec[keep]
        event_vec <- event_vec[keep]
        z_scores <- z_scores_survival(bulk_expression, time_vec, event_vec, verbose = verbose)
        score_label <- "survival_z"
      } else if (length(u) == 2) {
      # nocov end
        phenotype_type <- "binary"
      } else if (is.numeric(y) && length(u) > 2) {
        phenotype_type <- "continuous"
      } else {
        phenotype_type <- "binary"
      }
    }
    if (phenotype_type != "survival") {
      keep <- !is.na(y)
      bulk_expression <- bulk_expression[keep, , drop = FALSE]
      y <- y[keep]
      if (phenotype_type == "binary") {
        z_scores <- z_scores_binary(bulk_expression, y, verbose = verbose)
        score_label <- paste0(phenotype_column, "_z")
      } else {
        z_scores <- z_scores_continuous(bulk_expression, y, verbose = verbose)
        score_label <- paste0(phenotype_column, "_z")
      }
    }
  }

  out <- data.frame(z = z_scores, row.names = names(z_scores))
  names(out) <- score_label
  if (verbose) message(glue::glue("Derived reference with {length(z_scores)} genes"))
  return(out)
}


#' CPM without requiring edgeR
#'
#' @keywords internal
# nocov start - only used when normalize=TRUE and counts detected in derive_reference
edgeR_cpm_safe <- function(x) {
  if (nrow(x) == 0) return(x)
  lib_size <- rowSums(x, na.rm = TRUE)
  lib_size[lib_size == 0] <- 1
  sweep(x, 1, lib_size, "/") * 1e6
}
# nocov end


#' Z-scores from Cox PH per gene
#'
#' @keywords internal
z_scores_survival <- function(expr, time, event, verbose = TRUE) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for phenotype_type = 'survival'. Install with: install.packages('survival')")
  }
  expr <- as.matrix(expr)
  n_genes <- ncol(expr)
  z <- setNames(rep(NA_real_, n_genes), colnames(expr))
  for (j in seq_len(n_genes)) {
    fit <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ expr[, j]),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      coef_summary <- summary(fit)$coefficients
      if (!is.null(coef_summary)) {
        z[j] <- coef_summary[1, "z"]
      }
    }
  }
  na_before <- sum(is.na(z))
  if (verbose && na_before > 0) {
    # nocov start - verbose
    message(glue::glue("Cox PH: {na_before} genes had NA z-scores (convergence or low variation)"))
    # nocov end
  }
  z
}


#' Z-scores from logistic regression per gene (binary outcome)
#'
#' @keywords internal
z_scores_binary <- function(expr, group, verbose = TRUE) {
  expr <- as.matrix(expr)
  group <- as.factor(group)
  levs <- levels(group)
  if (length(levs) != 2) {
    stop("Binary phenotype must have exactly 2 levels; found ", length(levs))
  }
  y <- as.integer(group == levs[2])
  n_genes <- ncol(expr)
  z <- setNames(rep(NA_real_, n_genes), colnames(expr))
  for (j in seq_len(n_genes)) {
    fit <- tryCatch(
      stats::glm(y ~ expr[, j], family = stats::binomial),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(fit) && fit$converged) {
      coef_summary <- summary(fit)$coefficients
      if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
        z[j] <- coef_summary[2, "z value"]
      }
    }
  }
  if (verbose && sum(is.na(z)) > 0) {
    # nocov start - verbose
    message(glue::glue("Logistic regression: {sum(is.na(z))} genes had NA z-scores (convergence or separation)"))
    # nocov end
  }
  z
}


#' Z-scores from correlation per gene (Fisher z or test statistic)
#'
#' @keywords internal
z_scores_continuous <- function(expr, y, verbose = TRUE) {
  expr <- as.matrix(expr)
  y <- as.numeric(y)
  if (length(y) != nrow(expr)) stop("Length of y must match number of samples")
  r <- apply(expr, 2, function(x) {
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 4) return(NA)
    cor(x[ok], y[ok], use = "pairwise.complete.obs")
  })
  n <- nrow(expr)
  # Fisher z-transform: z = 0.5 * log((1+r)/(1-r)) * sqrt(n-3)
  r[r >= 1] <- 1 - 1e-6
  r[r <= -1] <- -1 + 1e-6
  z <- 0.5 * log((1 + r) / (1 - r)) * sqrt(n - 3)
  setNames(as.numeric(z), colnames(expr))
}
