#' Detect the format of an expression matrix
#'
#' Inspect a numeric expression matrix (genes x samples) and return a
#' classification of (1) the gene-ID style (HUGO symbols vs Ensembl IDs
#' vs mixed), (2) the likely expression format (raw counts / CPM or
#' TPM / log-normalized / z-scaled), and (3) whether the data looks
#' like single-cell or bulk based on sparsity and the number of
#' samples. The result is intended to power user-facing diagnostics
#' (e.g., the PhenoMapR Shiny app's 1. Data sidebar) and to drive
#' optional cleanup with \code{\link{clean_matrix_input}}.
#'
#' All heuristics are intentionally conservative: when in doubt the
#' detector returns \code{"unknown"} with a low confidence rating so
#' the caller can prompt the user for clarification.
#'
#' @section Detection heuristics:
#' \describe{
#'   \item{Gene-ID kind}{Fraction of rownames matching
#'     \code{^ENSG\\d} (human Ensembl gene IDs) vs the fraction
#'     matching a generic HUGO-symbol pattern (uppercase letters /
#'     digits / hyphens). >50\% Ensembl -> \code{"ensembl"};
#'     >50\% HUGO-like -> \code{"hugo"}; otherwise \code{"mixed"}.
#'     Mouse symbols (mixed-case) are also counted as HUGO-like.}
#'   \item{Format}{Subsamples up to \code{sample_cap} values:
#'     * \code{"z_scaled"} if a substantial fraction (>=10\%) of
#'       values are negative AND the per-column mean is close to 0.
#'     * \code{"raw_counts"} if >=99\% of values are non-negative
#'       integers.
#'     * \code{"cpm_or_tpm"} if column sums cluster near 1e6 and the
#'       matrix max is large (>50) -- TPM and CPM are
#'       indistinguishable without gene-length info.
#'     * \code{"log_normalized"} if values are non-negative, mostly
#'       non-integer, and max is moderate (typically <20).
#'     * \code{"unknown"} otherwise.}
#'   \item{Single-cell vs bulk}{Combines sparsity (fraction of zeros
#'     across the sampled values) and the number of samples
#'     (columns). Sparsity >=0.5 OR (>=200 samples AND sparsity
#'     >=0.3) -> single cell; <200 samples AND sparsity <0.3 ->
#'     bulk; everything else -> unclear.}
#' }
#'
#' @param x Matrix, \pkg{Matrix} sparse matrix, or data.frame with
#'   genes as rows and samples as columns. A data.frame is first
#'   coerced to a numeric matrix.
#' @param sample_cap Maximum number of matrix values to sample for the
#'   numeric heuristics (default 10000). Larger matrices are
#'   subsampled so the detector remains fast on very large objects.
#' @param verbose If \code{TRUE}, prints a one-line message
#'   summarising the detection result.
#'
#' @return A named list with elements: \code{gene_id_kind}
#'   ("hugo"/"ensembl"/"mixed"/"unknown"), \code{n_genes},
#'   \code{n_samples}, \code{n_ensembl}, \code{n_hugo_like},
#'   \code{prop_ensembl}, \code{prop_hugo_like}, \code{n_dup}
#'   (duplicate gene IDs), \code{dup_examples} (up to 5),
#'   \code{format} (one of \code{"raw_counts"}, \code{"cpm_or_tpm"},
#'   \code{"log_normalized"}, \code{"z_scaled"}, \code{"unknown"}),
#'   \code{format_label} (human-readable), \code{format_confidence}
#'   ("high"/"medium"/"low"), \code{sparsity},
#'   \code{sc_or_bulk} ("single_cell"/"bulk"/"unclear"),
#'   \code{sc_or_bulk_label}, \code{stats} (numeric vector of
#'   summary statistics), and \code{recommendations} (character
#'   vector of actionable suggestions, empty when nothing is
#'   recommended).
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(2000, 5), nrow = 100,
#'                  dimnames = list(paste0("G", 1:100), paste0("S", 1:20)))
#' detect_expression_format(counts)
#'
#' @export
detect_expression_format <- function(x, sample_cap = 10000L,
                                     verbose = FALSE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop("'x' must be a matrix, Matrix, or data.frame")
  }

  n_genes <- nrow(x)
  n_samples <- ncol(x)
  gene_ids <- rownames(x)

  out <- list(
    gene_id_kind = "unknown",
    n_genes = n_genes, n_samples = n_samples,
    n_ensembl = 0L, n_hugo_like = 0L,
    prop_ensembl = 0, prop_hugo_like = 0,
    n_dup = 0L, dup_examples = character(0),
    format = "unknown",
    format_label = "Unable to detect (insufficient data)",
    format_confidence = "low",
    sparsity = NA_real_,
    sc_or_bulk = "unclear",
    sc_or_bulk_label = "Cannot determine sc vs bulk",
    stats = c(min = NA_real_, max = NA_real_, mean = NA_real_,
              median = NA_real_, frac_integer = NA_real_,
              frac_negative = NA_real_, frac_zero = NA_real_),
    recommendations = character(0)
  )

  if (n_genes == 0L || n_samples == 0L) return(out)

  # ---- 1) Gene-ID style ---------------------------------------------------
  if (!is.null(gene_ids) && length(gene_ids)) {
    is_ensg <- grepl("^ENSG\\d", gene_ids)
    # HUGO-like: uppercase letters / digits / dashes / dots, length 1-15.
    # Mouse symbols (e.g. "Trp53") are mostly title-case and are also
    # accepted here so we don't mislabel mouse data as "unknown".
    is_hugo <- grepl("^[A-Za-z][A-Za-z0-9.\\-]{0,14}$", gene_ids) & !is_ensg
    n_e <- sum(is_ensg, na.rm = TRUE)
    n_h <- sum(is_hugo, na.rm = TRUE)
    p_e <- n_e / length(gene_ids)
    p_h <- n_h / length(gene_ids)
    out$n_ensembl <- as.integer(n_e)
    out$n_hugo_like <- as.integer(n_h)
    out$prop_ensembl <- p_e
    out$prop_hugo_like <- p_h
    out$gene_id_kind <- if (p_e > 0.5) "ensembl"
      else if (p_h > 0.5) "hugo"
      else if ((p_e + p_h) > 0.5) "mixed"
      else "unknown"

    dup_tbl <- table(gene_ids)
    dup_names <- names(dup_tbl)[dup_tbl > 1L]
    out$n_dup <- length(dup_names)
    if (out$n_dup) {
      out$dup_examples <- utils::head(dup_names, 5L)
    }
  }

  # ---- 2) Numeric heuristics ---------------------------------------------
  vals <- .subsample_matrix_values(x, sample_cap = sample_cap)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(out)

  mn   <- min(vals)
  mx   <- max(vals)
  med  <- stats::median(vals)
  mu   <- mean(vals)
  frac_int <- mean(abs(vals - round(vals)) < 1e-6)
  frac_neg <- mean(vals < 0)
  frac_zero <- mean(vals == 0)

  out$stats <- c(
    min = mn, max = mx, mean = mu, median = med,
    frac_integer = frac_int, frac_negative = frac_neg,
    frac_zero = frac_zero
  )
  out$sparsity <- frac_zero

  # Column-sum heuristic for CPM/TPM. Cheap for both dense and sparse.
  cs_mean <- tryCatch({
    if (inherits(x, "Matrix")) {
      Matrix::mean(Matrix::colSums(x, na.rm = TRUE))
    } else {
      mean(colSums(x, na.rm = TRUE))
    }
  }, error = function(e) NA_real_)
  cs_near_million <- isTRUE(is.finite(cs_mean) &&
                              cs_mean > 9e5 && cs_mean < 1.1e6)

  # ---- 3) Format classification ------------------------------------------
  if (frac_neg >= 0.10 && abs(mu) < 0.5 && mx < 50 && mn > -50) {
    out$format <- "z_scaled"
    out$format_label <- "Z-scaled (scale.data)"
    out$format_confidence <- if (frac_neg >= 0.30) "high" else "medium"
  } else if (frac_int >= 0.99 && frac_neg == 0 && mx >= 1) {
    out$format <- "raw_counts"
    out$format_label <- "Raw counts (integers)"
    out$format_confidence <- "high"
  } else if (cs_near_million && frac_neg == 0 && mx > 50) {
    out$format <- "cpm_or_tpm"
    out$format_label <- "CPM or TPM (per-column sums approx 1e6)"
    out$format_confidence <- "high"
  } else if (frac_int < 0.5 && frac_neg <= 0.001 && mx < 20 && mn >= 0) {
    out$format <- "log_normalized"
    out$format_label <- "Log-normalized (log1p / log2 scale)"
    out$format_confidence <- if (mx < 15) "high" else "medium"
  } else if (frac_int < 0.5 && frac_neg <= 0.001 && mx >= 20 && mx < 1e4) {
    out$format <- "cpm_or_tpm"
    out$format_label <- "CPM/TPM (no log transform applied)"
    out$format_confidence <- "medium"
  } else if (frac_int >= 0.95 && mx > 1e4) {
    out$format <- "raw_counts"
    out$format_label <- "Raw counts (high library depth)"
    out$format_confidence <- "high"
  }

  # ---- 4) Single-cell vs bulk --------------------------------------------
  if (is.finite(frac_zero)) {
    if (frac_zero >= 0.5 ||
        (n_samples >= 200 && frac_zero >= 0.3)) {
      out$sc_or_bulk <- "single_cell"
      out$sc_or_bulk_label <- sprintf(
        "Single-cell-like (%.0f%% zeros, %s samples)",
        frac_zero * 100, .fmt_int_safe(n_samples)
      )
    } else if (n_samples < 200 && frac_zero < 0.3) {
      out$sc_or_bulk <- "bulk"
      out$sc_or_bulk_label <- sprintf(
        "Bulk-like (%.0f%% zeros, %s samples)",
        frac_zero * 100, .fmt_int_safe(n_samples)
      )
    } else {
      out$sc_or_bulk <- "unclear"
      out$sc_or_bulk_label <- sprintf(
        "Unclear (%.0f%% zeros, %s samples)",
        frac_zero * 100, .fmt_int_safe(n_samples)
      )
    }
  }

  # ---- 5) Recommendations ------------------------------------------------
  recs <- character(0)
  if (out$gene_id_kind == "ensembl") {
    recs <- c(recs, paste0(
      "Gene IDs look like Ensembl (", out$n_ensembl, "/", n_genes,
      " ENSG-prefixed). PhenoMapR references use HUGO symbols -- map ",
      "to symbols (e.g. via biomaRt) before scoring."
    ))
  } else if (out$gene_id_kind == "mixed") {
    recs <- c(recs, paste0(
      "Mixed gene-ID styles detected (",
      out$n_ensembl, " ENSG + ", out$n_hugo_like,
      " HUGO-like out of ", n_genes, "). Consider cleaning to HUGO."
    ))
  }
  if (out$n_dup > 0L) {
    ex <- if (length(out$dup_examples)) {
      paste0(" (e.g. ",
             paste(out$dup_examples, collapse = ", "), ")")
    } else ""
    recs <- c(recs, paste0(
      out$n_dup, " duplicate gene ID(s) detected", ex,
      ". Cleaning will collapse them by mean per gene."
    ))
  }
  if (identical(out$format, "raw_counts")) {
    recs <- c(recs, paste0(
      "Raw counts detected. PhenoMapR's weighted-sum scoring expects ",
      "log-normalized expression. Click \"Clean & normalize\" to apply ",
      if (identical(out$sc_or_bulk, "single_cell"))
        "Seurat-style log-normalization (log1p of library-size-scaled counts)."
      else
        "log2(CPM+1) (bulk-style log normalization)."
    ))
  } else if (identical(out$format, "cpm_or_tpm")) {
    recs <- c(recs, paste0(
      "CPM/TPM detected (no log transform). Click \"Clean & normalize\" ",
      "to apply log2(x+1)."
    ))
  } else if (identical(out$format, "z_scaled")) {
    recs <- c(recs, paste0(
      "Z-scaled values detected. PhenoMapR scoring assumes non-negative ",
      "expression; consider re-supplying the un-scaled (data / counts) ",
      "layer instead."
    ))
  } else if (identical(out$format, "unknown")) {
    recs <- c(recs, paste0(
      "Could not confidently identify the expression format. ",
      "Inspect min/max/integer-fraction stats below and pre-process ",
      "before scoring if needed."
    ))
  }
  out$recommendations <- recs

  if (isTRUE(verbose)) {
    message(sprintf(
      "[detect_expression_format] gene IDs: %s; format: %s (%s); %s",
      out$gene_id_kind, out$format_label, out$format_confidence,
      out$sc_or_bulk_label
    ))
  }
  out
}


# Subsample matrix values without densifying a huge sparse matrix.
.subsample_matrix_values <- function(x, sample_cap = 10000L) {
  n_total <- as.numeric(nrow(x)) * as.numeric(ncol(x))
  if (n_total == 0) return(numeric(0))

  if (inherits(x, "sparseMatrix")) {
    # Cap densification by sampling random (row, col) coordinates.
    if (n_total <= sample_cap) {
      return(as.numeric(as.matrix(x)))
    }
    idx <- sample.int(n_total, size = sample_cap)
    rows <- ((idx - 1L) %% nrow(x)) + 1L
    cols <- ((idx - 1L) %/% nrow(x)) + 1L
    vals <- vapply(seq_along(idx),
                   function(k) x[rows[k], cols[k]],
                   numeric(1))
    return(vals)
  }

  if (n_total <= sample_cap) {
    return(as.numeric(x))
  }
  idx <- sample.int(n_total, size = sample_cap)
  as.numeric(x)[idx]
}


.fmt_int_safe <- function(x) {
  if (!is.finite(x)) return("NA")
  format(as.integer(x), big.mark = ",")
}


#' Clean and normalize an expression matrix
#'
#' Optionally clean gene IDs to approved HUGO symbols (collapsing
#' duplicate rows by mean) and apply a log-normalization that matches
#' the chosen analysis mode (single-cell or bulk). The function is a
#' standalone counterpart to the cleanup that
#' \code{\link{derive_reference_from_bulk}} performs internally on the
#' supplied bulk expression matrix; it can be used to pre-process the
#' expression matrix you pass to \code{\link{PhenoMap}} so that gene
#' IDs and value scale match what the reference signatures expect.
#'
#' @section Normalization:
#' Mode is resolved (when \code{mode = "auto"}) from a fresh call to
#' \code{\link{detect_expression_format}}, so callers don't have to
#' duplicate the detection logic. The resulting transformations are:
#' \itemize{
#'   \item single cell, raw counts: per-column library-size scaling to
#'         10,000 followed by \code{log1p}, mirroring Seurat's
#'         \code{LogNormalize}.
#'   \item bulk, raw counts: per-column scaling to 1e6 (CPM) followed
#'         by \code{log2(x + 1)}.
#'   \item raw CPM / TPM (any mode): \code{log2(x + 1)}.
#'   \item already log-normalized or z-scaled: no transformation; the
#'         function still cleans gene IDs if requested.
#' }
#'
#' @param x Matrix / Matrix / data.frame with genes as rows and
#'   samples as columns.
#' @param do_hugo Logical. Clean gene IDs to approved HUGO symbols via
#'   \pkg{HGNChelper} (default \code{TRUE}). A no-op when
#'   \pkg{HGNChelper} is not installed; the function returns the
#'   matrix unchanged on the gene-ID axis and adds a note to
#'   \code{steps}.
#' @param do_collapse_dups Logical. Collapse duplicate gene rows by
#'   mean per sample (default \code{TRUE}). Always honoured even when
#'   \code{do_hugo = FALSE}.
#' @param do_normalize Logical. Apply log-normalization when the
#'   detected format is raw counts or CPM/TPM (default \code{TRUE}).
#' @param mode \code{"auto"} / \code{"single_cell"} / \code{"bulk"}.
#'   When \code{"auto"} the helper calls
#'   \code{\link{detect_expression_format}} and uses its \code{sc_or_bulk}
#'   guess; ties default to bulk.
#' @param hugo_species \code{"human"} or \code{"mouse"} (passed to
#'   \code{HGNChelper::checkGeneSymbols}).
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with elements \code{matrix} (the cleaned matrix),
#'   \code{steps} (character vector of operations performed),
#'   \code{n_collapsed} (number of duplicate gene IDs that were
#'   collapsed by mean), \code{detection} (the
#'   \code{\link{detect_expression_format}} result computed before
#'   any cleanup), and \code{mode} (the resolved mode).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' counts <- matrix(rpois(2000, 5), nrow = 100,
#'                  dimnames = list(paste0("G", 1:100), paste0("S", 1:20)))
#' cleaned <- clean_matrix_input(counts, mode = "bulk")
#' str(cleaned$matrix); cleaned$steps
#' }
#'
#' @export
clean_matrix_input <- function(x,
                               do_hugo = TRUE,
                               do_collapse_dups = TRUE,
                               do_normalize = TRUE,
                               mode = c("auto", "single_cell", "bulk"),
                               hugo_species = c("human", "mouse"),
                               verbose = TRUE) {
  mode <- match.arg(mode)
  hugo_species <- match.arg(hugo_species)

  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop("'x' must be a matrix, Matrix, or data.frame")
  }
  if (is.null(rownames(x))) {
    stop("Matrix must have rownames (gene IDs)")
  }

  steps <- character(0)
  n_collapsed <- 0L
  detection <- detect_expression_format(x, verbose = FALSE)
  resolved_mode <- if (mode == "auto") {
    if (identical(detection$sc_or_bulk, "single_cell")) "single_cell"
    else "bulk"
  } else mode

  # ---- 1) Clean to HUGO ---------------------------------------------------
  if (isTRUE(do_hugo)) {
    if (requireNamespace("HGNChelper", quietly = TRUE)) {
      if (verbose) message("Cleaning gene symbols to approved HUGO IDs...")
      gene_names <- rownames(x)
      checked <- HGNChelper::checkGeneSymbols(
        gene_names, species = hugo_species, unmapped.as.na = FALSE
      )
      new_symbols <- if ("Suggested.Symbol" %in% names(checked)) {
        checked$Suggested.Symbol
      } else {
        checked[[ncol(checked)]]
      }
      na_idx <- is.na(new_symbols) | new_symbols == "" |
        new_symbols == "NA"
      new_symbols[na_idx] <- gene_names[na_idx]
      n_renamed <- sum(new_symbols != gene_names, na.rm = TRUE)
      rownames(x) <- new_symbols
      steps <- c(steps, sprintf(
        "HUGO cleanup: renamed %d gene symbol(s) via HGNChelper",
        n_renamed
      ))
    } else {
      steps <- c(steps, paste0(
        "HUGO cleanup skipped: HGNChelper not installed ",
        "(install.packages('HGNChelper'))"
      ))
    }
  }

  # ---- 2) Collapse duplicate gene rows by mean ----------------------------
  if (isTRUE(do_collapse_dups)) {
    gn <- rownames(x)
    if (anyDuplicated(gn)) {
      n_collapsed <- length(gn) - length(unique(gn))
      if (verbose) {
        message(sprintf(
          "Collapsing %d duplicate gene row(s) by mean expression...",
          n_collapsed
        ))
      }
      # Use rowsum (sum then divide by counts) for speed on dense or
      # sparse matrices.
      if (inherits(x, "Matrix")) {
        x_dense <- as.matrix(x)
      } else {
        x_dense <- x
      }
      sums <- rowsum(x_dense, group = gn, reorder = TRUE,
                     na.rm = TRUE)
      counts_per <- as.integer(table(gn)[rownames(sums)])
      x <- sweep(sums, 1, counts_per, "/")
      steps <- c(steps, sprintf(
        "Collapsed %d duplicate gene row(s) by mean", n_collapsed
      ))
    }
  }

  # ---- 3) Normalize ------------------------------------------------------
  if (isTRUE(do_normalize)) {
    fmt <- detection$format
    if (identical(fmt, "raw_counts")) {
      if (identical(resolved_mode, "single_cell")) {
        if (verbose) {
          message("Single-cell raw counts detected; applying ",
                  "Seurat-style log-normalization (log1p of ",
                  "library-size-scaled counts)...")
        }
        x <- .seurat_lognormalize(x, scale_factor = 1e4)
        steps <- c(steps, paste0(
          "Normalized: Seurat-style log1p(library_size-scaled to 1e4)"
        ))
      } else {
        if (verbose) {
          message("Bulk raw counts detected; applying log2(CPM+1)...")
        }
        x <- .bulk_log2_cpm(x)
        steps <- c(steps, "Normalized: log2(CPM + 1)")
      }
    } else if (identical(fmt, "cpm_or_tpm")) {
      if (verbose) {
        message("CPM/TPM detected (no log transform); applying log2(x+1)...")
      }
      x_dense <- if (inherits(x, "Matrix")) as.matrix(x) else x
      x <- log2(x_dense + 1)
      steps <- c(steps, "Normalized: log2(x + 1) on CPM/TPM")
    } else if (identical(fmt, "log_normalized")) {
      steps <- c(steps,
                 "Normalization skipped: matrix already log-normalized")
    } else if (identical(fmt, "z_scaled")) {
      steps <- c(steps,
                 "Normalization skipped: matrix already z-scaled")
    } else {
      steps <- c(steps,
                 "Normalization skipped: format could not be confidently detected")
    }
  }

  list(matrix = x, steps = steps, n_collapsed = as.integer(n_collapsed),
       detection = detection, mode = resolved_mode)
}


# Seurat-style LogNormalize: per-cell scale to `scale_factor` then log1p.
.seurat_lognormalize <- function(x, scale_factor = 1e4) {
  if (inherits(x, "sparseMatrix")) {
    col_sums <- Matrix::colSums(x, na.rm = TRUE)
    col_sums[col_sums == 0] <- 1
    scaled <- Matrix::t(Matrix::t(x) / col_sums) * scale_factor
    return(log1p(as.matrix(scaled)))
  }
  col_sums <- colSums(x, na.rm = TRUE)
  col_sums[col_sums == 0] <- 1
  scaled <- sweep(x, 2, col_sums, "/") * scale_factor
  log1p(scaled)
}


# Bulk log2(CPM + 1): per-column scale to 1e6 then log2(x + 1).
.bulk_log2_cpm <- function(x) {
  if (inherits(x, "sparseMatrix")) x <- as.matrix(x)
  col_sums <- colSums(x, na.rm = TRUE)
  col_sums[col_sums == 0] <- 1
  cpm <- sweep(x, 2, col_sums, "/") * 1e6
  log2(cpm + 1)
}
