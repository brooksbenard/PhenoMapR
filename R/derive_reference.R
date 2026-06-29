#' Derive Reference Z-Scores from Bulk Expression and Phenotype
#'
#' When you have bulk expression (**genes as rows, samples as columns**) and a
#' phenotype (binary, continuous, or survival), this function computes gene-level
#' association z-scores that can be used as a custom reference in \code{\link{PhenoMap}}.
#'
#' Steps: (1) ensure genes x samples format (transpose with message if heuristic
#' suggests the matrix was provided as samples x genes), (2) clean gene names to
#' approved HUGO symbols, (3) check if expression is normalized and normalize
#' if needed, (4) compute phenotype association z-scores per gene (Cox for
#' survival, logistic regression for binary, correlation for continuous).
#'
#' @param bulk_expression Matrix or data.frame. Bulk expression with **genes as
#'   rows** and **samples as columns**. If the matrix appears to be samples x
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
#'   (numeric with >2 unique -> continuous; 2 unique -> binary; or use survival
#'   if \code{survival_time} and \code{survival_event} are provided).
#' @param survival_time Character. Column name in \code{phenotype} for
#'   time-to-event (e.g. overall survival time). Required when
#'   \code{phenotype_type = "survival"}.
#' @param survival_event Character. Column name in \code{phenotype} for event
#'   indicator (0/1 or FALSE/TRUE). Required when
#'   \code{phenotype_type = "survival"}.
#' @param normalize Logical. If \code{TRUE}, run normalization when expression
#'   looks like counts (default \code{TRUE}). Set \code{FALSE} to skip.
#' @param platform One of \code{"auto"}, \code{"rnaseq"}, or \code{"microarray"}.
#'   Microarray inputs follow the PRECOG preprocessing strategy: optional
#'   probe-to-gene mapping, log2 transform when needed, removal of zero-variance
#'   genes and samples/genes with heavy missingness, per-study quantile
#'   normalization, and per-gene unit variance scaling across samples.
#'   RNA-seq inputs use log2(CPM+1) when the matrix looks like raw counts.
#' @param probe_annotation Optional data.frame with probe IDs and gene symbols
#'   (required for probe-level microarray matrices unless row names are already
#'   gene symbols). See \code{\link{collapse_probes_to_genes}}.
#' @param hugo_species Character. Species for HUGO symbol cleaning:
#'   \code{"human"} or \code{"mouse"} (default \code{"human"}).
#' @param binary_positive_reference For \code{phenotype_type} \code{"binary"} only:
#'   which level of the binary factor should correspond to the \strong{positive}
#'   outcome in logistic regression (\code{y = 1}), so that genes with
#'   \strong{positive} z-scores are those whose \strong{higher} expression is
#'   associated with that level. Use \code{"first"} when the first level of
#'   \code{factor(..., levels = c(...))} is your phenotype of interest (e.g.
#'   \code{mutated} vs \code{wt}); use \code{"second"} for the legacy convention
#'   (second level coded as \code{y = 1}). Default \code{"second"} preserves
#'   behaviour from previous versions when levels were implicit (e.g. alphabetical).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return A data.frame with genes as rownames and a single column of
#'   phenotype-association z-scores, suitable for \code{reference} in
#'   \code{\link{PhenoMap}}. When scoring with \code{\link{PhenoMap}}, a
#'   \strong{positive} weighted sum means higher expression of genes with
#'   \strong{positive} reference z is associated with the level you chose via
#'   \code{binary_positive_reference} (for binary phenotypes). Use
#'   \code{PhenoMap(..., reference_sign = -1)} if you need to flip the sign of
#'   the entire reference after the fact.
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
                                       platform = c("auto", "rnaseq", "microarray"),
                                       probe_annotation = NULL,
                                       hugo_species = c("human", "mouse"),
                                       binary_positive_reference = c("second", "first"),
                                       verbose = TRUE) {

  phenotype_type <- match.arg(phenotype_type)
  platform <- match.arg(platform)
  hugo_species <- match.arg(hugo_species)
  binary_positive_reference <- match.arg(binary_positive_reference)

  if (is.data.frame(bulk_expression)) {
    bulk_expression <- as.matrix(bulk_expression)
  }
  if (!is.matrix(bulk_expression)) {
    stop("'bulk_expression' must be a matrix or data.frame")
  }

  # Expected input: genes (rows) x samples (columns).
  #
  # Heuristic: there should be far more genes than samples. If instead the number
  # of columns is orders of magnitude larger than the number of rows, the matrix
  # is likely samples x genes and should be transposed.
  if (nrow(bulk_expression) > 0 && ncol(bulk_expression) > 0 &&
      !is.null(rownames(bulk_expression)) && !is.null(colnames(bulk_expression))) {
    if (ncol(bulk_expression) > 10 * nrow(bulk_expression)) {
      if (isTRUE(verbose)) {
        message(
          "Expression format transposed: expected genes as rows and samples as columns ",
          "(far more genes than samples)."
        )
      }
      bulk_expression <- t(bulk_expression)
    }
  }

  # At this point: genes x samples (or best-effort if dimnames missing).
  gene_names <- rownames(bulk_expression)
  sample_ids <- colnames(bulk_expression)
  if (is.null(gene_names)) gene_names <- paste0("G", seq_len(nrow(bulk_expression)))
  if (is.null(sample_ids)) sample_ids <- paste0("S", seq_len(ncol(bulk_expression)))
  rownames(bulk_expression) <- gene_names
  colnames(bulk_expression) <- sample_ids

  if (platform == "auto") {
    platform <- if (!is.null(probe_annotation) ||
                    .looks_like_microarray_probes(gene_names)) {
      "microarray"
    } else {
      "rnaseq"
    }
    if (verbose && platform == "microarray") {
      message("Inferred microarray platform from probe-like feature IDs.")
    }
  }

  if (platform == "microarray") {
    if (is.null(probe_annotation) && .looks_like_microarray_probes(gene_names)) {
      stop(
        "Microarray expression appears probe-level. Supply 'probe_annotation' ",
        "(GPL-style table with probe and gene symbol columns) or pre-collapse ",
        "probes to gene symbols before calling derive_reference_from_bulk()."
      )
    }
    if (!is.null(probe_annotation)) {
      if (verbose) message("Mapping microarray probes to gene symbols...")
      bulk_expression <- collapse_probes_to_genes(bulk_expression, probe_annotation)
      gene_names <- rownames(bulk_expression)
    }
    if (isTRUE(normalize)) {
      bulk_expression <- preprocess_microarray_expression(
        bulk_expression,
        verbose = verbose
      )
      gene_names <- rownames(bulk_expression)
    }
  }

  # Downstream modeling expects samples x genes
  bulk_expression <- t(bulk_expression)
  rownames(bulk_expression) <- sample_ids
  colnames(bulk_expression) <- gene_names

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
  n_samples_used <- length(common)

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

  # 2) Check normalization and optionally normalize (RNA-seq path only;
  # microarray preprocessing already ran on the genes x samples matrix).
  if (normalize && platform == "rnaseq") {
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
        z_scores <- z_scores_binary(
          bulk_expression,
          y,
          verbose = verbose,
          binary_positive_reference = binary_positive_reference
        )
        score_label <- paste0(phenotype_column, "_z")
      } else {
        z_scores <- z_scores_continuous(bulk_expression, y, verbose = verbose)
        score_label <- paste0(phenotype_column, "_z")
      }
    }
  }

  out <- data.frame(z = z_scores, row.names = names(z_scores))
  names(out) <- score_label
  attr(out, "phenotype_type") <- phenotype_type
  attr(out, "platform") <- platform
  attr(out, "n_samples") <- n_samples_used
  if (verbose) message(glue::glue("Derived reference with {length(z_scores)} genes"))
  return(out)
}


#' Collapse probe-level microarray expression to gene symbols
#'
#' Maps platform probe IDs to gene symbols using a GPL-style annotation table,
#' then averages expression for probes that map to the same symbol. Follows the
#' probe summarization strategy used in the original PRECOG resource.
#'
#' @param expr_probe_by_sample Matrix with probes as rows and samples as columns.
#' @param annot_df Annotation table (e.g. from GEO GPL) with probe and gene columns.
#'
#' @return Matrix with gene symbols as row names and the same sample columns.
#' @export
collapse_probes_to_genes <- function(expr_probe_by_sample, annot_df) {
  if (is.data.frame(expr_probe_by_sample)) {
    expr_probe_by_sample <- as.matrix(expr_probe_by_sample)
  }
  if (!is.matrix(expr_probe_by_sample)) {
    stop("'expr_probe_by_sample' must be a matrix or data.frame")
  }
  annot_df <- as.data.frame(annot_df)
  if (nrow(annot_df) == 0L) stop("'annot_df' has no rows")

  probe_col <- .pick_annotation_column(
    names(annot_df),
    c("ID", "ID_REF", "Probe.Set.ID", "ProbeID", "Probe ID", "Platform ID", "SPOT_ID"),
    regex = "^id(_ref)?$|probe"
  )
  if (is.na(probe_col)) probe_col <- names(annot_df)[1L]

  gene_col <- .pick_annotation_column(
    names(annot_df),
    c("Gene Symbol", "Gene.symbol", "Symbol", "SYMBOL", "GENE_SYMBOL",
      "GeneSymbol", "UniGene symbol", "Gene Symbol(s)"),
    regex = "genesymbol|gene.?symbol|^symbol$"
  )
  if (is.na(gene_col)) {
    stop(
      "Could not infer gene symbol column in probe annotation. Columns: ",
      paste(names(annot_df), collapse = ", ")
    )
  }

  annot <- annot_df[, c(probe_col, gene_col), drop = FALSE]
  names(annot) <- c("probe", "gene_symbol")
  annot$probe <- as.character(annot$probe)
  annot$gene_symbol <- as.character(annot$gene_symbol)
  annot$gene_symbol <- trimws(annot$gene_symbol)
  annot$gene_symbol <- sub("///.*$", "", annot$gene_symbol)
  annot$gene_symbol <- sub(";.*$", "", annot$gene_symbol)
  annot$gene_symbol <- sub(",.*$", "", annot$gene_symbol)

  expr_df <- data.frame(
    probe = rownames(expr_probe_by_sample),
    expr_probe_by_sample,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  merged <- merge(expr_df, annot, by = "probe", all.x = FALSE, all.y = FALSE)
  merged <- merged[!is.na(merged$gene_symbol) & nzchar(merged$gene_symbol), , drop = FALSE]
  if (nrow(merged) == 0L) stop("No probe-to-gene mappings remained after merge.")

  mat <- as.matrix(merged[, setdiff(names(merged), c("probe", "gene_symbol")), drop = FALSE])
  storage.mode(mat) <- "double"
  genes <- merged$gene_symbol
  sample_ids <- colnames(mat)
  gene_levels <- unique(genes)
  out <- matrix(NA_real_, nrow = length(gene_levels), ncol = ncol(mat),
                dimnames = list(gene_levels, sample_ids))
  for (g in gene_levels) {
    idx <- which(genes == g)
    if (length(idx) == 1L) {
      out[g, ] <- mat[idx, , drop = TRUE]
    } else {
      out[g, ] <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  out
}


#' PRECOG-style microarray preprocessing (genes x samples)
#'
#' @param expr_genes_by_sample Numeric matrix, genes as rows, samples as columns.
#' @param min_gene_sd Minimum per-gene standard deviation to retain a gene.
#' @param max_missing_frac Maximum allowed missing fraction per gene or sample.
#' @param verbose Logical.
#'
#' @return Preprocessed matrix (genes x samples).
#' @keywords internal
preprocess_microarray_expression <- function(expr_genes_by_sample,
                                             min_gene_sd = 1e-5,
                                             max_missing_frac = 0.8,
                                             verbose = TRUE) {
  x <- as.matrix(expr_genes_by_sample)
  storage.mode(x) <- "double"
  gene_ids <- rownames(x)
  sample_ids <- colnames(x)

  if (max(x, na.rm = TRUE) > 50 && stats::median(x, na.rm = TRUE) > 1) {
    if (verbose) message("Microarray values do not look log-scaled; applying log2(x + 1)...")
    x <- log2(pmax(x, 0) + 1)
  }

  gene_na_frac <- rowMeans(is.na(x))
  sample_na_frac <- colMeans(is.na(x))
  keep_genes <- gene_na_frac < max_missing_frac
  keep_samples <- sample_na_frac < max_missing_frac
  if (any(!keep_genes) || any(!keep_samples)) {
    if (verbose) {
      message(glue::glue(
        "Removing {sum(!keep_genes)} genes and {sum(!keep_samples)} samples with >={max_missing_frac * 100}% missing values"
      ))
    }
    x <- x[keep_genes, keep_samples, drop = FALSE]
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop("No genes or samples remain after missing-value filtering.")
  }

  gene_sd <- apply(x, 1, stats::sd, na.rm = TRUE)
  keep_var <- is.finite(gene_sd) & gene_sd >= min_gene_sd
  if (any(!keep_var)) {
    if (verbose) {
      message(glue::glue("Removing {sum(!keep_var)} zero-variance genes (SD < {min_gene_sd})"))
    }
    x <- x[keep_var, , drop = FALSE]
  }
  if (nrow(x) == 0L) stop("No variable genes remain after variance filtering.")

  if (verbose) message("Quantile-normalizing microarray samples...")
  x <- .quantile_normalize_genes_by_sample(x)

  if (verbose) message("Scaling each gene to unit variance across samples...")
  x <- t(scale(t(x)))
  x[!is.finite(x)] <- 0
  rownames(x) <- gene_ids[seq_len(nrow(x))]
  colnames(x) <- sample_ids[seq_len(ncol(x))]
  x
}


#' Derive a pan-cohort meta-z reference from multiple bulk studies
#'
#' Computes per-study gene z-scores with \code{\link{derive_reference_from_bulk}}
#' and combines them with a Stouffer weighted meta-analysis (weights =
#' \code{sqrt(n_samples)} per study), analogous to PRECOG / TCGA meta-z signatures.
#'
#' @param studies A named list. Each element must be a list with at least
#'   \code{bulk_expression} and \code{phenotype}. Optional per-study fields are
#'   forwarded to \code{derive_reference_from_bulk()} (e.g. \code{phenotype_type},
#'   \code{platform}, \code{probe_annotation}).
#' @param meta_z_label Character label for the output z-score column
#'   (default \code{"meta_z"}).
#' @param hugo_species Species for HGNC symbol validation (\code{"human"} or
#'   \code{"mouse"}); passed to each per-study call.
#' @param binary_positive_reference For binary phenotypes, which level is the
#'   positive reference (\code{"second"} or \code{"first"}); passed to each
#'   per-study call.
#' @param verbose Logical.
#'
#' @return A data.frame of combined meta-z scores suitable for \code{\link{PhenoMap}}.
#' @export
derive_meta_z_from_bulk_studies <- function(studies,
                                           meta_z_label = "meta_z",
                                           hugo_species = c("human", "mouse"),
                                           binary_positive_reference = c("second", "first"),
                                           verbose = TRUE) {
  hugo_species <- match.arg(hugo_species)
  binary_positive_reference <- match.arg(binary_positive_reference)
  if (!is.list(studies) || length(studies) == 0L) {
    stop("'studies' must be a non-empty list of study specifications")
  }
  if (is.null(names(studies)) || any(!nzchar(names(studies)))) {
    names(studies) <- paste0("study_", seq_along(studies))
  }

  z_by_study <- list()
  n_by_study <- list()
  for (nm in names(studies)) {
    spec <- studies[[nm]]
    if (!is.list(spec) || is.null(spec$bulk_expression) || is.null(spec$phenotype)) {
      stop("Each study must be a list with 'bulk_expression' and 'phenotype'")
    }
    call_args <- list(
      bulk_expression = spec$bulk_expression,
      phenotype = spec$phenotype,
      sample_id_column = if (is.null(spec$sample_id_column)) NULL else spec$sample_id_column,
      phenotype_column = if (is.null(spec$phenotype_column)) NULL else spec$phenotype_column,
      phenotype_type = if (is.null(spec$phenotype_type)) "auto" else spec$phenotype_type,
      survival_time = if (is.null(spec$survival_time)) NULL else spec$survival_time,
      survival_event = if (is.null(spec$survival_event)) NULL else spec$survival_event,
      normalize = if (is.null(spec$normalize)) TRUE else isTRUE(spec$normalize),
      platform = if (is.null(spec$platform)) "auto" else spec$platform,
      probe_annotation = if (is.null(spec$probe_annotation)) NULL else spec$probe_annotation,
      hugo_species = hugo_species,
      binary_positive_reference = binary_positive_reference,
      verbose = isTRUE(verbose)
    )
    ref <- do.call(derive_reference_from_bulk, call_args)
    z_col <- colnames(ref)[1L]
    z_vec <- setNames(as.numeric(ref[[z_col]]), rownames(ref))
    z_by_study[[nm]] <- z_vec
    n_by_study[[nm]] <- attr(ref, "n_samples")
    if (is.null(n_by_study[[nm]]) || is.na(n_by_study[[nm]])) {
      n_by_study[[nm]] <- ncol(spec$bulk_expression)
    }
  }

  all_genes <- unique(unlist(lapply(z_by_study, names), use.names = FALSE))
  meta_z <- setNames(rep(NA_real_, length(all_genes)), all_genes)
  for (g in all_genes) {
    zs <- vapply(names(z_by_study), function(st) {
      z <- z_by_study[[st]][g]
      if (is.na(z)) NA_real_ else z
    }, numeric(1))
    ns <- vapply(names(z_by_study), function(st) {
      if (is.na(zs[st])) 0 else sqrt(as.numeric(n_by_study[[st]]))
    }, numeric(1))
    ok <- is.finite(zs) & ns > 0
    if (!any(ok)) next
    meta_z[g] <- sum(zs[ok] * ns[ok], na.rm = TRUE) / sqrt(sum(ns[ok]^2, na.rm = TRUE))
  }
  meta_z <- meta_z[is.finite(meta_z)]
  if (!length(meta_z)) stop("Meta-z combination produced no finite gene z-scores.")

  out <- data.frame(z = as.numeric(meta_z), row.names = names(meta_z))
  names(out) <- meta_z_label
  attr(out, "phenotype_type") <- "meta_z"
  attr(out, "n_studies") <- length(studies)
  if (verbose) {
    message(glue::glue(
      "Derived meta-z reference from {length(studies)} studies ({length(meta_z)} genes)"
    ))
  }
  out
}


#' @keywords internal
.pick_annotation_column <- function(nm, candidates, regex = NULL) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) return(hit[1L])
  if (!is.null(regex)) {
    idx <- grep(regex, nm, ignore.case = TRUE, value = FALSE)
    if (length(idx)) return(nm[idx[1L]])
  }
  NA_character_
}


#' @keywords internal
.looks_like_microarray_probes <- function(ids) {
  ids <- as.character(ids)
  if (!length(ids)) return(FALSE)
  probeish <- grepl(
    "^(AFFX-|ILMN_|A_23_|probe_|[^A-Za-z0-9])|^[0-9]{6,}$",
    ids,
    ignore.case = TRUE
  )
  mean(probeish) >= 0.5
}


#' @keywords internal
.quantile_normalize_genes_by_sample <- function(expr_genes_by_sample) {
  mat <- as.matrix(expr_genes_by_sample)
  if (ncol(mat) < 2L) return(mat)
  if (requireNamespace("preprocessCore", quietly = TRUE)) {
    return(t(preprocessCore::normalize.quantiles(t(mat))))
  }
  ranked <- apply(mat, 2, rank, ties.method = "average", na.last = "keep")
  ref <- sort(mat[, 1L], na.last = TRUE)
  ref <- ref[!is.na(ref)]
  if (!length(ref)) return(mat)
  out <- mat
  for (j in seq_len(ncol(mat))) {
    r <- ranked[, j]
    ok <- is.finite(r)
    if (!any(ok)) next
    out[ok, j] <- stats::approx(
      x = seq_along(ref) / length(ref),
      y = ref,
      xout = r[ok] / max(r[ok], na.rm = TRUE),
      rule = 2
    )$y
  }
  out
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
z_scores_binary <- function(expr,
                            group,
                            verbose = TRUE,
                            binary_positive_reference = c("second", "first")) {
  binary_positive_reference <- match.arg(binary_positive_reference)
  expr <- as.matrix(expr)
  if (!is.factor(group)) {
    group <- as.factor(group)
  } else {
    group <- droplevels(group)
  }
  levs <- levels(group)
  if (length(levs) != 2) {
    stop("Binary phenotype must have exactly 2 levels; found ", length(levs))
  }
  if (binary_positive_reference == "first") {
    y <- as.integer(group == levs[1])
  } else {
    y <- as.integer(group == levs[2])
  }
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
