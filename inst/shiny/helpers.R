# =============================================================================
# PhenoMapR Shiny app — input parsing, validation, and small UI helpers.
#
# This file is sourced by app.R via source("helpers.R", local = TRUE). All
# functions are kept small and side-effect-free; the app's reactive state lives
# in app.R.
# =============================================================================

# ---- small UI / format helpers ----------------------------------------------

.fmt_int <- function(x) format(as.integer(x), big.mark = ",")
.fmt_num <- function(x, digits = 3) format(round(as.numeric(x), digits), nsmall = digits)

# Cleanly say "x cells" / "x samples" depending on what the user uploaded.
.fmt_n_units <- function(n, unit_singular, unit_plural = NULL) {
  if (is.null(unit_plural)) unit_plural <- paste0(unit_singular, "s")
  paste(.fmt_int(n), if (n == 1L) unit_singular else unit_plural)
}

# Build a small read-only HTML table from a 2+ column data.frame using shiny
# tag helpers. Used inside renderUI() callbacks where we can't (cleanly) embed
# a renderTable output. Safe to insert into a tagList().
tag_table <- function(df, class = "table table-sm table-bordered tag-table") {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    return(shiny::tags$p(shiny::tags$em("(no rows)")))
  }
  shiny::tags$table(
    class = class,
    shiny::tags$thead(
      shiny::tags$tr(
        lapply(colnames(df), shiny::tags$th)
      )
    ),
    shiny::tags$tbody(
      lapply(seq_len(nrow(df)), function(i) {
        shiny::tags$tr(
          lapply(colnames(df), function(cn) {
            shiny::tags$td(as.character(df[i, cn]))
          })
        )
      })
    )
  )
}

# Friendly description block at the top of each tab.
help_block <- function(...) {
  shiny::tags$div(
    class = "help-block",
    shiny::tagList(...)
  )
}

# ---- input file parsers ------------------------------------------------------

# Read an expression upload. Returns a list with:
#   $object       : the loaded R object (matrix, Matrix, Seurat, SCE, data.frame, etc.)
#   $kind         : "matrix" | "seurat" | "sce" | "spatial" | "anndata"
#   $n_genes
#   $n_samples
#   $sample_ids   : character vector of column names
#   $gene_ids     : character vector of row names
#   $notes        : character vector of human-readable status messages
#
# Supported file types (by extension, case-insensitive):
#   .rds                  any object readRDS() supports (Seurat, SCE, matrix, ...)
#   .h5  / .h5ad          10x HDF5 (counts) or AnnData via reticulate
#   .tsv / .csv / .txt    gene × sample matrix; first column = gene IDs
#   .mtx (+ .tsv genes)   not handled here — users should preprocess to .rds
parse_expression_upload <- function(file_path, file_name) {
  notes <- character(0)
  if (!file.exists(file_path)) {
    stop("Uploaded file is missing on disk: ", file_path)
  }
  ext <- tolower(tools::file_ext(file_name))

  obj <- switch(
    ext,
    rds = readRDS(file_path),
    h5 = .read_h5_counts(file_path),
    h5ad = .read_h5ad(file_path),
    tsv = .read_tabular_matrix(file_path, sep = "\t"),
    csv = .read_tabular_matrix(file_path, sep = ","),
    txt = .read_tabular_matrix(file_path, sep = "\t"),
    stop(
      "Unsupported expression file extension: '", ext, "'. ",
      "Accepted: .rds, .h5, .h5ad, .tsv, .csv, .txt"
    )
  )
  if (ext %in% c("h5", "h5ad")) {
    notes <- c(notes, sprintf("Loaded %s from %s.", toupper(ext), file_name))
  } else if (ext == "rds") {
    notes <- c(notes, sprintf("Loaded R object from %s.", file_name))
  } else {
    notes <- c(notes, sprintf("Read tabular matrix from %s.", file_name))
  }
  summary <- summarize_expression_object(obj)
  c(summary, list(object = obj, notes = notes))
}

# Compute kind / dims / sample / gene ids from a loaded expression object.
summarize_expression_object <- function(obj) {
  if (inherits(obj, "Seurat")) {
    if ("Spatial" %in% names(obj@assays)) {
      kind <- "spatial"
    } else {
      kind <- "seurat"
    }
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
  } else if (inherits(obj, "SingleCellExperiment")) {
    kind <- "sce"
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
  } else if (inherits(obj, "SpatialExperiment")) {
    kind <- "spatial"
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
  } else if (is.matrix(obj) || inherits(obj, "Matrix") || is.data.frame(obj)) {
    kind <- "matrix"
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
  } else if (inherits(obj, "python.builtin.object")) {
    kind <- "anndata"
    if (requireNamespace("reticulate", quietly = TRUE)) {
      n_samples <- tryCatch(as.integer(reticulate::py_to_r(obj$n_obs)),
                            error = function(e) NA_integer_)
      n_genes   <- tryCatch(as.integer(reticulate::py_to_r(obj$n_vars)),
                            error = function(e) NA_integer_)
      sample_ids <- tryCatch(PhenoMapR:::.anndata_obs_names(obj, n = n_samples),
                             error = function(e) NULL)
      gene_ids   <- tryCatch(PhenoMapR:::.anndata_var_names(obj, n = n_genes),
                             error = function(e) NULL)
    } else {
      n_samples <- NA_integer_
      n_genes <- NA_integer_
      sample_ids <- NULL
      gene_ids <- NULL
    }
  } else {
    stop(
      "Unsupported expression object class: ", paste(class(obj), collapse = "/"),
      ". Provide a matrix, Matrix, Seurat, SingleCellExperiment, SpatialExperiment, or AnnData."
    )
  }
  list(
    kind = kind,
    n_genes = as.integer(n_genes),
    n_samples = as.integer(n_samples),
    sample_ids = sample_ids,
    gene_ids = gene_ids
  )
}

.read_tabular_matrix <- function(path, sep) {
  # First column treated as gene IDs. Tolerate header presence / absence.
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
  } else {
    dt <- utils::read.table(
      path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  if (ncol(dt) < 2L) {
    stop("Tabular matrix must have at least one ID column plus one sample column.")
  }
  id_col <- dt[[1L]]
  mat <- as.matrix(dt[, -1L, drop = FALSE])
  rownames(mat) <- as.character(id_col)
  storage.mode(mat) <- "double"
  mat
}

.read_h5_counts <- function(path) {
  # Try Seurat's Read10X_h5 for 10X-style HDF5; fall back to raw HDF5 read.
  if (requireNamespace("Seurat", quietly = TRUE)) {
    obj <- tryCatch(Seurat::Read10X_h5(path), error = function(e) NULL)
    if (!is.null(obj)) {
      if (is.list(obj)) obj <- obj[[1L]]
      return(obj)
    }
  }
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    # Heuristic: try matrix at "/matrix" (10X v3)
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    if ("matrix" %in% names(h5)) {
      grp <- h5[["matrix"]]
      data <- grp[["data"]]$read()
      indices <- grp[["indices"]]$read()
      indptr <- grp[["indptr"]]$read()
      shape <- grp[["shape"]]$read()
      features <- grp[["features"]][["name"]]$read()
      barcodes <- grp[["barcodes"]]$read()
      m <- Matrix::sparseMatrix(
        i = indices + 1L,
        p = indptr,
        x = data,
        dims = shape,
        index1 = TRUE
      )
      rownames(m) <- as.character(features)
      colnames(m) <- as.character(barcodes)
      return(m)
    }
  }
  stop(
    "Could not parse HDF5 file. Install Seurat (10X HDF5) or hdf5r, or export to .rds."
  )
}

.read_h5ad <- function(path) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Reading .h5ad requires reticulate + a Python with anndata installed.")
  }
  ad <- tryCatch(reticulate::import("anndata", convert = FALSE),
    error = function(e) {
      stop(
        "Could not import the Python `anndata` module via reticulate.\n",
        "Install with: reticulate::py_install('anndata') ",
        "or point reticulate at a venv where anndata is available."
      )
    }
  )
  ad$read_h5ad(path)
}

# Cell metadata parser (optional). Returns data.frame with cell IDs in 1st col
# unless the file already has a column named after `cell_id_col`.
parse_metadata_upload <- function(file_path, file_name) {
  ext <- tolower(tools::file_ext(file_name))
  df <- switch(
    ext,
    rds = readRDS(file_path),
    tsv = utils::read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    csv = utils::read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unsupported metadata file extension: '", ext, "' (use .tsv/.csv/.rds).")
  )
  if (!is.data.frame(df)) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  df
}

# Custom reference parser (data.frame or named numeric vector → 1-col df).
parse_reference_upload <- function(file_path, file_name) {
  ext <- tolower(tools::file_ext(file_name))
  obj <- switch(
    ext,
    rds = readRDS(file_path),
    tsv = utils::read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1L),
    csv = utils::read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1L),
    txt = utils::read.delim(file_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1L),
    stop("Unsupported reference file extension: '", ext, "' (use .tsv/.csv/.rds).")
  )
  if (is.numeric(obj) && !is.null(names(obj))) {
    obj <- data.frame(z_score = as.numeric(obj), row.names = names(obj))
  }
  if (!is.data.frame(obj)) obj <- as.data.frame(obj)
  if (ncol(obj) < 1L) stop("Reference must have at least one z-score column.")
  if (is.null(rownames(obj))) {
    stop("Reference must have gene symbols as row names.")
  }
  obj
}

# ---- expression -> metadata extraction --------------------------------------

extract_object_metadata <- function(obj) {
  if (inherits(obj, "Seurat")) {
    md <- obj@meta.data
    md$.cell_id <- rownames(md)
    return(md)
  }
  if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SpatialExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(NULL)
    md <- as.data.frame(SummarizedExperiment::colData(obj))
    md$.cell_id <- colnames(obj)
    return(md)
  }
  if (inherits(obj, "python.builtin.object")) {
    md <- tryCatch(PhenoMapR:::.anndata_obs_df(obj), error = function(e) NULL)
    if (is.null(md)) return(NULL)
    md$.cell_id <- rownames(md)
    return(md)
  }
  NULL
}

# ---- safe wrappers around PhenoMapR API -------------------------------------

# Lookup cancer types for a built-in reference, returning a character vector or
# a single error message string.
safe_list_cancer_types <- function(reference) {
  tryCatch(
    {
      ct <- PhenoMapR::list_cancer_types(reference)
      if (is.list(ct)) ct <- ct[[1L]]
      sort(as.character(ct))
    },
    error = function(e) {
      character(0)
    }
  )
}

# Wrap the long-running PhenoMap() call inside a progress notification. The
# heavy lifting itself doesn't expose a per-step hook, so we just show a busy
# spinner and an elapsed-time counter.
run_phenomap_with_progress <- function(expression, reference, cancer_type,
                                       z_score_cutoff, pseudobulk, group_by,
                                       assay, slot, reference_sign,
                                       progress = NULL) {
  if (!is.null(progress)) progress$set(message = "Scoring expression…", value = 0.3)
  t0 <- Sys.time()
  scores <- PhenoMapR::PhenoMap(
    expression = expression,
    reference = reference,
    cancer_type = if (is.null(cancer_type) || !nzchar(cancer_type)) NULL else cancer_type,
    z_score_cutoff = z_score_cutoff,
    pseudobulk = pseudobulk,
    group_by = if (pseudobulk && nzchar(group_by %||% "")) group_by else NULL,
    assay = if (nzchar(assay %||% "")) assay else NULL,
    slot = slot,
    reference_sign = reference_sign,
    verbose = TRUE
  )
  if (!is.null(progress)) progress$set(value = 1.0, message = "Done")
  attr(scores, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  scores
}

# ---- embeddings & expression matrix extraction ------------------------------

# Return the names of available 2D embeddings on a single-cell object.
# Returns a character vector ("umap", "tsne", "pca", ...) or character(0).
list_available_embeddings <- function(obj) {
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    nms <- names(obj@reductions)
    # Surface tissue / spot coordinates from spatial Seurat objects (Visium
    # etc.) as a synthetic "spatial" reduction so the Visualization tab can
    # overlay PhenoMapR scores on the tissue without any extra UI plumbing.
    if (length(methods::slot(obj, "images") %||% list()) > 0L) {
      nms <- c(nms, "spatial")
    }
    return(nms %||% character(0))
  }
  if ((inherits(obj, "SingleCellExperiment") || inherits(obj, "SpatialExperiment")) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    nms <- SingleCellExperiment::reducedDimNames(obj)
    if (inherits(obj, "SpatialExperiment") &&
        requireNamespace("SpatialExperiment", quietly = TRUE)) {
      coords <- tryCatch(SpatialExperiment::spatialCoords(obj),
                         error = function(e) NULL)
      if (!is.null(coords) && NROW(coords) > 0L && NCOL(coords) >= 2L) {
        nms <- unique(c(nms, "spatial"))
      }
    }
    return(nms %||% character(0))
  }
  if (inherits(obj, "python.builtin.object")) {
    keys <- tryCatch(PhenoMapR:::.anndata_obsm_keys(obj),
                     error = function(e) character(0))
    # AnnData spatial coordinates conventionally live in obsm["spatial"];
    # if present, expose them under the same "spatial" label as Seurat /
    # SpatialExperiment so the UI is consistent across object kinds.
    if ("spatial" %in% keys) {
      keys <- unique(c(setdiff(keys, "spatial"), "spatial"))
    }
    return(keys)
  }
  character(0)
}

# Pull a 2D embedding off of obj. Returns a data.frame with columns
# cell_id, dim1, dim2, dim1_name, dim2_name.
extract_embedding <- function(obj, name) {
  emb <- NULL
  is_spatial <- identical(name, "spatial")
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    if (is_spatial) {
      # Pull tissue coordinates from the first image slot. `imagecol`
      # holds the x position; `imagerow` holds the y but is image-space
      # (origin top-left). We negate it later in the plot via
      # `scale_y_reverse()` rather than mutating values here, so the raw
      # coordinates stay in their native frame and any axis text is
      # accurate.
      imgs <- methods::slot(obj, "images") %||% list()
      first <- imgs[[1L]]
      coords <- tryCatch(first@coordinates, error = function(e) NULL)
      if (!is.null(coords) && all(c("imagerow", "imagecol") %in% colnames(coords))) {
        emb <- as.matrix(coords[, c("imagecol", "imagerow"), drop = FALSE])
        rownames(emb) <- rownames(coords)
        colnames(emb) <- c("x", "y")
      }
    } else {
      emb <- Seurat::Embeddings(obj, reduction = name)
    }
  } else if ((inherits(obj, "SingleCellExperiment") || inherits(obj, "SpatialExperiment")) &&
             requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    if (is_spatial && inherits(obj, "SpatialExperiment") &&
        requireNamespace("SpatialExperiment", quietly = TRUE)) {
      emb <- tryCatch(as.matrix(SpatialExperiment::spatialCoords(obj)),
                      error = function(e) NULL)
      if (!is.null(emb) && is.null(rownames(emb))) {
        rownames(emb) <- colnames(obj)
      }
      if (!is.null(emb) && (is.null(colnames(emb)) || !length(colnames(emb)))) {
        colnames(emb) <- c("x", "y")
      }
    } else {
      emb <- SingleCellExperiment::reducedDim(obj, name)
    }
  } else if (inherits(obj, "python.builtin.object")) {
    emb <- tryCatch(PhenoMapR:::.anndata_obsm_array(obj, name),
                    error = function(e) NULL)
  }
  if (is.null(emb)) return(NULL)
  emb <- as.matrix(emb)
  if (ncol(emb) < 2L) return(NULL)
  d1_name <- colnames(emb)[1L]
  d2_name <- colnames(emb)[2L]
  if (is.null(d1_name) || !nzchar(d1_name)) d1_name <- paste0(name, "_1")
  if (is.null(d2_name) || !nzchar(d2_name)) d2_name <- paste0(name, "_2")
  ids <- rownames(emb)
  if (is.null(ids)) {
    # Fall back to colnames of the parent object
    ids <- tryCatch(colnames(obj), error = function(e) NULL)
    if (is.null(ids) || length(ids) != nrow(emb)) {
      ids <- paste0("Cell_", seq_len(nrow(emb)))
    }
  }
  data.frame(
    cell_id = as.character(ids),
    dim1 = as.numeric(emb[, 1L]),
    dim2 = as.numeric(emb[, 2L]),
    dim1_name = d1_name,
    dim2_name = d2_name,
    is_spatial = is_spatial,
    stringsAsFactors = FALSE
  )
}

# Detect 2D-embedding column pairs in a tabular metadata data.frame so the
# UMAP / embedding tab can offer them automatically (saving the user from
# uploading the same file twice when the expression input is a plain matrix
# / data.frame). Returns a named list:
#   list(UMAP = list(dim1 = "UMAP_1", dim2 = "UMAP_2"), tSNE = list(...))
# We pair numeric columns whose names share a common stem and end in either
# 1/2 (Seurat / R convention) or 0/1 (Python / Scanpy convention).
detect_metadata_embeddings <- function(meta_df) {
  if (is.null(meta_df) || !is.data.frame(meta_df) || !ncol(meta_df)) {
    return(list())
  }
  num_cols <- vapply(meta_df, is.numeric, logical(1L))
  cols <- colnames(meta_df)[num_cols]
  if (length(cols) < 2L) return(list())

  rx <- "^(.+?)[ _\\.-]?([012])$"
  m <- regmatches(cols, regexec(rx, cols))
  stems <- vapply(m,
                  function(x) if (length(x) == 3L) x[2L] else NA_character_,
                  character(1L))
  idx <- vapply(m,
                function(x) if (length(x) == 3L) x[3L] else NA_character_,
                character(1L))
  ok <- !is.na(stems) & nzchar(stems) & !is.na(idx)

  pairs <- list()
  for (s in unique(stems[ok])) {
    sub <- which(ok & stems == s)
    sub_idx <- idx[sub]
    sub_cols <- cols[sub]
    cand <- NULL
    if ("1" %in% sub_idx && "2" %in% sub_idx) {
      cand <- list(dim1 = sub_cols[match("1", sub_idx)],
                   dim2 = sub_cols[match("2", sub_idx)])
    } else if ("0" %in% sub_idx && "1" %in% sub_idx) {
      cand <- list(dim1 = sub_cols[match("0", sub_idx)],
                   dim2 = sub_cols[match("1", sub_idx)])
    }
    if (!is.null(cand)) {
      pairs[[s]] <- cand
    }
  }
  pairs
}

# Build a Shiny-app embedding data.frame (cell_id / dim1 / dim2 / dim1_name /
# dim2_name) from a metadata data.frame and a pair of coordinate columns.
# `cell_id_col` defaults to whichever column the Data tab is using.
# ---- patient / sample column detection -------------------------------------
#
# The Dataset-overview tab auto-counts distinct patients and samples without
# asking the user to wire up yet another dropdown. We look for column names
# matching the conventional single-cell / clinical vocabularies:
#   * patient — patient[_id], donor[_id], subject[_id], case[_id], individual[_id]
#   * sample  — sample[_id], library[_id], aliquot, specimen, orig.ident
#
# Both helpers return the first matching column name (NULL when nothing
# obvious is in the metadata).
.detect_meta_column <- function(meta_df, patterns) {
  if (!is.data.frame(meta_df) || !ncol(meta_df)) return(NULL)
  cols <- colnames(meta_df)
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = TRUE, perl = TRUE, value = TRUE)
    if (length(hits)) return(hits[1L])
  }
  NULL
}

detect_patient_column <- function(meta_df) {
  .detect_meta_column(
    meta_df,
    patterns = c(
      "^(patient|patient[._]?id)$",
      "^(donor|donor[._]?id)$",
      "^(subject|subject[._]?id)$",
      "^(case|case[._]?id)$",
      "^(individual|individual[._]?id)$"
    )
  )
}

detect_sample_column <- function(meta_df) {
  .detect_meta_column(
    meta_df,
    patterns = c(
      "^(sample|sample[._]?id)$",
      "^(library|library[._]?id)$",
      "^(aliquot|specimen)$",
      "^orig\\.ident$"
    )
  )
}

# Count distinct non-NA values in a metadata column. Returns NA_integer_ when
# the column is missing or empty so callers can render "—".
count_distinct_meta <- function(meta_df, col) {
  if (is.null(col) || is.null(meta_df) || !nzchar(col) ||
      !(col %in% colnames(meta_df))) {
    return(NA_integer_)
  }
  v <- meta_df[[col]]
  v <- v[!is.na(v) & nzchar(as.character(v))]
  if (!length(v)) return(NA_integer_)
  length(unique(v))
}

extract_embedding_from_metadata <- function(meta_df, dim1_col, dim2_col,
                                            cell_id_col = NULL) {
  if (!is.data.frame(meta_df)) return(NULL)
  if (is.null(dim1_col) || is.null(dim2_col) ||
      !(dim1_col %in% colnames(meta_df)) ||
      !(dim2_col %in% colnames(meta_df))) {
    return(NULL)
  }
  if (is.null(cell_id_col) || !nzchar(cell_id_col) ||
      !(cell_id_col %in% colnames(meta_df))) {
    cell_id_col <- if (".cell_id" %in% colnames(meta_df)) {
      ".cell_id"
    } else {
      colnames(meta_df)[1L]
    }
  }
  ids <- as.character(meta_df[[cell_id_col]])
  if (length(ids) == 0L) return(NULL)
  data.frame(
    cell_id = ids,
    dim1 = as.numeric(meta_df[[dim1_col]]),
    dim2 = as.numeric(meta_df[[dim2_col]]),
    dim1_name = dim1_col,
    dim2_name = dim2_col,
    stringsAsFactors = FALSE
  )
}

# Return a genes × cells expression matrix from a Seurat / SCE / matrix input.
# Falls back to NULL if the object kind isn't recognised.
#
# `gene_subset` (optional character vector): when supplied and the input is
# an AnnData object, only those genes are pulled into R — drastically
# cutting peak memory for multi-GB h5ad files. Ignored for in-memory objects
# (matrix / Seurat / SCE) where row subsetting is already cheap downstream.
extract_expression_matrix <- function(obj, assay = NULL, slot = "data",
                                      gene_subset = NULL) {
  if (is.matrix(obj) || inherits(obj, "Matrix")) return(obj)
  if (is.data.frame(obj)) return(as.matrix(obj))
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    use_assay <- assay
    if (is.null(use_assay) || !nzchar(use_assay)) {
      use_assay <- if ("Spatial" %in% names(obj@assays)) "Spatial" else "RNA"
    }
    layer_name <- if (slot %in% c("data", "counts", "scale.data")) slot else "data"
    m <- tryCatch(
      Seurat::GetAssayData(obj, assay = use_assay, layer = layer_name),
      error = function(e) tryCatch(
        Seurat::GetAssayData(obj, assay = use_assay, slot = slot),
        error = function(e2) NULL
      )
    )
    return(m)
  }
  if ((inherits(obj, "SingleCellExperiment") || inherits(obj, "SpatialExperiment")) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    use_assay <- assay
    if (is.null(use_assay) || !nzchar(use_assay)) {
      anames <- SummarizedExperiment::assayNames(obj)
      use_assay <- if ("logcounts" %in% anames) "logcounts" else anames[1L]
    }
    return(SummarizedExperiment::assay(obj, use_assay))
  }
  if (inherits(obj, "python.builtin.object")) {
    return(tryCatch(
      PhenoMapR:::.anndata_X_to_genes_cells(obj, gene_subset = gene_subset),
      error = function(e) NULL
    ))
  }
  NULL
}

# ---- per-cell-type Wilcoxon brackets ---------------------------------------

# Compute per-cell-type Wilcoxon p-values between two sources and the geometry
# needed for geom_signif-style brackets / annotations. `df` must contain
# `Score`, `Cell_type`, `Source` columns (Source has at most 2 distinct values
# per cell type for the bracket to make sense). Scaled scores are used so the
# y-positioning is in the same coordinate system as a plot of scale(Score).
# ---- ggchicklet2 thin wrappers ---------------------------------------------
#
# https://github.com/brooksbenard/ggchicklet2 — drop-in rounded geoms for
# ggplot2 (`geom_chicklet_bar`, `geom_chicklet_histogram`,
# `geom_chicklet_boxplot`, `geom_chicklet`). When the package is installed
# we use the rounded variants so the app's bar / box / histogram plots
# match the manuscript figures; otherwise we fall back to the standard
# ggplot2 geoms so the app keeps working.
.use_chicklet <- function() {
  requireNamespace("ggchicklet2", quietly = TRUE)
}

# Rounded `geom_col()` — pre-aggregated heights (stat = "identity").
.geom_rounded_col <- function(..., radius = grid::unit(3, "pt"),
                              color = "white") {
  if (.use_chicklet()) {
    ggchicklet2::geom_chicklet_bar(
      stat = "identity", radius = radius, color = color, ...
    )
  } else {
    ggplot2::geom_col(color = color, ...)
  }
}

# Rounded `geom_col(position = "stack", ...)` for stacked bar charts.
.geom_rounded_stack <- function(..., radius = grid::unit(3, "pt"),
                                color = "white") {
  if (.use_chicklet()) {
    ggchicklet2::geom_chicklet_bar(
      stat = "identity", position = "stack",
      radius = radius, color = color, ...
    )
  } else {
    ggplot2::geom_col(position = "stack", color = color, ...)
  }
}

# Rounded `geom_histogram()`.
.geom_rounded_histogram <- function(..., radius = grid::unit(2, "pt"),
                                    color = "white") {
  if (.use_chicklet()) {
    ggchicklet2::geom_chicklet_histogram(
      radius = radius, color = color, ...
    )
  } else {
    ggplot2::geom_histogram(color = color, ...)
  }
}

# Rounded `geom_boxplot()`.
.geom_rounded_boxplot <- function(..., radius = grid::unit(4, "pt")) {
  if (.use_chicklet()) {
    ggchicklet2::geom_chicklet_boxplot(radius = radius, ...)
  } else {
    ggplot2::geom_boxplot(...)
  }
}

# ---- PhenoMapR brand palette ---------------------------------------------
#
# Curated qualitative palette anchored on the app's brand primary (#264653)
# and secondary (#2A9D8F), extended with complementary earth + jewel tones.
# Used as the default discrete fill in the Data-tab composition plots so
# they don't look like default ggplot2 hue swatches. Falls back to
# `colorRampPalette` when more colors than the curated set are needed.
pm_brand_colors <- c(
  "#2A9D8F",  # brand teal
  "#E76F51",  # coral
  "#E9C46A",  # saffron
  "#264653",  # brand charcoal
  "#6A4C93",  # purple
  "#3A86FF",  # azure
  "#8AC926",  # lime
  "#F4A261",  # tangerine
  "#8DA0CB",  # periwinkle
  "#A65A4A"   # terracotta
)

# Pick `n` colors from the brand palette, interpolating with
# `colorRampPalette` when `n` exceeds the curated count so the palette
# still feels coherent for high-cardinality factors.
pm_brand_palette <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) return(character(0))
  if (n <= length(pm_brand_colors)) {
    return(unname(pm_brand_colors[seq_len(n)]))
  }
  grDevices::colorRampPalette(pm_brand_colors)(n)
}

# Drop-in ggplot2 discrete fill scale that picks `n` colors lazily based on
# the number of factor levels seen at training time. Usage:
#   ggplot(...) + scale_fill_phenomapr_d()
#
# Note: the `scale_name` arg of `discrete_scale()` was deprecated in
# ggplot2 3.5.0, so it's omitted here.
scale_fill_phenomapr_d <- function(..., na.value = "#BBBBBB") {
  ggplot2::discrete_scale(
    aesthetics = "fill",
    palette = function(n) pm_brand_palette(n),
    na.value = na.value,
    ...
  )
}

scale_color_phenomapr_d <- function(..., na.value = "#BBBBBB") {
  ggplot2::discrete_scale(
    aesthetics = "colour",
    palette = function(n) pm_brand_palette(n),
    na.value = na.value,
    ...
  )
}

# ---- pairwise Wilcoxon tests between cell types ---------------------------
#
# Used by the "Score by cell type [and source]" plot when the user has *only*
# mapped a cell-type column. Computes a pairwise Wilcoxon test between every
# pair of cell types on the *scaled* score (so the y-positions match what the
# plot draws) and returns the significant ones in a tidy data.frame ready to
# hand to ggsignif::geom_signif(manual = TRUE) — `xmin`, `xmax`, `y_pos`,
# `label` columns plus the cell-type levels as an attribute. Brackets are
# stacked at increasing y so they don't overlap.
celltype_pairwise_pvalues <- function(df,
                                      score_col = "Score",
                                      cell_type_col = "Cell_type",
                                      max_brackets = 12L,
                                      p_threshold = 0.05,
                                      cell_levels = NULL) {
  req_cols <- c(score_col, cell_type_col)
  if (!all(req_cols %in% colnames(df))) return(NULL)
  d <- df[, req_cols]
  colnames(d) <- c("Score", "Cell_type")
  d <- d[is.finite(d$Score) & !is.na(d$Cell_type), ]
  if (nrow(d) < 4L) return(NULL)
  d$Score_scaled <- as.numeric(scale(d$Score))

  # Caller can pass a pre-computed cell-type ordering (e.g. ascending median
  # score) so the bracket xmin/xmax positions line up with the plot axis.
  observed <- unique(as.character(d$Cell_type))
  if (is.null(cell_levels)) {
    cell_levels <- sort(observed)
  } else {
    cell_levels <- intersect(as.character(cell_levels), observed)
    if (!length(cell_levels)) cell_levels <- sort(observed)
  }
  if (length(cell_levels) < 2L) return(NULL)

  y_top <- max(d$Score_scaled, na.rm = TRUE)
  if (!is.finite(y_top)) y_top <- 0

  # All unordered cell-type pairs.
  combos <- utils::combn(cell_levels, 2L, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    a <- pair[[1L]]; b <- pair[[2L]]
    sa <- d$Score_scaled[d$Cell_type == a]
    sb <- d$Score_scaled[d$Cell_type == b]
    if (length(sa) < 2L || length(sb) < 2L) {
      return(NULL)
    }
    p <- tryCatch(
      stats::wilcox.test(sa, sb, alternative = "two.sided",
                         exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    if (is.na(p)) return(NULL)
    label <- if (p < 0.001) "***" else if (p < 0.01) "**" else
             if (p < 0.05)  "*"   else "ns"
    data.frame(
      cell_a = a, cell_b = b, p_val = p, label = label,
      n_a = length(sa), n_b = length(sb),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) return(NULL)

  # Keep only significant pairs to avoid bracket clutter (typical cohorts
  # have many cell types, so all-pairs annotation overwhelms the plot).
  out <- out[out$p_val < p_threshold, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out <- out[order(out$p_val), , drop = FALSE]
  if (nrow(out) > max_brackets) {
    out <- out[seq_len(max_brackets), , drop = FALSE]
  }

  out$xmin <- match(out$cell_a, cell_levels)
  out$xmax <- match(out$cell_b, cell_levels)

  # Stack brackets vertically (most significant nearest the boxes).
  step <- max(0.18, abs(y_top) * 0.07)
  out$y_pos <- y_top * 1.05 + step * (seq_len(nrow(out)) - 1L)

  attr(out, "cell_levels") <- cell_levels
  attr(out, "y_top") <- y_top
  out
}

celltype_source_pvalues <- function(df,
                                    score_col = "Score",
                                    cell_type_col = "Cell_type",
                                    source_col = "Source",
                                    cell_levels = NULL,
                                    significant_only = TRUE,
                                    p_threshold = 0.05) {
  req_cols <- c(score_col, cell_type_col, source_col)
  if (!all(req_cols %in% colnames(df))) return(NULL)
  d <- df[, req_cols]
  colnames(d) <- c("Score", "Cell_type", "Source")
  d <- d[is.finite(d$Score) & !is.na(d$Cell_type) & !is.na(d$Source), ]
  if (nrow(d) < 4L) return(NULL)

  d$Score_scaled <- as.numeric(scale(d$Score))
  y_top <- max(d$Score_scaled, na.rm = TRUE)
  if (!is.finite(y_top)) y_top <- 0

  observed <- unique(as.character(d$Cell_type))
  if (is.null(cell_levels)) {
    cell_levels <- sort(observed)
  } else {
    cell_levels <- intersect(as.character(cell_levels), observed)
    if (!length(cell_levels)) cell_levels <- sort(observed)
  }
  src_levels <- sort(unique(as.character(d$Source)))

  pvals <- lapply(cell_levels, function(ct) {
    sub <- d[d$Cell_type == ct, , drop = FALSE]
    srcs <- unique(sub$Source)
    if (length(srcs) < 2L) {
      return(data.frame(
        Cell_type = ct, p_val = NA_real_, label = "ns",
        n_a = sum(sub$Source == srcs[1L]),
        n_b = NA_integer_,
        source_a = srcs[1L], source_b = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    # If more than 2 sources, use the first two encountered for the contrast.
    a <- srcs[1L]; b <- srcs[2L]
    p <- tryCatch(
      stats::wilcox.test(
        sub$Score_scaled[sub$Source == a],
        sub$Score_scaled[sub$Source == b],
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    lab <- if (is.na(p)) "ns"
           else if (p < 0.001) "***"
           else if (p < 0.01)  "**"
           else if (p < 0.05)  "*"
           else "ns"
    data.frame(
      Cell_type = ct, p_val = p, label = lab,
      n_a = sum(sub$Source == a),
      n_b = sum(sub$Source == b),
      source_a = a, source_b = b,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pvals)
  out <- out[!is.na(out$p_val), , drop = FALSE]
  if (nrow(out) == 0L) return(NULL)
  # Drop non-significant brackets so the plot only annotates what's worth
  # calling out (callers can opt back in with significant_only = FALSE).
  if (isTRUE(significant_only)) {
    out <- out[out$p_val < p_threshold, , drop = FALSE]
    if (nrow(out) == 0L) return(NULL)
  }
  out$x_pos <- match(out$Cell_type, cell_levels)
  out$xmin  <- out$x_pos - 0.2
  out$xmax  <- out$x_pos + 0.2
  out$y_pos <- y_top * 1.05
  attr(out, "cell_levels") <- cell_levels
  attr(out, "src_levels") <- src_levels
  attr(out, "y_top") <- y_top
  out
}

# ---- assemble a unified per-cell data frame --------------------------------
#
# Pull score + phenotype group + cell type + source into a single tibble so
# downstream plots can pivot on whichever combination of columns is present.
# Returns NULL when there are no scores.
build_cell_table <- function(scores,
                              groups = NULL,
                              metadata = NULL,
                              cell_id_col = NULL,
                              cell_type_col = NULL,
                              source_col = NULL) {
  if (is.null(scores)) return(NULL)
  s <- as.data.frame(scores)
  s$cell_id <- rownames(scores)
  score_name <- setdiff(colnames(s), "cell_id")[1L]
  out <- data.frame(
    cell_id = s$cell_id,
    score   = as.numeric(s[[score_name]]),
    stringsAsFactors = FALSE
  )
  attr(out, "score_name") <- score_name

  if (!is.null(groups) && "cell_id" %in% colnames(groups)) {
    grp_col <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1L]
    if (!is.na(grp_col)) {
      g <- data.frame(
        cell_id = as.character(groups$cell_id),
        group   = as.character(groups[[grp_col]]),
        stringsAsFactors = FALSE
      )
      out <- dplyr::left_join(out, g, by = "cell_id")
    }
  }

  if (!is.null(metadata)) {
    md <- as.data.frame(metadata)
    id_col <- cell_id_col
    if (is.null(id_col) || !nzchar(id_col) || id_col == "(none)" ||
        !id_col %in% colnames(md)) {
      id_col <- if (".cell_id" %in% colnames(md)) ".cell_id" else colnames(md)[1L]
    }
    md$cell_id <- as.character(md[[id_col]])

    if (!is.null(cell_type_col) && nzchar(cell_type_col) && cell_type_col != "(none)" &&
        cell_type_col %in% colnames(md)) {
      out$.cell_type <- md[match(out$cell_id, md$cell_id), cell_type_col]
      names(out)[names(out) == ".cell_type"] <- "cell_type"
    }
    if (!is.null(source_col) && nzchar(source_col) && source_col != "(none)" &&
        source_col %in% colnames(md)) {
      out$.source <- md[match(out$cell_id, md$cell_id), source_col]
      names(out)[names(out) == ".source"] <- "source"
    }
  }
  out
}

# ---- small utility: %||% ----------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ============================================================================
# Server-side file picker (single widget for every file input in the app)
# ============================================================================
#
# Every file-input control in the PhenoMapR Shiny app uses a single
# `shinyFiles::shinyFilesButton` styled to look like the standard
# `shiny::fileInput()` "Browse..." control. The picker browses the
# filesystem visible to the R process running the app, so:
#   - locally, that's the user's own machine (same as a fileInput would
#     have been), and
#   - remotely (PhenoMapR::run_app(host = "0.0.0.0")), it's the server's
#     disk, which is what users with large RDS / h5ad files actually want.
#
# Roots are resolved at runtime: home, the current working directory, the
# filesystem root (or available drives on Windows), plus any additional
# paths advertised via PHENOMAPR_SHINY_ROOTS (comma-separated). Override
# the default by passing `roots = ...` to `phenomapr_file_pick()`.
#
# If `shinyFiles` is not installed (it's a Suggests dependency, not an
# Imports), both helpers fall back to the standard browser-side
# `fileInput`, so the app keeps working with reduced functionality.

.has_shinyFiles <- function() {
  isTRUE(requireNamespace("shinyFiles", quietly = TRUE))
}

# Resolve the server-side roots that shinyFiles will expose. Returns a named
# character vector suitable for shinyFiles::shinyFileChoose(roots = ...).
phenomapr_app_server_roots <- function(extra = NULL) {
  roots <- c()
  home <- tryCatch(path.expand("~"), error = function(e) NULL)
  if (!is.null(home) && dir.exists(home)) roots <- c(roots, home = home)
  cwd <- tryCatch(normalizePath(getwd(), mustWork = FALSE),
                  error = function(e) NULL)
  if (!is.null(cwd) && dir.exists(cwd) && (is.null(home) || cwd != home)) {
    roots <- c(roots, working_dir = cwd)
  }
  env <- Sys.getenv("PHENOMAPR_SHINY_ROOTS", unset = "")
  if (nzchar(env)) {
    paths <- trimws(strsplit(env, ",", fixed = TRUE)[[1]])
    paths <- paths[nzchar(paths) & dir.exists(paths)]
    if (length(paths)) {
      names(paths) <- paste0("custom_", seq_along(paths))
      roots <- c(roots, paths)
    }
  }
  if (!is.null(extra) && length(extra)) {
    roots <- c(roots, extra)
  }
  if (.Platform$OS.type == "unix") {
    roots <- c(roots, root = "/")
  } else {
    drives <- tryCatch(
      vapply(LETTERS, function(L) paste0(L, ":/"), character(1)),
      error = function(e) "C:/"
    )
    drives <- drives[vapply(drives, dir.exists, logical(1))]
    if (length(drives)) {
      names(drives) <- tolower(sub(":/$", "_drive", drives))
      roots <- c(roots, drives)
    }
  }
  roots[!duplicated(names(roots))]
}

# UI: a single shinyFiles "Browse..." button + a filename display that mimic
# the look of shiny::fileInput(). Falls back to fileInput when shinyFiles is
# not installed.
phenomapr_file_input <- function(id, label = NULL, accept = NULL,
                                 width = "100%",
                                 button_label = "Browse\u2026",
                                 ...) {
  if (!.has_shinyFiles()) {
    return(shiny::fileInput(id, label, accept = accept, width = width, ...))
  }
  picker_id <- paste0(id, "_server")
  shiny::tags$div(
    class = "form-group shiny-input-container phenomapr-file-input",
    style = if (!is.null(width)) sprintf("width: %s;", width) else NULL,
    if (!is.null(label) && length(label) > 0L && !identical(label, "")) {
      shiny::tags$label(class = "control-label", `for` = picker_id, label)
    },
    shiny::tags$div(
      class = "phenomapr-file-input-row",
      shinyFiles::shinyFilesButton(
        picker_id,
        label = button_label,
        title = "Select a file",
        multiple = FALSE,
        icon = shiny::icon("folder-open"),
        class = "btn btn-default btn-file phenomapr-file-input-btn"
      ),
      shiny::uiOutput(
        paste0(picker_id, "_chosen"),
        inline = TRUE,
        class = "phenomapr-file-input-name"
      )
    )
  )
}

# Server: wire the file picker for `id`. Returns a reactive whose value is
# NULL (nothing picked) or list(datapath, name, source = "server"|"upload").
# Same shape regardless of whether shinyFiles is available.
phenomapr_file_pick <- function(id, input, output, session, roots = NULL,
                                accept = NULL) {
  if (is.null(roots)) roots <- phenomapr_app_server_roots()
  picker_id <- paste0(id, "_server")

  if (.has_shinyFiles()) {
    filetypes <- if (!is.null(accept) && length(accept)) {
      exts <- sub("^\\.", "", tolower(accept))
      exts <- unique(exts[nzchar(exts)])
      if (length(exts)) exts else NULL
    } else NULL

    shinyFiles::shinyFileChoose(
      input,
      picker_id,
      roots = roots,
      filetypes = filetypes,
      session = session
    )

    output[[paste0(picker_id, "_chosen")]] <- shiny::renderUI({
      sel <- input[[picker_id]]
      if (is.null(sel) || identical(sel, integer(0))) {
        return(shiny::tags$span(
          class = "phenomapr-file-input-name-empty",
          "No file selected"
        ))
      }
      path_df <- tryCatch(
        shinyFiles::parseFilePaths(roots, sel),
        error = function(e) NULL
      )
      if (is.null(path_df) || nrow(path_df) == 0L) {
        return(shiny::tags$span(
          class = "phenomapr-file-input-name-empty",
          "No file selected"
        ))
      }
      shiny::tags$span(
        class = "phenomapr-file-input-name-set",
        title = path_df$datapath[1],
        basename(path_df$datapath[1])
      )
    })

    return(shiny::reactive({
      sel <- input[[picker_id]]
      if (is.null(sel) || identical(sel, integer(0))) return(NULL)
      path_df <- tryCatch(
        shinyFiles::parseFilePaths(roots, sel),
        error = function(e) NULL
      )
      if (is.null(path_df) || nrow(path_df) == 0L) return(NULL)
      list(
        datapath = as.character(path_df$datapath[1]),
        name = as.character(path_df$name[1] %||% basename(path_df$datapath[1])),
        source = "server"
      )
    }))
  }

  # Fallback: shinyFiles not installed -> the UI rendered a plain fileInput
  # under `id`, so the picker reactive simply forwards that.
  shiny::reactive({
    up <- input[[id]]
    if (is.null(up) || nrow(up) == 0L) return(NULL)
    list(
      datapath = as.character(up$datapath[1]),
      name = as.character(up$name[1]),
      source = "upload"
    )
  })
}
