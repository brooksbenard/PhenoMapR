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
#
# Also returns:
#   - assays_avail   : character vector of assay names (Seurat/SCE/SpatialExp);
#                      length 0 for plain matrices / AnnData (AnnData uses
#                      `X` + named layers, surfaced via `layers_avail`).
#   - default_assay  : single string with the natural default assay
#                      (Seurat::DefaultAssay or assayNames(.)[1]); NA otherwise.
#   - layers_avail   : character vector of available layer names within the
#                      default assay (Seurat: "counts" / "data" / "scale.data";
#                      SCE: same as assays_avail; AnnData: c("X", layers...)).
# These are used by the Score tab to (a) show a "Detected: ..." status block,
# (b) pre-fill the assay text input with the detected default, and (c) pick
# a sensible default for the slot/layer radio (prefer "data" if available,
# otherwise "counts").
summarize_expression_object <- function(obj) {
  assays_avail   <- character(0)
  default_assay  <- NA_character_
  layers_avail   <- character(0)

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
    assays_avail <- tryCatch(names(obj@assays),
                             error = function(e) character(0))
    default_assay <- tryCatch(
      if (requireNamespace("Seurat", quietly = TRUE))
        Seurat::DefaultAssay(obj)
      else if (length(assays_avail))
        assays_avail[1L]
      else NA_character_,
      error = function(e) if (length(assays_avail)) assays_avail[1L] else NA_character_
    )
    layers_avail <- tryCatch(
      .seurat_assay_layers(obj, default_assay),
      error = function(e) character(0)
    )
  } else if (inherits(obj, "SingleCellExperiment") ||
             inherits(obj, "SpatialExperiment")) {
    kind <- if (inherits(obj, "SpatialExperiment")) "spatial" else "sce"
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      assays_avail <- tryCatch(
        SummarizedExperiment::assayNames(obj),
        error = function(e) character(0)
      )
    }
    if (length(assays_avail)) {
      default_assay <- assays_avail[1L]
      # In SCE/SpatialExperiment each entry of `assays_avail` is itself a
      # full matrix (counts, logcounts, ...) -- they ARE the layers.
      layers_avail <- assays_avail
    }
  } else if (is.matrix(obj) || inherits(obj, "Matrix") || is.data.frame(obj)) {
    kind <- "matrix"
    n_samples <- ncol(obj)
    n_genes <- nrow(obj)
    sample_ids <- colnames(obj)
    gene_ids <- rownames(obj)
    # Plain matrix has no concept of assays / layers.
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
      layers_avail <- tryCatch(
        .anndata_layers_avail(obj),
        error = function(e) "X"
      )
      # AnnData has no "assay" concept; leave default_assay NA.
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
    gene_ids = gene_ids,
    assays_avail  = assays_avail,
    default_assay = default_assay,
    layers_avail  = layers_avail
  )
}

# Inspect a Seurat assay and return the subset of {counts, data, scale.data}
# that are present and non-empty. Supports both Assay (Seurat v4) and Assay5
# (Seurat v5) storage.
.seurat_assay_layers <- function(obj, assay_name) {
  if (is.null(obj) || is.na(assay_name) || !nzchar(assay_name)) return(character(0))
  if (!(assay_name %in% names(obj@assays))) return(character(0))
  a <- obj@assays[[assay_name]]
  candidates <- c("counts", "data", "scale.data")
  # Assay5 (Seurat v5): layers are named entries in @layers.
  if (inherits(a, "Assay5")) {
    have <- tryCatch(
      if (requireNamespace("SeuratObject", quietly = TRUE))
        SeuratObject::Layers(a) else names(a@layers),
      error = function(e) character(0)
    )
    return(intersect(candidates, have))
  }
  # Classic Assay (Seurat v3/v4): @counts / @data / @scale.data slots.
  out <- character(0)
  for (s in candidates) {
    val <- tryCatch(methods::slot(a, s), error = function(e) NULL)
    if (is.null(val)) next
    if (!is.null(dim(val)) && all(dim(val) > 0L)) out <- c(out, s)
  }
  out
}

# List the AnnData layers a user could realistically score from. Always
# includes the primary X matrix; adds named entries from adata.layers.keys()
# when present.
.anndata_layers_avail <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return("X")
  keys <- tryCatch({
    builtins <- reticulate::import_builtins(convert = FALSE)
    as.character(reticulate::py_to_r(builtins$list(obj$layers$keys())))
  }, error = function(e) character(0))
  unique(c("X", keys))
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
    res <- phenomapr_anndata_obs_df(obj)
    md  <- res$df
    if (is.null(md) || ncol(md) == 0L) {
      # Diagnostic warnings are intentionally not attached here (NULL
      # cannot carry attributes). Callers that need the failure reason
      # can re-call phenomapr_anndata_obs_df() on the same object.
      return(NULL)
    }
    md$.cell_id <- rownames(md)
    attr(md, "anndata_obs_warnings") <- res$warnings
    return(md)
  }
  NULL
}

# Resilient extractor for AnnData.obs that handles round-trip dtype
# pitfalls (categoricals, pandas extension dtypes, anndataR-written h5ad
# files) via three layered fallbacks. Returns
# list(df = <data.frame|NULL>, warnings = <character>) so the caller can
# surface a clear message when extraction degraded or failed.
phenomapr_anndata_obs_df <- function(obj) {
  warnings_acc <- character(0)
  push_w <- function(msg) {
    warnings_acc[[length(warnings_acc) + 1L]] <<- msg
  }
  result <- function(df) list(df = df, warnings = warnings_acc)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    push_w("reticulate is not installed; cannot read AnnData.obs.")
    return(result(NULL))
  }
  if (!reticulate::py_has_attr(obj, "obs")) {
    push_w("This AnnData object has no `obs` slot.")
    return(result(NULL))
  }

  finalize <- function(df) {
    if (is.null(df)) return(NULL)
    if (!is.data.frame(df)) {
      df <- tryCatch(as.data.frame(df, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (is.null(df)) return(NULL)
    }
    rn <- rownames(df)
    needs_rn <- is.null(rn) ||
      identical(rn, as.character(seq_len(nrow(df)))) ||
      anyDuplicated(rn) > 0L
    if (needs_rn) {
      ids <- tryCatch(
        as.character(reticulate::py_to_r(obj$obs_names$tolist())),
        error = function(e) NULL
      )
      if (!is.null(ids) && length(ids) == nrow(df)) {
        rownames(df) <- make.unique(ids)
      }
    }
    df
  }

  # Path 1: direct py_to_r() of the entire pandas DataFrame (fast path).
  df <- tryCatch(reticulate::py_to_r(obj$obs), error = function(e) {
    push_w(paste0("Direct py_to_r(adata.obs) failed: ", conditionMessage(e)))
    NULL
  })
  if (!is.null(df) && is.data.frame(df) && ncol(df) > 0L) {
    return(result(finalize(df)))
  }
  if (!is.null(df) && is.data.frame(df) && ncol(df) == 0L) {
    push_w("Direct py_to_r() returned a 0-column data.frame; trying decategorize fallback.")
  }

  # Path 2: cast every pandas `category` column to `object` (plain str) on
  # a copy of obs and convert that. Categorical -> R conversion is the
  # most common failure mode for AnnData files written from R via
  # anndataR / SCE -> AnnData converters that mark Seurat factor columns
  # as pandas categoricals.
  df <- tryCatch(
    .phenomapr_anndata_obs_decategorize(obj),
    error = function(e) {
      push_w(paste0("Decategorize fallback failed: ", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(df) && is.data.frame(df) && ncol(df) > 0L) {
    return(result(finalize(df)))
  }

  # Path 3: rebuild column-by-column using astype(object).tolist() and
  # tolerate per-column failures by filling that column with NA. This is
  # the slowest path but the most resilient one.
  df <- tryCatch(
    .phenomapr_anndata_obs_columnwise(obj, push_w),
    error = function(e) {
      push_w(paste0("Columnwise fallback failed: ", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(df) && is.data.frame(df) && ncol(df) > 0L) {
    return(result(finalize(df)))
  }

  push_w("All AnnData.obs extraction paths failed; no metadata recovered.")
  result(NULL)
}

.phenomapr_anndata_obs_decategorize <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  py_code <- paste(
    "def _phenomapr_decategorize(df):",
    "    out = df.copy()",
    "    for c in out.columns:",
    "        s = out[c]",
    "        dt = str(getattr(s, 'dtype', ''))",
    "        if dt == 'category':",
    "            out[c] = s.astype(object)",
    "        elif dt.startswith('Int') or dt.startswith('UInt'):",
    "            out[c] = s.astype('float').astype(object)",
    "        elif dt in ('boolean', 'string'):",
    "            out[c] = s.astype(object)",
    "    return out",
    sep = "\n"
  )
  py_env <- tryCatch(
    reticulate::py_run_string(py_code, convert = FALSE),
    error = function(e) NULL
  )
  if (is.null(py_env)) return(NULL)
  fn <- tryCatch(py_env$`_phenomapr_decategorize`, error = function(e) NULL)
  if (is.null(fn)) return(NULL)
  plain <- tryCatch(fn(obj$obs), error = function(e) NULL)
  if (is.null(plain)) return(NULL)
  reticulate::py_to_r(plain)
}

.phenomapr_anndata_obs_columnwise <- function(obj, push_w = function(msg) {}) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  cols <- tryCatch(
    as.character(reticulate::py_to_r(obj$obs$columns$tolist())),
    error = function(e) {
      push_w(paste0("Could not list obs.columns: ", conditionMessage(e)))
      character(0)
    }
  )
  if (!length(cols)) return(NULL)
  n_obs <- tryCatch(
    as.integer(reticulate::py_to_r(obj$n_obs)),
    error = function(e) NA_integer_
  )
  if (is.na(n_obs)) {
    push_w("Could not determine adata.n_obs.")
    return(NULL)
  }
  out <- vector("list", length(cols))
  names(out) <- cols
  for (col in cols) {
    val <- tryCatch({
      series <- reticulate::py_get_item(obj$obs, col)
      reticulate::py_to_r(series$astype("object")$tolist())
    }, error = function(e) {
      push_w(sprintf("obs column '%s' failed: %s", col, conditionMessage(e)))
      NULL
    })
    if (is.null(val)) {
      val <- rep(NA, n_obs)
    } else if (is.list(val)) {
      val <- vapply(val, function(x) {
        if (is.null(x)) NA_character_ else as.character(x)
      }, character(1L))
    }
    if (length(val) != n_obs) {
      push_w(sprintf("obs column '%s' had unexpected length %d (n_obs=%d).",
                     col, length(val), n_obs))
      val <- rep(NA, n_obs)
    }
    out[[col]] <- val
  }
  data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
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

# Wrap the long-running PhenoMap() call so the caller sees an elapsed-time
# attribute. The visible "we are working" indicator is now driven by the
# centered busy modal in the calling observer (see phenomapr_busy_show()),
# so this helper has no further UI side-effects.
run_phenomap_with_progress <- function(expression, reference, cancer_type,
                                       z_score_cutoff, pseudobulk, group_by,
                                       assay, slot, reference_sign) {
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
#   cell_id, dim1, dim2, dim1_name, dim2_name, is_spatial, sample
# The `sample` column is the per-cell tissue-section / library / core
# label for spatial objects with multiple slides (Seurat @images named
# entries, SpatialExperiment::colData()$sample_id, AnnData obs columns
# named library_id / sample_id / sample / slide_id / section). For
# non-spatial or single-section spatial inputs every row is tagged with a
# single section name -- this lets the Visualization tab unconditionally
# read `unique(emb$sample)` to decide whether to surface a section
# switcher in the sidebar.
extract_embedding <- function(obj, name) {
  emb <- NULL
  is_spatial <- identical(name, "spatial")
  sample_vec <- NULL   # length = nrow(emb), labels each spot's section
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    if (is_spatial) {
      # Visium / Slide-seq Seurat objects store one entry per tissue
      # section in @images, each carrying its own coords frame. We stack
      # ALL of them and tag each spot with its image name so the UI can
      # let the user flip through sections -- previously we silently used
      # only the first image, which dropped every other slide.
      #
      # `imagecol` holds the x position; `imagerow` holds the y but is
      # image-space (origin top-left). We reverse y later in the plot via
      # scale_y_reverse() rather than mutating values here, so raw
      # coordinates stay in their native frame and any axis text remains
      # accurate.
      imgs <- methods::slot(obj, "images") %||% list()
      if (length(imgs)) {
        img_names <- names(imgs)
        if (is.null(img_names) || any(!nzchar(img_names))) {
          img_names <- ifelse(nzchar(img_names %||% ""),
                              img_names,
                              paste0("image_", seq_along(imgs)))
        }
        parts <- lapply(seq_along(imgs), function(i) {
          coords <- tryCatch(imgs[[i]]@coordinates, error = function(e) NULL)
          if (is.null(coords)) return(NULL)
          if (!all(c("imagerow", "imagecol") %in% colnames(coords))) return(NULL)
          m <- as.matrix(coords[, c("imagecol", "imagerow"), drop = FALSE])
          colnames(m) <- c("x", "y")
          rownames(m) <- rownames(coords)
          list(emb = m, sample = rep(img_names[i], nrow(m)))
        })
        parts <- parts[!vapply(parts, is.null, logical(1L))]
        if (length(parts)) {
          emb <- do.call(rbind, lapply(parts, `[[`, "emb"))
          sample_vec <- unlist(lapply(parts, `[[`, "sample"), use.names = FALSE)
        }
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
      # SpatialExperiment carries `sample_id` on its colData -- use it
      # verbatim. Multi-sample objects (joinSamples, mergeSpatial
      # workflows etc.) all share one spatialCoords frame, so we lean on
      # the colData column to attach per-spot section labels.
      sample_vec <- tryCatch({
        cd <- SummarizedExperiment::colData(obj)
        if (!is.null(cd) && "sample_id" %in% colnames(cd)) {
          as.character(cd$sample_id)
        } else NULL
      }, error = function(e) NULL)
    } else {
      emb <- SingleCellExperiment::reducedDim(obj, name)
    }
  } else if (inherits(obj, "python.builtin.object")) {
    emb <- tryCatch(PhenoMapR:::.anndata_obsm_array(obj, name),
                    error = function(e) NULL)
    if (is_spatial) {
      # AnnData multi-section spatial datasets keep one big obsm["spatial"]
      # but disambiguate sections via an obs column (Squidpy convention
      # is `library_id`; Scanpy users sometimes use `sample_id` /
      # `sample` / `slide_id` / `section`). We sniff for any of those and
      # ride along on the resulting vector.
      sample_vec <- tryCatch(
        .anndata_spatial_sample_vec(obj),
        error = function(e) NULL
      )
    }
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
    ids <- tryCatch(colnames(obj), error = function(e) NULL)
    if (is.null(ids) || length(ids) != nrow(emb)) {
      ids <- paste0("Cell_", seq_len(nrow(emb)))
    }
  }
  # Normalize / fill the section label so the downstream df always has a
  # `sample` column of the right length. NULL or wrong-length -> single
  # synthetic section.
  if (is.null(sample_vec) || length(sample_vec) != nrow(emb)) {
    sample_vec <- rep(if (is_spatial) "section_1" else NA_character_,
                      nrow(emb))
  } else {
    # Promote anything non-character (factor / pandas categorical) to
    # plain character so the select input choices stay clean.
    sample_vec <- as.character(sample_vec)
    sample_vec[is.na(sample_vec) | !nzchar(sample_vec)] <- "(unknown)"
  }
  data.frame(
    cell_id = as.character(ids),
    dim1 = as.numeric(emb[, 1L]),
    dim2 = as.numeric(emb[, 2L]),
    dim1_name = d1_name,
    dim2_name = d2_name,
    is_spatial = is_spatial,
    sample = sample_vec,
    stringsAsFactors = FALSE
  )
}

# Probe an AnnData object's `obs` for the column conventionally used to
# disambiguate tissue sections in a multi-sample spatial dataset.
# Tries (in order): library_id, sample_id, sample, slide_id, section.
# Returns a character vector the same length as adata.n_obs, or NULL.
.anndata_spatial_sample_vec <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  candidates <- c("library_id", "sample_id", "sample", "slide_id",
                  "section", "section_id", "core", "core_id")
  obs_cols <- tryCatch({
    builtins <- reticulate::import_builtins(convert = FALSE)
    as.character(reticulate::py_to_r(builtins$list(obj$obs$columns)))
  }, error = function(e) character(0))
  hit <- intersect(candidates, obs_cols)
  if (!length(hit)) return(NULL)
  col <- hit[1L]
  v <- tryCatch({
    # Decategorize before py_to_r so pandas Categorical doesn't trip
    # reticulate when the user has fewer categories than reticulate
    # expects (we saw this with AnnData multi-slide datasets that store
    # library_id as Categorical with NA codes).
    builtins <- reticulate::import_builtins(convert = FALSE)
    series <- obj$obs[col]
    series <- series$astype("object")
    as.character(reticulate::py_to_r(builtins$list(series$tolist())))
  }, error = function(e) NULL)
  v
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
    `data-pfi-picker-id` = picker_id,
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
      # Wrap the uiOutput in an explicit span so the flexbox child rules
      # land on a stable element regardless of how shiny::uiOutput()
      # versions render the container's class attribute.
      shiny::tags$span(
        class = "phenomapr-file-input-name",
        shiny::uiOutput(paste0(picker_id, "_chosen"), inline = TRUE)
      )
    ),
    # Indeterminate progress bar mimicking the green "uploading" sliver
    # under shiny::fileInput(). Sliver animates whenever this file-input
    # is mid-load; toggled by the JS in phenomapr_busy_assets() in
    # response to shiny:inputchanged on `*_server` -> shiny:idle. The
    # nested inner div is what actually slides; the outer div is the
    # track.
    shiny::tags$div(
      class = "phenomapr-file-input-progress",
      shiny::tags$div(class = "phenomapr-file-input-progress-bar")
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

# ============================================================================
# Centered "busy" overlay with a 2-second activation grace period
# ============================================================================
#
# `phenomapr_busy_show()` schedules a centered, foggy-backdrop overlay with a
# spinner + message. The overlay only appears if the work has been running
# for at least `delay_seconds` (default 2s), so quick file reads / cheap
# computations don't flash a popup at the user. `phenomapr_busy_hide()`
# cancels the pending show OR removes the overlay if it's already visible.
#
# Implementation: the actual setTimeout / overlay toggle lives in the
# browser. The two custom messages we send here (`phenomapr-busy-show` and
# `phenomapr-busy-hide`) are handled by JS injected via
# phenomapr_busy_assets() in the UI head. JS handles the timing because R
# is single-threaded and won't fire `later::later()` callbacks while a
# compute-bound observer is still running.
#
# Usage pattern (mirroring the old shiny::Progress$new() pattern):
#
#   phenomapr_busy_show("Computing PhenoMap scores...", "Cancer: LUAD")
#   on.exit(phenomapr_busy_hide(), add = TRUE)
#   res <- PhenoMap(...)
#
# Calling `phenomapr_busy_show()` while a previous show is still pending or
# visible replaces the message (the JS handler clears its previous timer).

phenomapr_busy_show <- function(message,
                                detail = NULL,
                                hint = "This may take a moment...",
                                delay_seconds = 2,
                                session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) return(invisible(NULL))
  # Note: the popup's *visibility* is now driven entirely by Shiny's
  # shiny:busy / shiny:idle duration on the client side. This R-side
  # call only sets the message that gets displayed if/when the popup
  # appears (i.e. when the server has been continuously busy for
  # > delay_seconds). Existing observer-side usage:
  #
  #   phenomapr_busy_show("Computing PhenoMap scores...")
  #   on.exit(phenomapr_busy_hide(), add = TRUE)
  #
  # remains correct: the show updates the in-flight message text, the
  # on.exit hide reverts it back to a neutral default once the work is
  # done. Fast operations cannot trigger the popup at all because Shiny
  # goes idle before the 2s show timer fires.
  session$sendCustomMessage("phenomapr-busy-show", list(
    message  = as.character(message %||% ""),
    detail   = if (is.null(detail) || !nzchar(detail))   "" else as.character(detail),
    hint     = if (is.null(hint)   || !nzchar(hint))     "" else as.character(hint),
    delay_ms = as.integer(round(as.numeric(delay_seconds) * 1000))
  ))
  invisible(NULL)
}

phenomapr_busy_hide <- function(session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) return(invisible(NULL))
  session$sendCustomMessage("phenomapr-busy-hide", list())
  invisible(NULL)
}

# Inject the busy-overlay markup, custom-message handlers, and DOM bootstrap
# into the page <head>. Call from page_navbar()'s `header = tags$head(...)`.
#
# Implementation notes:
#   - Counter-based instead of token-based: every "show" increments an
#     activeOps counter, every "hide" decrements it. The overlay is
#     allowed to render only while activeOps > 0; the moment it drops to
#     zero we cancel any pending timer AND hide synchronously. This is
#     immune to message ordering issues (e.g. messages arriving with
#     >2s of latency between them) because the counter, not a timer
#     race, gates rendering.
#   - shiny:idle safety net: if Shiny tells us it's idle but our counter
#     is still positive (meaning some show was queued without a matching
#     hide -- a server-side bug, never a normal flow), we forcibly reset
#     the counter and hide. This guarantees the overlay can never get
#     stuck on screen.
#   - Append `?busyDebug=1` to the URL to enable console.debug logging
#     of every show/hide message and every state transition. Useful when
#     diagnosing future timing reports without touching the code.
phenomapr_busy_assets <- function() {
  shiny::tags$script(shiny::HTML(
'
(function () {
  var DEBUG = false;
  try { DEBUG = window.location.search.indexOf("busyDebug=1") !== -1; }
  catch (e) { DEBUG = false; }
  function dbg() {
    if (!DEBUG) return;
    try { console.debug.apply(console, ["[phenomapr-busy]"].concat([].slice.call(arguments))); }
    catch (e) {}
  }

  // ---- Tunables ------------------------------------------------------
  // The popup can become visible via TWO independent paths -- whichever
  // criterion is met first wins. Both paths converge on the same hide
  // logic so there is no way for the popup to get "stuck":
  //
  //   Path A: shiny:busy duration  (event-driven)
  //     pollVisibility() reveals the popup the first time Shiny has
  //     been continuously busy for >= SHOW_DELAY_MS. Robust for any
  //     long server op that produces normal busy/idle events.
  //
  //   Path B: R-side phenomapr_busy_show() timer  (message-driven)
  //     When R calls phenomapr_busy_show("..."), the client schedules
  //     a setTimeout(SHOW_DELAY_MS) that reveals the popup if it has
  //     not been cancelled by then. Catches long ops on environments
  //     where shiny:busy events are unreliable.
  //
  // Hides ALWAYS happen (in addition to events) via:
  //   * shiny:idle after IDLE_GRACE_MS of stable idleness
  //   * R-side phenomapr_busy_hide() (clears the R timer + hides the
  //     popup if visible -- this is the critical bug fix for the
  //     "popup stuck for a minute" regression).
  //   * User dismiss (close button / backdrop click)
  var SHOW_DELAY_MS      = 2000;   // both paths use the same threshold
  var IDLE_GRACE_MS      = 250;    // tolerate brief idle bounces between chained reactives
  var POLL_INTERVAL_MS   = 100;    // how often the busy-duration polling fires
  var ABS_VISIBLE_MAX_MS = 30000;  // belt-and-suspenders: popup auto-hides this long after appearing
                                   // (only meant to cover truly broken sessions; not a normal hide path)
  var FILE_INPUT_SUFFIX  = "_server"; // shinyFiles inputs created by phenomapr_file_input

  // ---- DOM helpers ---------------------------------------------------
  function ensureOverlay() {
    if (document.getElementById("phenomapr-busy-overlay")) return;
    var overlay = document.createElement("div");
    overlay.id = "phenomapr-busy-overlay";
    overlay.className = "phenomapr-busy-overlay";
    overlay.setAttribute("role", "alertdialog");
    overlay.setAttribute("aria-modal", "true");
    overlay.setAttribute("aria-live", "polite");
    overlay.innerHTML = ""
      + "<div class=\\"phenomapr-busy-card\\">"
      + "  <button type=\\"button\\" class=\\"phenomapr-busy-close\\" id=\\"phenomapr-busy-close\\" aria-label=\\"Dismiss\\" title=\\"Dismiss\\">&times;</button>"
      + "  <div class=\\"phenomapr-busy-spinner\\"></div>"
      + "  <div class=\\"phenomapr-busy-text\\">"
      + "    <div class=\\"phenomapr-busy-message\\" id=\\"phenomapr-busy-message\\"></div>"
      + "    <div class=\\"phenomapr-busy-detail\\" id=\\"phenomapr-busy-detail\\"></div>"
      + "    <div class=\\"phenomapr-busy-hint\\" id=\\"phenomapr-busy-hint\\"></div>"
      + "  </div>"
      + "</div>";
    document.body.appendChild(overlay);
  }
  function setText(id, val) {
    var el = document.getElementById(id);
    if (!el) return;
    if (val && String(val).length) {
      el.textContent = val;
      el.style.display = "block";
    } else {
      el.textContent = "";
      el.style.display = "none";
    }
  }
  function applyOverlayText() {
    ensureOverlay();
    setText("phenomapr-busy-message", currentMessage.message || "Working...");
    setText("phenomapr-busy-detail",  currentMessage.detail  || "");
    setText("phenomapr-busy-hint",    currentMessage.hint    || "This may take a moment...");
  }
  function renderShow() {
    applyOverlayText();
    var ov = document.getElementById("phenomapr-busy-overlay");
    if (ov && !ov.classList.contains("is-visible")) {
      ov.classList.add("is-visible");
      isVisible = true;
      shownAt = Date.now();
      dbg("overlay shown", currentMessage);
    }
  }
  function renderHide() {
    var ov = document.getElementById("phenomapr-busy-overlay");
    if (ov && ov.classList.contains("is-visible")) {
      ov.classList.remove("is-visible");
      isVisible = false;
      shownAt = null;
      dbg("overlay hidden");
    }
  }

  // ---- State ---------------------------------------------------------
  // The popup is driven entirely by:
  //   busyStartedAt    -- timestamp (Date.now()) of the most recent
  //                       transition from idle to busy, or null when
  //                       Shiny is (effectively) idle. The polling loop
  //                       below checks (Date.now() - busyStartedAt) and
  //                       reveals the popup once that crosses
  //                       SHOW_DELAY_MS. Short idle bounces inside one
  //                       composite op tolerate themselves via the
  //                       grace timer instead of resetting this stamp.
  //   idleGraceTimer   -- setTimeout pending the IDLE_GRACE_MS debounce
  //                       on shiny:idle. When it fires we declare
  //                       Shiny truly idle, clear busyStartedAt, and
  //                       hide the popup if it was visible.
  //   isVisible        -- mirrors the .is-visible CSS class so we do
  //                       not repeatedly call renderShow/renderHide.
  //   dismissedThisRun -- true between user-dismiss and the next
  //                       genuine idle->busy transition. Suppresses
  //                       re-showing the popup for the SAME op after a
  //                       manual dismiss.
  var currentMessage    = { message: "Working...", detail: "", hint: "" };
  var defaultMessage    = { message: "Working...", detail: "", hint: "" };
  var busyStartedAt     = null;     // Path A: when shiny:busy first fired this run
  var rSideShowTimer    = null;     // Path B: pending setTimeout from R show
  var idleGraceTimer    = null;
  var isVisible         = false;
  var shownAt           = null;
  var dismissedThisRun  = false;
  // suppressUntilIdle: when R has explicitly said hide (or the user
  // dismissed), block Path A (shiny:busy-duration) from re-revealing
  // the popup until the next confirmed idle. Without this, downstream
  // reactives that fire shiny:busy after R hide (e.g. plot/table
  // rendering after a long observer returns) would re-show the popup
  // because busyStartedAt has not yet been cleared by an idle grace.
  // Cleared on idle grace, disconnect, or an explicit R show.
  var suppressUntilIdle = false;

  function clearIdleGraceTimer() {
    if (idleGraceTimer !== null) {
      clearTimeout(idleGraceTimer);
      idleGraceTimer = null;
    }
  }
  function clearRSideShowTimer() {
    if (rSideShowTimer !== null) {
      clearTimeout(rSideShowTimer);
      rSideShowTimer = null;
    }
  }
  function resetMessage() {
    currentMessage = {
      message: defaultMessage.message,
      detail:  defaultMessage.detail,
      hint:    defaultMessage.hint
    };
  }

  // ---- Shiny busy/idle handlers --------------------------------------
  // shiny:busy fires whenever Shiny starts processing input changes.
  // We record the transition time but never schedule a setTimeout to
  // reveal the popup -- the polling loop below does that based on
  // sustained busy duration.
  function onShinyBusy() {
    clearIdleGraceTimer();
    if (busyStartedAt === null) {
      busyStartedAt = Date.now();
      dbg("busy started at " + busyStartedAt);
    }
  }
  // shiny:idle fires when Shiny finishes the current cycle. We do
  // NOT act immediately -- short bounces are normal in chained reactives.
  // After IDLE_GRACE_MS of stable idleness we declare the op truly
  // over: clear busyStartedAt, cancel any pending R-side show, hide
  // the popup if visible, reset the dismissed flag for the next op.
  function onShinyIdle() {
    if (idleGraceTimer !== null) return;
    idleGraceTimer = setTimeout(function () {
      idleGraceTimer = null;
      busyStartedAt = null;
      clearRSideShowTimer();
      if (isVisible) renderHide();
      clearAllFileInputLoading();
      dismissedThisRun = false;
      suppressUntilIdle = false;
      resetMessage();
      dbg("idle confirmed: popup hidden, state reset");
    }, IDLE_GRACE_MS);
  }
  function onShinyDisconnected() {
    clearIdleGraceTimer();
    clearRSideShowTimer();
    busyStartedAt = null;
    if (isVisible) renderHide();
    clearAllFileInputLoading();
    dismissedThisRun = false;
    suppressUntilIdle = false;
    resetMessage();
    dbg("disconnected: forced hide");
  }

  // ---- Polling loop --------------------------------------------------
  // Runs every POLL_INTERVAL_MS. The job is purely to flip the popup
  // ON when shiny:busy has been continuously true for >= SHOW_DELAY_MS.
  // Hiding is event-driven (idle grace timer) -- the poll also serves
  // as a watchdog against the popup getting visually stuck longer than
  // ABS_VISIBLE_MAX_MS for any reason.
  function pollVisibility() {
    var now = Date.now();
    if (busyStartedAt !== null &&
        !isVisible &&
        !dismissedThisRun &&
        !suppressUntilIdle &&
        (now - busyStartedAt) >= SHOW_DELAY_MS) {
      renderShow();
    }
    if (isVisible && shownAt !== null &&
        (now - shownAt) > ABS_VISIBLE_MAX_MS) {
      dbg("watchdog: visible > " + ABS_VISIBLE_MAX_MS + "ms; force-hide");
      renderHide();
      busyStartedAt = null;
    }
  }

  // ---- User-driven dismiss -------------------------------------------
  function userDismiss(reason) {
    dbg("user dismissed", reason);
    clearRSideShowTimer();
    if (isVisible) renderHide();
    busyStartedAt = null;
    dismissedThisRun = true;
    suppressUntilIdle = true;
  }

  // ---- R-side message API (Path B) -----------------------------------
  // phenomapr_busy_show("...") updates the popup text AND schedules a
  // setTimeout(SHOW_DELAY_MS) that reveals the popup if it has not
  // been cancelled by then. The matching on.exit phenomapr_busy_hide()
  // unconditionally cancels the timer AND hides the popup if visible.
  //
  // CRITICAL BUG FIX: handleHide must hide the popup, not just cancel
  // the timer. Without this, an environment-specific race could let
  // the R-side setTimeout fire AFTER the work was done -- the popup
  // would appear, no shiny:idle would be coming to clear it, and it
  // would sit on screen until ABS_VISIBLE_MAX_MS. Hiding from
  // handleHide closes that window: even in the worst-case ordering,
  // R-side hide guarantees the popup goes away promptly.
  function maybeShowR() {
    rSideShowTimer = null;
    if (!dismissedThisRun && !isVisible) {
      dbg("R-side show timer fired");
      renderShow();
    }
  }
  function handleShow(p) {
    currentMessage = {
      message: (p && p.message) || "Working...",
      detail:  (p && p.detail)  || "",
      hint:    (p && p.hint)    || "This may take a moment..."
    };
    // An explicit R show overrides any previous suppress -- the
    // server has just begun a new piece of work and the popup is
    // welcome again.
    suppressUntilIdle = false;
    if (isVisible) applyOverlayText();
    if (rSideShowTimer === null && !dismissedThisRun && !isVisible) {
      var d = (p && typeof p.delay_ms === "number" && p.delay_ms >= 0)
        ? p.delay_ms : SHOW_DELAY_MS;
      rSideShowTimer = setTimeout(maybeShowR, d);
      dbg("R-side show: timer scheduled in " + d + "ms");
    }
  }
  function handleHide() {
    clearRSideShowTimer();
    // (1) Cancel the pending R-side timer if it has not fired.
    // (2) Hide the popup immediately if visible -- closes the prior
    //     "stuck for a minute" bug where a stale timer could fire
    //     after the work was already done.
    // (3) Reset busyStartedAt and arm suppressUntilIdle so Path A
    //     does NOT immediately re-show the popup when downstream
    //     reactives (plot/table rendering after the originating
    //     observer returns) trigger more shiny:busy events. Without
    //     this, the popup would re-flash within ~100ms and stay
    //     until the next clean idle, which the user perceives as
    //     "the popup lingers several seconds after the work is done".
    if (isVisible) renderHide();
    busyStartedAt = null;
    suppressUntilIdle = true;
    resetMessage();
    dbg("R-side hide: timer cleared + popup hidden + Path A suppressed until next idle");
  }

  // ---- File-input progress feedback ----------------------------------
  // Two-phase visual loader under the Browse row:
  //
  //   .is-loading           -> CSS fills the bar 0% -> 95% over 30s
  //                            (see phenomapr-file-input-fill keyframe).
  //   .is-loading-complete  -> bar leaps to 100% and fades; CSS handles
  //                            the transition.
  //
  // The pipeline:
  //   1. shiny:inputchanged on `<id>_server` (a phenomapr_file_input
  //      picker change) -> add `.is-loading` to that file-input element.
  //   2. shiny:idle (after IDLE_GRACE_MS) -> swap `.is-loading` for
  //      `.is-loading-complete` on every still-loading file-input.
  //   3. After FILE_INPUT_COMPLETE_FADE_MS -> remove
  //      `.is-loading-complete` so the file-input is back to its
  //      default invisible-bar state.
  var FILE_INPUT_COMPLETE_FADE_MS = 900;

  function clearAllFileInputLoading() {
    // Move every currently-loading file-input into the completing phase
    // (bar -> 100% + fade). Then schedule a final cleanup that strips
    // the helper classes entirely so the next file pick starts fresh.
    var loading = document.querySelectorAll(".phenomapr-file-input.is-loading");
    for (var i = 0; i < loading.length; i++) {
      loading[i].classList.remove("is-loading");
      loading[i].classList.add("is-loading-complete");
    }
    if (loading.length > 0) {
      setTimeout(function () {
        var completing = document.querySelectorAll(
          ".phenomapr-file-input.is-loading-complete"
        );
        for (var j = 0; j < completing.length; j++) {
          completing[j].classList.remove("is-loading-complete");
        }
      }, FILE_INPUT_COMPLETE_FADE_MS);
    }
  }
  function markFileInputLoading(pickerId) {
    if (!pickerId) return;
    var el = document.querySelector(
      ".phenomapr-file-input[data-pfi-picker-id=\\"" + pickerId + "\\"]"
    );
    if (el) {
      // If a previous load is still in its fade-out phase, strip it so
      // the new bar starts from 0% rather than 100%.
      el.classList.remove("is-loading-complete");
      // Re-trigger the keyframe animation by stripping and re-adding
      // the loading class. Without the void access browsers may not
      // restart the animation if .is-loading was already present.
      el.classList.remove("is-loading");
      // Force a reflow so the animation restarts on re-add.
      // (Reading offsetWidth is a well-known reflow trick.)
      void el.offsetWidth;
      el.classList.add("is-loading");
      dbg("file-input loading start: " + pickerId);
    }
  }
  function onInputChanged(evt) {
    // Shiny fires this on every input value update. We only react to
    // names matching the phenomapr_file_input picker convention
    // (`<id>_server`) so other input changes do not accidentally flash
    // the bar.
    var nm = (evt && evt.detail && evt.detail.name) ||
             (evt && evt.name);
    if (!nm) return;
    if (nm.length <= FILE_INPUT_SUFFIX.length) return;
    if (nm.lastIndexOf(FILE_INPUT_SUFFIX) !==
        nm.length - FILE_INPUT_SUFFIX.length) return;
    markFileInputLoading(nm);
  }

  // ---- Polling loop --------------------------------------------------
  // pollVisibility() (defined above) inspects busyStartedAt every
  // POLL_INTERVAL_MS and flips the popup ON once the busy duration
  // crosses SHOW_DELAY_MS. It also force-hides if the popup somehow
  // stays visible past ABS_VISIBLE_MAX_MS (belt-and-suspenders).
  setInterval(pollVisibility, POLL_INTERVAL_MS);

  // ---- Wire up event listeners ---------------------------------------
  // Native fallback layer. Shiny normally dispatches via jQuery, and
  // jquery/jquery#2476 means $(document).trigger() does not reliably
  // invoke native addEventListener handlers for namespaced custom
  // events. But many Shiny+browser combos do bubble these events to
  // document, so we register on the native API too -- whichever layer
  // works first wins; both firing is harmless because all handlers are
  // idempotent.
  if (document && document.addEventListener) {
    document.addEventListener("shiny:busy",         onShinyBusy,         false);
    document.addEventListener("shiny:idle",         onShinyIdle,         false);
    document.addEventListener("shiny:disconnected", onShinyDisconnected, false);
    document.addEventListener("shiny:inputchanged", onInputChanged,      false);
  }

  function attachDismissHandlers() {
    var ov = document.getElementById("phenomapr-busy-overlay");
    if (ov && !ov._dismissBound) {
      ov._dismissBound = true;
      ov.addEventListener("click", function (e) {
        if (e.target === ov) userDismiss("backdrop-click");
      }, false);
    }
    var btn = document.getElementById("phenomapr-busy-close");
    if (btn && !btn._dismissBound) {
      btn._dismissBound = true;
      btn.addEventListener("click", function (e) {
        e.stopPropagation();
        userDismiss("close-button");
      }, false);
    }
  }

  function bind() {
    if (typeof Shiny === "undefined" || !Shiny.addCustomMessageHandler) {
      setTimeout(bind, 50);
      return;
    }
    ensureOverlay();
    attachDismissHandlers();
    Shiny.addCustomMessageHandler("phenomapr-busy-show", handleShow);
    Shiny.addCustomMessageHandler("phenomapr-busy-hide", handleHide);

    var $j = window.jQuery || window.$;
    if ($j) {
      $j(document).on("shiny:busy",         onShinyBusy);
      $j(document).on("shiny:idle",         onShinyIdle);
      $j(document).on("shiny:disconnected", onShinyDisconnected);
      $j(document).on("shiny:inputchanged", onInputChanged);
    } else {
      dbg("jQuery not available at bind(); relying on native listeners + watchdog");
    }
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", bind);
  } else {
    bind();
  }
})();
'
  ))
}
