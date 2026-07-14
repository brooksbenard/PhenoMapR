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

# User-facing phenotype tail labels (red = +, blue = −), shared app-wide.
.phenomapr_phenotype_plus_class  <- "phenomapr-phenotype-plus"
.phenomapr_phenotype_minus_class <- "phenomapr-phenotype-minus"

phenomapr_phenotype_plus <- function(text = "Most Phenotype +") {
  shiny::tags$span(class = .phenomapr_phenotype_plus_class, text)
}

phenomapr_phenotype_minus <- function(text = "Most Phenotype \u2212") {
  shiny::tags$span(class = .phenomapr_phenotype_minus_class, text)
}

phenomapr_phenotype_plus_html <- function(text = "Most Phenotype +") {
  htmltools::HTML(sprintf(
    '<span class="%s">%s</span>',
    .phenomapr_phenotype_plus_class,
    htmltools::htmlEscape(text)
  ))
}

phenomapr_phenotype_minus_html <- function(text = "Most Phenotype \u2212") {
  htmltools::HTML(sprintf(
    '<span class="%s">%s</span>',
    .phenomapr_phenotype_minus_class,
    htmltools::htmlEscape(text)
  ))
}

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
  # Cell-content rendering rule:
  #   * If the cell is already a Shiny tag / tagList / HTML object,
  #     embed it verbatim (so callers can pass shiny::HTML("<strong>..."),
  #     icon(), tagList(...), etc. without it being escaped into text).
  #   * Otherwise stringify with as.character() -- Shiny escapes the
  #     result as a plain text node, which is the safe default.
  render_cell <- function(val) {
    if (inherits(val, c("shiny.tag", "shiny.tag.list", "html"))) val
    else as.character(val)
  }
  # df[[cn]] returns the underlying vector OR list for the named
  # column; `[[i]]` then yields a scalar character (regular column)
  # or the original object (list-column). Using this two-step
  # accessor instead of `df[i, cn]` avoids the data.frame's
  # row-indexing path, which wraps list-column results in another
  # list and breaks the inherits() check above.
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
            shiny::tags$td(render_cell(df[[cn]][[i]]))
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
  # Two-stage strategy:
  #   1) PREFER the pure-R hdf5r route. h5ad files are plain HDF5,
  #      so hdf5r can pull the X matrix and the obs / var index
  #      columns out itself and return a (genes x cells) dgCMatrix
  #      with HUGO-style rownames and cell-id colnames, plus an
  #      attached `phenomapr_obs` data.frame (read out of /obs via
  #      hdf5r as well). This flows straight into the matrix path
  #      the rest of the app already handles AND keeps cell-level
  #      metadata available without ever touching Python -- which
  #      is what most users actually want (no reticulate / Python
  #      install, no "Could not import anndata" errors, no
  #      "find_phenotype_markers failed: 'expression' must be a
  #      matrix" downstream because the markers path doesn't accept
  #      python.builtin.object directly).
  #   2) FALL BACK to the Python anndata route (via reticulate)
  #      only when hdf5r is NOT installed. Returns a fully-fledged
  #      AnnData object that the AnnData-aware helpers
  #      (process_anndata, .anndata_X_to_genes_cells, the obsm UMAP
  #      detection, etc.) consume. Slightly richer (obsm /
  #      additional layers) but requires a Python toolchain.
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    return(.read_h5ad_hdf5r(path))
  }
  if (requireNamespace("reticulate", quietly = TRUE)) {
    ad <- tryCatch(reticulate::import("anndata", convert = FALSE),
                   error = function(e) NULL)
    if (!is.null(ad)) {
      return(tryCatch(ad$read_h5ad(path),
                      error = function(e) {
                        stop(
                          "Reading .h5ad via Python anndata failed (",
                          conditionMessage(e),
                          "). Install the `hdf5r` R package for a ",
                          "Python-free reader: install.packages(\"hdf5r\")."
                        )
                      }))
    }
  }
  .read_h5ad_hdf5r(path)  # triggers a clean "install hdf5r" error.
}

# Python-free h5ad reader. Pulls the (cells x genes) X matrix out of
# the HDF5 file via hdf5r and returns it transposed to the
# (genes x cells) orientation PhenoMapR uses everywhere else, with
# `var/_index` (or `var/index`) supplying gene rownames and
# `obs/_index` (or `obs/index`) supplying cell colnames. Handles the
# two on-disk X layouts AnnData writes:
#   - sparse CSR / CSC: a group containing data + indices + indptr
#   - dense:           a 2-D dataset
.read_h5ad_hdf5r <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop(
      "Reading .h5ad files requires the `hdf5r` R package. ",
      "Install with: install.packages(\"hdf5r\")."
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required to read .h5ad")
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  if (!"X" %in% names(h5)) {
    stop("h5ad file is missing the /X dataset (no expression matrix).")
  }
  X_node <- h5[["X"]]

  # ---- 1) Read X and build a genes x cells matrix directly ------------
  #
  # We always end up wanting a (n_vars x n_obs) = (genes x cells)
  # matrix to match PhenoMapR's convention. AnnData stores X on disk
  # as (n_obs x n_vars), in one of three layouts:
  #
  #   * CSR (encoding-type = "csr_matrix", default for sparse):
  #       indices are column indices into n_vars,
  #       indptr  partitions rows of n_obs.
  #     The same data/indices/indptr arrays interpreted as a CSC of the
  #     TRANSPOSED matrix give us a genes x cells dgCMatrix directly,
  #     with no triplet/CsparseMatrix intermediate. This is the same
  #     CSR-as-CSC-of-transpose trick used in
  #     .anndata_X_to_genes_cells() and avoids the ~3x peak-RAM blow-up
  #     of the previous repr = "T" -> CsparseMatrix path (which was
  #     blowing past R's default 16 GB limit for a 2.6 GB h5ad).
  #
  #   * CSC (encoding-type = "csc_matrix"): build the (n_obs x n_vars)
  #     dgCMatrix as-is and Matrix::t() it. The transpose is O(nnz)
  #     and doesn't need a triplet intermediate.
  #
  #   * Dense H5D dataset: hdf5r reads it as a regular matrix in
  #     (n_obs, n_vars) orientation; t() to get genes x cells.
  X_genes_x_cells <- if (inherits(X_node, "H5Group")) {
    data    <- X_node[["data"]]$read()
    indices <- as.integer(X_node[["indices"]]$read())
    indptr  <- as.integer(X_node[["indptr"]]$read())
    shape   <- tryCatch(
      as.integer(hdf5r::h5attr(X_node, "shape")),
      error = function(e) {
        if ("shape" %in% names(X_node)) as.integer(X_node[["shape"]]$read())
        else NULL
      }
    )
    enc <- tryCatch(
      hdf5r::h5attr(X_node, "encoding-type"),
      error = function(e) NA_character_
    )
    if (is.null(shape) || length(shape) != 2L) {
      stop("h5ad sparse /X has no readable `shape` attribute.")
    }
    if (identical(as.character(enc), "csc_matrix")) {
      # CSC of (n_obs x n_vars). Build as cells x genes, then transpose
      # to genes x cells. dgCMatrix transpose is cheap (O(nnz)).
      m_csc <- Matrix::sparseMatrix(
        i = indices + 1L, p = indptr, x = data,
        dims = c(shape[1L], shape[2L]),
        index1 = TRUE
      )
      Matrix::t(m_csc)
    } else {
      # CSR of (n_obs x n_vars) -- same arrays are CSC of the
      # transposed (n_vars x n_obs) matrix, which is exactly what we
      # want for genes x cells. Treat unspecified encoding the same
      # way (CSR is AnnData's default).
      Matrix::sparseMatrix(
        i = indices + 1L,        # column indices in the original CSR
                                 #   == row indices in the CSC of transpose
        p = indptr,
        x = data,
        dims = c(shape[2L], shape[1L]),  # n_vars x n_obs = genes x cells
        index1 = TRUE
      )
    }
  } else {
    # Dense H5D dataset in (n_obs, n_vars) orientation.
    arr <- X_node$read()
    if (!is.matrix(arr)) arr <- as.matrix(arr)
    t(arr)
  }

  # ---- 2) Extract gene names from /var --------------------------------
  read_index <- function(grp_name) {
    if (!grp_name %in% names(h5)) return(NULL)
    grp <- h5[[grp_name]]
    if (!inherits(grp, "H5Group")) return(NULL)
    candidates <- c("_index", "index")
    pick <- candidates[candidates %in% names(grp)]
    if (!length(pick)) {
      # AnnData also encodes the index column name as an attr.
      idx_attr <- tryCatch(hdf5r::h5attr(grp, "_index"),
                           error = function(e) NULL)
      if (!is.null(idx_attr) && idx_attr %in% names(grp)) pick <- idx_attr
    }
    if (!length(pick)) return(NULL)
    as.character(grp[[pick[[1L]]]]$read())
  }
  var_names <- read_index("var")
  obs_names <- read_index("obs")

  # X_genes_x_cells is already in (n_vars x n_obs) = genes x cells
  # orientation -- so var_names go onto rows, obs_names onto columns.
  if (!is.null(var_names) && length(var_names) == nrow(X_genes_x_cells)) {
    rownames(X_genes_x_cells) <- var_names
  }
  if (!is.null(obs_names) && length(obs_names) == ncol(X_genes_x_cells)) {
    colnames(X_genes_x_cells) <- obs_names
  }

  # ---- 3) Read /obs columns (cell-level metadata) ---------------------
  # Pull obs out via hdf5r so users get the AnnData metadata without
  # ever touching Python. Quietly skip on any failure -- the matrix
  # itself is the load-bearing return value, metadata is a bonus.
  obs_df <- tryCatch(
    .read_h5ad_obs_hdf5r(h5, expected_n = ncol(X_genes_x_cells)),
    error = function(e) NULL
  )
  if (!is.null(obs_df) && nrow(obs_df) == ncol(X_genes_x_cells)) {
    # Always replace the synthetic 1:N rownames `as.data.frame()`
    # auto-generates with the real cell IDs from /obs/_index, so
    # downstream extract_object_metadata() can populate `.cell_id`
    # correctly.
    if (!is.null(obs_names) && length(obs_names) == nrow(obs_df)) {
      rownames(obs_df) <- obs_names
    }
  }

  # Stash the obs data.frame on the returned matrix so
  # extract_object_metadata() can surface it back to the Shiny app
  # without re-opening the .h5ad file. attr<- works on both regular
  # matrices and dgCMatrix (it attaches via the S4 slot, but Matrix's
  # print/operator methods don't choke on the extra attribute).
  if (!is.null(obs_df)) {
    attr(X_genes_x_cells, "phenomapr_obs") <- obs_df
  }

  # ---- 4) Read /obsm entries (cell-level embeddings) ------------------
  # When we switched the default loader to hdf5r (Python-free), we
  # stopped going through the Python AnnData object that previously
  # provided list_available_embeddings()/extract_embedding() with live
  # access to obsm["X_umap"], obsm["spatial"], obsm["segmentation"],
  # etc. Without this step the Visualization tab would show an empty
  # Reduction dropdown for every h5ad upload -- a silent regression
  # that lost UMAP / spatial / segmentation views for users.
  #
  # We pre-extract every 2D obsm array into a named list and stash it
  # as the `phenomapr_obsm` attribute. list_available_embeddings()
  # surfaces those keys for plain matrices, and extract_embedding()
  # round-trips the corresponding cells x 2 matrix back as a tidy
  # data.frame for plotting. Quietly skip on any failure -- losing an
  # obsm entry should never break the matrix itself.
  obsm_list <- tryCatch(
    .read_h5ad_obsm_hdf5r(h5, expected_n = ncol(X_genes_x_cells),
                          obs_names = obs_names),
    error = function(e) NULL
  )
  if (!is.null(obsm_list) && length(obsm_list)) {
    attr(X_genes_x_cells, "phenomapr_obsm") <- obsm_list
  }
  if (is.data.frame(obs_df)) {
    seg_mat <- .build_segmentation_obsm_from_obs(obs_df, cell_ids = obs_names)
    if (!is.null(seg_mat)) {
      obsm_list <- attr(X_genes_x_cells, "phenomapr_obsm") %||% list()
      if (!length(intersect(.anndata_segmentation_obsm_keys, names(obsm_list)))) {
        obsm_list$segmentation <- seg_mat
        attr(X_genes_x_cells, "phenomapr_obsm") <- obsm_list
      }
    }
  }

  X_genes_x_cells
}

# Pure-R reader for /obsm in an open hdf5r::H5File. Returns a named
# list of (n_obs x k) numeric matrices for every entry under /obsm
# whose second dimension is at least 2 (we only need 2D embeddings;
# higher-rank arrays like obsm["X_pca"] of width 50 also pass through
# with their full width, leaving the first two columns to be picked up
# by extract_embedding()). Rownames on each matrix are populated from
# /obs/_index when available so cell IDs round-trip through the
# Visualization tab. Quietly returns NULL when /obsm is missing or
# unreadable so the matrix path keeps working in stripped-down h5ads.
.read_h5ad_obsm_hdf5r <- function(h5, expected_n, obs_names = NULL) {
  if (!"obsm" %in% names(h5)) return(NULL)
  grp <- h5[["obsm"]]
  if (!inherits(grp, "H5Group")) return(NULL)
  keys <- names(grp)
  if (!length(keys)) return(NULL)
  out <- list()
  for (k in keys) {
    arr <- tryCatch(grp[[k]]$read(), error = function(e) NULL)
    if (is.null(arr)) next
    if (!is.numeric(arr) && !is.integer(arr)) next
    if (!is.matrix(arr)) {
      # 1D obsm entries (rare but legal) aren't useful as 2D embeddings.
      arr <- tryCatch(as.matrix(arr), error = function(e) NULL)
      if (is.null(arr)) next
    }
    # AnnData stores obsm as (n_obs, k); hdf5r reads HDF5 row-major
    # arrays back as transposed, so we end up with a (k, n_obs)
    # matrix. Detect-and-flip rather than blindly transposing so we
    # don't break any writer that already lays its arrays out the
    # R-friendly way.
    if (nrow(arr) != expected_n && ncol(arr) == expected_n) {
      arr <- t(arr)
    }
    if (nrow(arr) != expected_n || ncol(arr) < 2L) next
    if (!is.null(obs_names) && length(obs_names) == nrow(arr)) {
      rownames(arr) <- obs_names
    }
    if (is.null(colnames(arr)) || !length(colnames(arr))) {
      colnames(arr) <- paste0(k, "_", seq_len(ncol(arr)))
    }
    out[[k]] <- arr
  }
  if (!length(out)) NULL else out
}

# Pure-R reader for /obs in an open hdf5r::H5File. Walks every column
# entry under /obs (excluding _index / index, which we already use for
# rownames) and converts it into a data.frame column. AnnData v0.8+
# encodes most columns as plain HDF5 datasets, but categorical columns
# arrive as a sub-group with /categories + /codes. We handle both.
.read_h5ad_obs_hdf5r <- function(h5, expected_n) {
  if (!"obs" %in% names(h5)) return(NULL)
  grp <- h5[["obs"]]
  if (!inherits(grp, "H5Group")) return(NULL)

  # Skip the index dataset(s) used elsewhere for rownames.
  idx_attr <- tryCatch(hdf5r::h5attr(grp, "_index"),
                       error = function(e) NULL)
  skip <- unique(c("_index", "index", idx_attr))
  cols <- setdiff(names(grp), skip)
  if (!length(cols)) return(NULL)

  decode_col <- function(name) {
    node <- grp[[name]]
    if (inherits(node, "H5Group")) {
      # Categorical: AnnData stores `categories` (the levels) and
      # `codes` (0-based integer indices into categories, with -1
      # meaning NA). Older / non-anndataR writers also use a
      # `values` dataset for the integer codes.
      if (all(c("categories", "codes") %in% names(node))) {
        cats  <- as.character(node[["categories"]]$read())
        codes <- as.integer(node[["codes"]]$read())
        x <- ifelse(codes < 0L, NA_character_, cats[codes + 1L])
        return(factor(x, levels = cats))
      }
      return(NULL)
    }
    if (inherits(node, "H5D")) {
      v <- tryCatch(node$read(), error = function(e) NULL)
      if (is.null(v)) return(NULL)
      # hdf5r returns variable-length strings as a list-of-character
      # in some encodings; collapse to a plain character vector.
      if (is.list(v) && all(vapply(v, is.character, logical(1)))) {
        v <- vapply(v, function(z) if (length(z)) z[[1L]] else NA_character_,
                    character(1))
      }
      return(v)
    }
    NULL
  }

  parts <- lapply(cols, decode_col)
  names(parts) <- cols
  parts <- parts[vapply(parts,
                        function(v) !is.null(v) && length(v) == expected_n,
                        logical(1))]
  if (!length(parts)) return(NULL)
  as.data.frame(parts, stringsAsFactors = FALSE)
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
  # Matrix / dgCMatrix / data.frame returned by .read_h5ad_hdf5r() and
  # similar HDF5-driven loaders carry the AnnData /obs frame on a
  # `phenomapr_obs` attribute. Surface it as the cell-metadata table so
  # users get the exact same column set they would have gotten via the
  # Python reticulate/anndata path -- without needing Python.
  obs_attr <- attr(obj, "phenomapr_obs")
  if (is.data.frame(obs_attr) && nrow(obs_attr) > 0L) {
    md <- obs_attr
    if (is.null(rownames(md)) || any(is.na(rownames(md)))) {
      cn <- if (is.null(colnames(obj))) NULL else colnames(obj)
      if (!is.null(cn) && length(cn) == nrow(md)) rownames(md) <- cn
    }
    md$.cell_id <- rownames(md) %||% (
      if (!is.null(colnames(obj))) colnames(obj) else NULL
    )
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
# obsm key conventions for AnnData cell-segmentation centroids. The first
# match wins. We accept the bare name plus the Scanpy / Squidpy "X_*"
# variants so users don't have to rename their AnnData obsm just to get
# the segmentation embedding into the dropdown.
.anndata_segmentation_obsm_keys <- c(
  "segmentation",
  "X_segmentation",
  "X_segmentation_centroid",
  "segmentation_centroid",
  "segmentation_centroids"
)

# CosMx / imaging platforms often store slide-level coordinates in
# obsm["spatial"] while the per-cell segmentation centroids in FOV pixel
# space live in obs (e.g. x_FOV_px / y_FOV_px). Return the first matching
# column pair, or NULL when none are present.
.cosmx_segmentation_coord_columns <- function(meta_df) {
  if (!is.data.frame(meta_df) || !ncol(meta_df)) return(NULL)
  pairs <- list(
    c("x_FOV_px", "y_FOV_px"),
    c("x_fov_px", "y_fov_px"),
    c("CenterX_local_px", "CenterY_local_px"),
    c("centroid_x", "centroid_y"),
    c("x_centroid", "y_centroid")
  )
  for (p in pairs) {
    if (!all(p %in% colnames(meta_df))) next
    x <- suppressWarnings(as.numeric(meta_df[[p[1L]]]))
    y <- suppressWarnings(as.numeric(meta_df[[p[2L]]]))
    if (sum(is.finite(x) & is.finite(y)) >= 2L) {
      return(list(dim1_col = p[1L], dim2_col = p[2L]))
    }
  }
  NULL
}

# Build a (cells x 2) segmentation-centroid matrix from obs columns.
.build_segmentation_obsm_from_obs <- function(meta_df, cell_ids = NULL) {
  cols <- .cosmx_segmentation_coord_columns(meta_df)
  if (is.null(cols)) return(NULL)
  x <- suppressWarnings(as.numeric(meta_df[[cols$dim1_col]]))
  y <- suppressWarnings(as.numeric(meta_df[[cols$dim2_col]]))
  n <- length(x)
  if (n < 2L || sum(is.finite(x) & is.finite(y)) < 2L) return(NULL)
  m <- cbind(x = x, y = y)
  colnames(m) <- c(cols$dim1_col, cols$dim2_col)
  ids <- cell_ids
  if (is.null(ids) && !is.null(rownames(meta_df)) &&
      length(rownames(meta_df)) == n) {
    ids <- rownames(meta_df)
  }
  if (!is.null(ids) && length(ids) == n) rownames(m) <- as.character(ids)
  m
}

# Segmentation is available from obsm keys and/or CosMx-style obs columns.
.segmentation_embedding_available <- function(obsm_keys, meta_df = NULL) {
  if (length(intersect(.anndata_segmentation_obsm_keys, obsm_keys))) {
    return(TRUE)
  }
  if (is.data.frame(meta_df) && !is.null(.cosmx_segmentation_coord_columns(meta_df))) {
    return(TRUE)
  }
  FALSE
}

# Column names used to tag tissue sections / FOVs on spatial embeddings.
.spatial_sample_column_candidates <- function() {
  c("library_id", "sample_id", "sample", "slide_id", "section", "section_id",
    "core", "core_id", "fov", "fov_id", "fov_uid", "FOV")
}

.spatial_sample_vec_from_obs <- function(obs_df) {
  if (!is.data.frame(obs_df) || !ncol(obs_df)) return(NULL)
  hit <- intersect(.spatial_sample_column_candidates(), colnames(obs_df))
  if (!length(hit)) return(NULL)
  as.character(obs_df[[hit[1L]]])
}

# True iff the Seurat image carries a non-empty Seurat 5 Boundaries() set
# (Xenium / CosMx / MERSCOPE / Visium HD FOV objects). VisiumV1 spot
# grids and any pre-Seurat-5 image classes return character(0), so this
# test cleanly distinguishes "real cell-segmentation" from the legacy
# spot-lattice we already surface as the "spatial" reduction.
.has_segmentation_boundaries <- function(img) {
  b <- tryCatch(Seurat::Boundaries(img), error = function(e) character(0))
  length(b) > 0L
}

list_available_embeddings <- function(obj) {
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    nms <- names(obj@reductions)
    imgs <- methods::slot(obj, "images") %||% list()
    # Surface tissue / spot coordinates from spatial Seurat objects (Visium
    # etc.) as a synthetic "spatial" reduction so the Visualization tab can
    # overlay PhenoMapR scores on the tissue without any extra UI plumbing.
    if (length(imgs) > 0L) {
      nms <- c(nms, "spatial")
      # Add a sibling "segmentation" reduction whenever any image slot
      # carries cell-level segmentation boundaries (Seurat 5 FOV: Xenium,
      # CosMx, MERSCOPE, Visium HD). We keep "spatial" and "segmentation"
      # as separate dropdown entries because they answer different
      # questions: "spatial" plots whatever the upstream pipeline calls
      # the canonical spot/cell xy (often the spot grid or registered
      # centroid) while "segmentation" plots the cell-segmentation
      # centroids directly out of @boundaries$centroids -- this is the
      # only embedding that lines up 1:1 with the segmented cells in the
      # tissue image, so it's the right view for inspecting per-cell
      # phenotype scores in their morphological context.
      if (any(vapply(imgs, .has_segmentation_boundaries, logical(1L)))) {
        nms <- c(nms, "segmentation")
      }
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
    keys <- .collapse_anndata_embedding_keys(keys)
    if (!"segmentation" %in% keys) {
      md <- extract_object_metadata(obj)
      if (.segmentation_embedding_available(character(0), md)) {
        keys <- unique(c(keys, "segmentation"))
      }
    }
    return(keys)
  }
  # Plain matrix / Matrix / dgCMatrix: typically the result of an h5ad
  # loaded via .read_h5ad_hdf5r(). The obsm entries live on the
  # `phenomapr_obsm` attribute -- expose them with the same key
  # transformations as the live-AnnData branch above, so a "spatial"
  # entry stays "spatial" and "X_segmentation" / "X_segmentation_centroid"
  # collapse to a single "segmentation" choice.
  if (is.matrix(obj) || inherits(obj, "Matrix")) {
    obsm_attr <- attr(obj, "phenomapr_obsm")
    keys <- if (is.list(obsm_attr) && length(obsm_attr)) {
      .collapse_anndata_embedding_keys(names(obsm_attr))
    } else {
      character(0)
    }
    if (!"segmentation" %in% keys) {
      obs_attr <- attr(obj, "phenomapr_obs")
      if (.segmentation_embedding_available(names(obsm_attr %||% list()),
                                            obs_attr)) {
        keys <- unique(c(keys, "segmentation"))
      }
    }
    return(keys)
  }
  character(0)
}

# Both the live-AnnData branch and the hdf5r-loaded-matrix branch want
# the same dropdown-friendly transforms applied to the raw obsm keys:
# bring "spatial" to the end of the list, and collapse the various
# "X_segmentation" variants into a single "segmentation" entry. Pulled
# out into a helper so both branches stay in lock-step.
.collapse_anndata_embedding_keys <- function(keys) {
  if (!length(keys)) return(character(0))
  if ("spatial" %in% keys) {
    keys <- unique(c(setdiff(keys, "spatial"), "spatial"))
  }
  seg_hit <- intersect(.anndata_segmentation_obsm_keys, keys)
  if (length(seg_hit)) {
    keys <- unique(c(setdiff(keys, .anndata_segmentation_obsm_keys),
                     "segmentation"))
  }
  keys
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
  # Treat "segmentation" as a sibling spatial-style embedding: same
  # plotting frame (tissue coords, fixed aspect ratio, reversed Y), same
  # per-section sample tagging. The only difference is *which* per-cell
  # xy table we read off the object.
  is_segmentation <- identical(name, "segmentation")
  sample_vec <- NULL   # length = nrow(emb), labels each spot's section
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    if (is_segmentation) {
      # Walk every image slot and harvest per-cell centroids straight
      # out of @boundaries$centroids (Xenium / CosMx / MERSCOPE) or the
      # `centroids` boundary on Visium HD FOV objects. We use
      # `which = "centroids"` rather than `which = "segmentation"`
      # because the latter returns polygon vertices (one row per
      # vertex per cell) and we want one row per cell. If a particular
      # image only has the polygon boundary, we collapse polygons to
      # their per-cell median coordinate so the dropdown still
      # produces a usable per-cell embedding.
      imgs <- methods::slot(obj, "images") %||% list()
      img_names <- names(imgs)
      if (is.null(img_names) || any(!nzchar(img_names))) {
        img_names <- ifelse(nzchar(img_names %||% ""),
                            img_names,
                            paste0("image_", seq_along(imgs)))
      }
      parts <- lapply(seq_along(imgs), function(i) {
        img <- imgs[[i]]
        if (!.has_segmentation_boundaries(img)) return(NULL)
        bnames <- tryCatch(Seurat::Boundaries(img), error = function(e) character(0))
        coords <- NULL
        if ("centroids" %in% bnames) {
          coords <- tryCatch(
            Seurat::GetTissueCoordinates(obj, image = img_names[i],
                                          which = "centroids"),
            error = function(e) NULL
          )
        }
        if (is.null(coords) && length(bnames)) {
          # Fall back to polygon segmentation: median xy per cell so we
          # still get a per-cell point cloud for the dropdown.
          poly <- tryCatch(
            Seurat::GetTissueCoordinates(obj, image = img_names[i],
                                          which = bnames[1L]),
            error = function(e) NULL
          )
          if (!is.null(poly) && all(c("x", "y", "cell") %in% colnames(poly))) {
            coords <- stats::aggregate(
              cbind(x = poly$x, y = poly$y) ~ cell, data = poly,
              FUN = function(z) stats::median(as.numeric(z))
            )
            coords <- data.frame(
              x    = as.numeric(coords$x),
              y    = as.numeric(coords$y),
              cell = as.character(coords$cell),
              stringsAsFactors = FALSE
            )
          }
        }
        if (is.null(coords)) return(NULL)
        if (!all(c("x", "y") %in% colnames(coords))) return(NULL)
        m <- as.matrix(coords[, c("x", "y"), drop = FALSE])
        if ("cell" %in% colnames(coords)) {
          rownames(m) <- as.character(coords$cell)
        }
        list(emb = m, sample = rep(img_names[i], nrow(m)))
      })
      parts <- parts[!vapply(parts, is.null, logical(1L))]
      if (length(parts)) {
        emb <- do.call(rbind, lapply(parts, `[[`, "emb"))
        sample_vec <- unlist(lapply(parts, `[[`, "sample"), use.names = FALSE)
      }
    } else if (is_spatial) {
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
    if (is_segmentation) {
      # The dropdown collapses several conventional obsm keys
      # (segmentation, X_segmentation, X_segmentation_centroid, ...)
      # into a single "segmentation" entry. Resolve back to the actual
      # obsm key by trying them in order; first hit wins.
      keys <- tryCatch(PhenoMapR:::.anndata_obsm_keys(obj),
                       error = function(e) character(0))
      seg_hit <- intersect(.anndata_segmentation_obsm_keys, keys)
      if (length(seg_hit)) {
        emb <- tryCatch(PhenoMapR:::.anndata_obsm_array(obj, seg_hit[1L]),
                        error = function(e) NULL)
      }
      if (is.null(emb)) {
        md <- extract_object_metadata(obj)
        if (!is.null(md)) {
          emb <- .build_segmentation_obsm_from_obs(md, cell_ids = md$.cell_id)
        }
      }
    } else {
      emb <- tryCatch(PhenoMapR:::.anndata_obsm_array(obj, name),
                      error = function(e) NULL)
    }
    if (is_spatial || is_segmentation) {
      # AnnData multi-section spatial datasets keep one big obsm["spatial"]
      # (and similarly one obsm["segmentation"]) but disambiguate sections
      # via an obs column (Squidpy convention is `library_id`; Scanpy
      # users sometimes use `sample_id` / `sample` / `slide_id` /
      # `section`). We sniff for any of those and ride along on the
      # resulting vector so the section switcher works for segmentation
      # the same way it does for spatial.
      sample_vec <- tryCatch(
        .anndata_spatial_sample_vec(obj),
        error = function(e) NULL
      )
    }
  } else if (is.matrix(obj) || inherits(obj, "Matrix")) {
    # h5ads loaded via .read_h5ad_hdf5r() return a plain dgCMatrix
    # (genes x cells) with the original /obsm group stashed on the
    # `phenomapr_obsm` attribute. We resolve `name` against that list
    # the same way we do for live AnnData, including the segmentation
    # alias collapsing.
    obsm_attr <- attr(obj, "phenomapr_obsm")
    if (is.list(obsm_attr) && length(obsm_attr)) {
      pick <- if (is_segmentation) {
        seg_hit <- intersect(.anndata_segmentation_obsm_keys, names(obsm_attr))
        if (length(seg_hit)) seg_hit[1L] else NA_character_
      } else if (name %in% names(obsm_attr)) {
        name
      } else {
        NA_character_
      }
      if (!is.na(pick)) {
        emb <- obsm_attr[[pick]]
      } else if (is_segmentation) {
        obs_attr <- attr(obj, "phenomapr_obs")
        if (is.data.frame(obs_attr)) {
          emb <- .build_segmentation_obsm_from_obs(obs_attr,
                                                    cell_ids = colnames(obj))
        }
      }
    }
    if (is_spatial || is_segmentation) {
      # The hdf5r loader doesn't currently surface a per-cell tissue
      # section vector, but the obs data.frame is on `phenomapr_obs` -
      # sniff the same column conventions the AnnData branch uses so
      # multi-library spatial h5ads still get the section switcher.
      obs_attr <- attr(obj, "phenomapr_obs")
      if (is.data.frame(obs_attr)) {
        sample_vec <- .spatial_sample_vec_from_obs(obs_attr)
      }
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
    sample_vec <- rep(
      if (is_spatial || is_segmentation) "section_1" else NA_character_,
      nrow(emb)
    )
  } else {
    # Promote anything non-character (factor / pandas categorical) to
    # plain character so the select input choices stay clean.
    sample_vec <- as.character(sample_vec)
    sample_vec[is.na(sample_vec) | !nzchar(sample_vec)] <- "(unknown)"
  }
  # Segmentation centroids live in the same tissue coordinate frame as
  # the "spatial" reduction, so the Visualization tab needs the same
  # coord_fixed() + reversed-Y layout. We therefore mark segmentation
  # rows as `is_spatial = TRUE`. The original embedding name is
  # preserved verbatim (dim1_name / dim2_name) so plot titles and
  # downloads still differentiate "spatial" from "segmentation".
  data.frame(
    cell_id = as.character(ids),
    dim1 = as.numeric(emb[, 1L]),
    dim2 = as.numeric(emb[, 2L]),
    dim1_name = d1_name,
    dim2_name = d2_name,
    is_spatial = is_spatial || is_segmentation,
    sample = sample_vec,
    stringsAsFactors = FALSE
  )
}

# ---- spatial polygon / cell-mask extraction ---------------------------------
#
# Seurat 5 FOV objects (Xenium / CosMx / MERSCOPE / Visium HD) store
# per-cell segmentation boundaries in @images[[*]]@boundaries. The
# Visualization tab can render those as filled polygons when the user
# selects the "segmentation" reduction. Plain matrices / h5ad loads can
# also carry a long-format `phenomapr_polygons` attribute:
#   data.frame(cell_id, x, y, sample) with multiple rows per cell.

.seurat_polygon_boundary_name <- function(bnames) {
  bnames <- as.character(bnames)
  bnames <- bnames[nzchar(bnames)]
  if (!length(bnames)) return(NA_character_)
  poly_names <- setdiff(bnames, "centroids")
  if (!length(poly_names)) return(NA_character_)
  if ("segmentation" %in% poly_names) return("segmentation")
  poly_names[1L]
}

.seurat_has_polygon_boundaries <- function(img) {
  bnames <- tryCatch(Seurat::Boundaries(img), error = function(e) character(0))
  !is.na(.seurat_polygon_boundary_name(bnames))
}

.normalize_spatial_polygon_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  req_cols <- c("cell_id", "x", "y")
  if (!all(req_cols %in% colnames(df))) return(NULL)
  out <- data.frame(
    cell_id = as.character(df$cell_id),
    x = suppressWarnings(as.numeric(df$x)),
    y = suppressWarnings(as.numeric(df$y)),
    stringsAsFactors = FALSE
  )
  if ("sample" %in% colnames(df)) {
    out$sample <- as.character(df$sample)
  } else {
    out$sample <- "section_1"
  }
  ok <- is.finite(out$x) & is.finite(out$y) & !is.na(out$cell_id) & nzchar(out$cell_id)
  out <- out[ok, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  # Keep only cells with at least three vertices (a degenerate "polygon").
  n_vert <- as.integer(stats::ave(seq_len(nrow(out)), out$cell_id, FUN = length))
  out <- out[n_vert >= 3L, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out
}

.filter_spatial_polygons_section <- function(poly_df, section = NULL) {
  if (is.null(poly_df) || !nrow(poly_df)) return(poly_df)
  if (is.null(section) || !nzchar(section) ||
      identical(section, "__all__") || !"sample" %in% colnames(poly_df)) {
    return(poly_df)
  }
  poly_df[!is.na(poly_df$sample) & poly_df$sample == section, , drop = FALSE]
}

.extract_seurat_spatial_polygons <- function(obj, section = NULL) {
  if (!inherits(obj, "Seurat") || !requireNamespace("Seurat", quietly = TRUE)) {
    return(NULL)
  }
  imgs <- methods::slot(obj, "images") %||% list()
  if (!length(imgs)) return(NULL)
  img_names <- names(imgs)
  if (is.null(img_names) || any(!nzchar(img_names))) {
    img_names <- paste0("image_", seq_along(imgs))
  }
  parts <- lapply(seq_along(imgs), function(i) {
    if (!is.null(section) && nzchar(section) && !identical(section, "__all__") &&
        img_names[i] != section) {
      return(NULL)
    }
    img <- imgs[[i]]
    if (!.seurat_has_polygon_boundaries(img)) return(NULL)
    bnames <- tryCatch(Seurat::Boundaries(img), error = function(e) character(0))
    which_use <- .seurat_polygon_boundary_name(bnames)
    if (is.na(which_use)) return(NULL)
    coords <- tryCatch(
      Seurat::GetTissueCoordinates(obj, image = img_names[i], which = which_use),
      error = function(e) NULL
    )
    if (is.null(coords) || !nrow(coords)) return(NULL)
    if (!all(c("x", "y", "cell") %in% colnames(coords))) return(NULL)
    data.frame(
      cell_id = as.character(coords$cell),
      x = as.numeric(coords$x),
      y = as.numeric(coords$y),
      sample = img_names[i],
      stringsAsFactors = FALSE
    )
  })
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (!length(parts)) return(NULL)
  .normalize_spatial_polygon_df(do.call(rbind, parts))
}

#' @keywords internal
spatial_polygons_available <- function(obj) {
  if (is.null(obj)) return(FALSE)
  poly <- extract_spatial_polygons(obj, section = NULL)
  !is.null(poly) && nrow(poly) > 0L
}

#' @keywords internal
extract_spatial_polygons <- function(obj, section = NULL) {
  poly <- NULL
  if (inherits(obj, "Seurat")) {
    poly <- .extract_seurat_spatial_polygons(obj, section = section)
  } else if (is.matrix(obj) || inherits(obj, "Matrix")) {
    attr_poly <- attr(obj, "phenomapr_polygons")
    if (is.data.frame(attr_poly)) {
      poly <- .normalize_spatial_polygon_df(attr_poly)
      poly <- .filter_spatial_polygons_section(poly, section)
    }
  } else if (inherits(obj, "python.builtin.object")) {
    poly <- .extract_anndata_spatial_polygons(obj, section = section)
  }
  if (is.null(poly) || !nrow(poly)) return(NULL)
  poly
}

.extract_anndata_spatial_polygons <- function(obj, section = NULL) {
  if (!inherits(obj, "python.builtin.object")) return(NULL)
  NULL
}

#' Join per-vertex polygon rows with per-cell coloring columns for ggplot.
#' @keywords internal
join_spatial_polygon_colors <- function(poly_df, cell_df, value_col) {
  if (is.null(poly_df) || !nrow(poly_df) || is.null(cell_df) ||
      !value_col %in% colnames(cell_df)) {
    return(NULL)
  }
  keep <- unique(c("cell_id", value_col))
  meta <- cell_df[, keep, drop = FALSE]
  out <- dplyr::left_join(poly_df, meta, by = "cell_id")
  if (!value_col %in% colnames(out)) return(NULL)
  out
}

# Per-FOV pixel coordinates (e.g. CosMx x_FOV_px / y_FOV_px) are local to each
# imaging field. For "All sections (combined)" on the segmentation reduction,
# plot in the same global slide frame as the spatial reduction.
#' @keywords internal
.spatial_coords_are_fov_local <- function(plot_df) {
  if (is.null(plot_df) || !nrow(plot_df)) return(FALSE)
  n1 <- unique(plot_df$dim1_name)[1L] %||% ""
  n2 <- unique(plot_df$dim2_name)[1L] %||% ""
  grepl("FOV_px|fov_px", paste(n1, n2, collapse = " "), ignore.case = TRUE)
}

#' Swap FOV-local embedding coordinates for global spatial coordinates.
#'
#' Used when the user selects segmentation + all sections combined on CosMx-
#' style data: each FOV shares a 0..N pixel range, so we re-use the object's
#' global \code{spatial} embedding positions (slide-mm / tissue frame) while
#' keeping the segmentation reduction selected.
#'
#' @param obj Expression object passed to \code{extract_embedding()}.
#' @param emb_df Embedding data.frame from the segmentation reduction.
#' @return \code{emb_df} with \code{dim1}/\code{dim2} replaced when possible.
#' @keywords internal
apply_global_spatial_coords_for_combined <- function(obj, emb_df) {
  if (is.null(emb_df) || !nrow(emb_df)) return(emb_df)
  if (!isTRUE(.spatial_coords_are_fov_local(emb_df))) return(emb_df)

  global_emb <- tryCatch(extract_embedding(obj, "spatial"), error = function(e) NULL)
  if (is.null(global_emb) || !nrow(global_emb)) return(emb_df)

  idx <- match(emb_df$cell_id, global_emb$cell_id)
  if (!any(is.finite(idx))) return(emb_df)

  out <- emb_df
  out$dim1 <- global_emb$dim1[idx]
  out$dim2 <- global_emb$dim2[idx]
  out$dim1_name <- global_emb$dim1_name[1L] %||% out$dim1_name
  out$dim2_name <- global_emb$dim2_name[1L] %||% out$dim2_name
  if ("sample" %in% colnames(global_emb)) {
    out$sample <- global_emb$sample[idx]
  }
  out
}

# Probe an AnnData object's `obs` for the column conventionally used to
# disambiguate tissue sections in a multi-sample spatial dataset.
# Tries (in order): library_id, sample_id, sample, slide_id, section.
# Returns a character vector the same length as adata.n_obs, or NULL.
.anndata_spatial_sample_vec <- function(obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  candidates <- .spatial_sample_column_candidates()
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
    # Use the Seurat 4 / 5 compatibility shim so we don't trip the
    # "slot is now defunct" error on SeuratObject >= 5.0 (the shim
    # picks layer= or slot= based on what the installed
    # Seurat::GetAssayData actually accepts).
    m <- tryCatch(
      PhenoMapR:::.get_assay_data_compat(obj, assay = use_assay, slot = layer_name),
      error = function(e) NULL
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
#
# The plot-radius slider in the navbar populates the global option
# `phenomapr.plot_radius_pt` on change (see app.R server). All wrappers
# read that as their default `radius`, so a single slider drives every
# rounded-geom call in the app without having to thread the value down
# manually at each site. Callers can still override per-call by passing
# `radius = grid::unit(...)` explicitly.
.use_chicklet <- function() {
  requireNamespace("ggchicklet2", quietly = TRUE)
}

# Default radius (in points) for the bar / stack / boxplot wrappers.
# Falls back to 3 pt when nothing has set the option (e.g., in unit
# tests sourcing helpers.R in isolation). On-screen plots are built
# once at the default radius; per-export tweaks are layered on at
# preview / save time via .apply_chicklet_radius() below.
.plot_radius_unit <- function(fallback_pt = 3) {
  r <- getOption("phenomapr.plot_radius_pt", default = fallback_pt)
  if (is.null(r) || !is.numeric(r) || !is.finite(r) || r < 0) {
    r <- fallback_pt
  }
  grid::unit(r, "pt")
}

# Walk a ggplot's layers and replace the `radius` parameter on every
# ggchicklet2 geom in-place (these store it in `geom_params$radius`
# as a `grid::unit`). Used by the plot-download modal so users can
# tweak the corner sharpness per export without re-running the
# underlying renderPlot. Returns the (possibly modified) plot.
#
# Safe to call on any ggplot:
#   - Non-ggplot inputs are returned unchanged.
#   - Plots without ggchicklet2 layers are no-ops.
#   - Invalid radii (negative, NA, non-numeric) are ignored.
.apply_chicklet_radius <- function(p, radius_pt) {
  if (!inherits(p, c("ggplot", "patchwork"))) return(p)
  if (!is.numeric(radius_pt) || length(radius_pt) != 1L ||
      !is.finite(radius_pt) || radius_pt < 0) {
    return(p)
  }
  new_unit <- grid::unit(radius_pt, "pt")

  apply_to_one <- function(gg) {
    if (!inherits(gg, "ggplot")) return(gg)
    if (!length(gg$layers)) return(gg)
    for (i in seq_along(gg$layers)) {
      layer <- gg$layers[[i]]
      if (!is.null(layer$geom_params) &&
          "radius" %in% names(layer$geom_params)) {
        layer$geom_params$radius <- new_unit
      }
      if (!is.null(layer$stat_params) &&
          "radius" %in% names(layer$stat_params)) {
        layer$stat_params$radius <- new_unit
      }
    }
    gg
  }

  if (inherits(p, "patchwork")) {
    # patchwork stores child plots in p$patches$plots and the
    # top-level plot itself on p; walk both.
    p <- apply_to_one(p)
    plots <- tryCatch(p$patches$plots, error = function(e) NULL)
    if (length(plots)) {
      p$patches$plots <- lapply(plots, apply_to_one)
    }
    return(p)
  }
  apply_to_one(p)
}

# Rounded `geom_col()` — pre-aggregated heights (stat = "identity").
.geom_rounded_col <- function(..., radius = .plot_radius_unit(3),
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
.geom_rounded_stack <- function(..., radius = .plot_radius_unit(3),
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
.geom_rounded_histogram <- function(..., radius = .plot_radius_unit(2),
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
.geom_rounded_boxplot <- function(..., radius = .plot_radius_unit(4)) {
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

# ---- Marker heatmap color schemes (Shiny + plot_phenotype_markers) ---------

#' @keywords internal
marker_heatmap_color_mode_choices <- function() {
  c(
    "Package default" = "default",
    "ColorBrewer" = "brewer",
    "Viridis" = "viridis",
    "Custom colors" = "manual"
  )
}

#' @keywords internal
marker_heatmap_brewer_choices <- function(kind = c("diverging", "qualitative", "both")) {
  kind <- match.arg(kind)
  pal <- PhenoMapR::list_marker_heatmap_color_palettes()
  names <- switch(
    kind,
    diverging = pal$brewer$diverging,
    qualitative = pal$brewer$qualitative,
    both = c(pal$brewer$diverging, pal$brewer$qualitative)
  )
  if (!length(names)) return(c("(RColorBrewer unavailable)" = ""))
  stats::setNames(names, names)
}

#' @keywords internal
marker_heatmap_viridis_choices <- function() {
  pal <- PhenoMapR::list_marker_heatmap_color_palettes()
  if (!length(pal$viridis)) return(c("(viridis unavailable)" = ""))
  stats::setNames(pal$viridis, pal$viridis)
}

#' @keywords internal
.marker_heatmap_resolve_color_spec <- function(input, component, manual_colors_fn) {
  mode <- input[[paste0("hm_colors_", component, "_mode")]] %||% "default"
  if (identical(mode, "default")) return("default")
  if (identical(mode, "brewer")) {
    name <- input[[paste0("hm_colors_", component, "_brewer")]] %||% ""
    if (!nzchar(name)) return("default")
    return(paste0("brewer:", name))
  }
  if (identical(mode, "viridis")) {
    name <- input[[paste0("hm_colors_", component, "_viridis")]] %||% ""
    if (!nzchar(name)) return("default")
    return(paste0("viridis:", name))
  }
  list(source = "manual", colors = manual_colors_fn(input))
}

#' Build `color_schemes` list for plot_phenotype_markers() from Shiny inputs.
#' @keywords internal
marker_heatmap_color_schemes_from_input <- function(input) {
  list(
    phenotype = .marker_heatmap_resolve_color_spec(input, "phenotype", function(inp) {
      c(
        inp$hm_colors_phenotype_fav %||% "#2166AC",
        inp$hm_colors_phenotype_other %||% "#F7F7F7",
        inp$hm_colors_phenotype_adv %||% "#B2182B"
      )
    }),
    score = .marker_heatmap_resolve_color_spec(input, "score", function(inp) {
      c(
        inp$hm_colors_score_low %||% "#2166AC",
        inp$hm_colors_score_mid %||% "#FFFFFF",
        inp$hm_colors_score_high %||% "#B2182B"
      )
    }),
    expression = .marker_heatmap_resolve_color_spec(input, "expr", function(inp) {
      c(
        inp$hm_colors_expr_low %||% "#1A1A1A",
        inp$hm_colors_expr_mid %||% "#FFFFFF",
        inp$hm_colors_expr_high %||% "#67001F"
      )
    }),
    celltype = .marker_heatmap_resolve_color_spec(input, "celltype", function(inp) {
      txt <- inp$hm_colors_celltype_manual %||% ""
      cols <- trimws(unlist(strsplit(txt, "[,;\\s]+")))
      cols <- cols[nzchar(cols)]
      if (!length(cols)) {
        cols <- c("#2A9D8F", "#E76F51", "#264653", "#6A4C93", "#3A86FF")
      }
      cols
    })
  )
}

#' UI builder for one marker-heatmap color component (mode + palette pickers).
#' @keywords internal
marker_heatmap_color_controls <- function(component,
                                          label,
                                          brewer_kind = c("diverging", "qualitative", "both"),
                                          default_brewer = "RdBu",
                                          default_viridis = "viridis",
                                          manual_ui = NULL) {
  brewer_kind <- match.arg(brewer_kind)
  brewer_choices <- marker_heatmap_brewer_choices(brewer_kind)
  viridis_choices <- marker_heatmap_viridis_choices()
  tagList(
    selectInput(
      paste0("hm_colors_", component, "_mode"),
      label,
      choices = marker_heatmap_color_mode_choices(),
      selected = "default"
    ),
    conditionalPanel(
      sprintf("input.hm_colors_%s_mode == 'brewer'", component),
      selectInput(
        paste0("hm_colors_", component, "_brewer"),
        paste(label, "(ColorBrewer palette)"),
        choices = brewer_choices,
        selected = if (default_brewer %in% brewer_choices) default_brewer else brewer_choices[1L]
      )
    ),
    conditionalPanel(
      sprintf("input.hm_colors_%s_mode == 'viridis'", component),
      selectInput(
        paste0("hm_colors_", component, "_viridis"),
        paste(label, "(viridis palette)"),
        choices = viridis_choices,
        selected = if (default_viridis %in% viridis_choices) default_viridis else viridis_choices[1L]
      )
    ),
    conditionalPanel(
      sprintf("input.hm_colors_%s_mode == 'manual'", component),
      manual_ui
    )
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

# ---- single global ANOVA across cell types -------------------------------
#
# Used by the "Score by cell type" plot when no source column is
# mapped. Replaces the previous "pairwise Wilcoxon between every cell
# type pair" annotation (which produced a forest of brackets) with a
# single global F-test answering: "Does score vary across cell types
# at all?". Returns a data.frame with one row carrying the
# bracket-style geometry (xmin / xmax spanning all cell types, y_pos
# slightly above the data) plus the F-test p-value and a Wilkinson
# significance label.
celltype_anova_pvalue <- function(df,
                                  score_col = "Score",
                                  cell_type_col = "Cell_type",
                                  cell_levels = NULL) {
  req_cols <- c(score_col, cell_type_col)
  if (!all(req_cols %in% colnames(df))) return(NULL)
  d <- df[, req_cols]
  colnames(d) <- c("Score", "Cell_type")
  d <- d[is.finite(d$Score) & !is.na(d$Cell_type), ]
  if (nrow(d) < 4L) return(NULL)
  d$Score_scaled <- as.numeric(scale(d$Score))

  observed <- unique(as.character(d$Cell_type))
  if (is.null(cell_levels)) {
    cell_levels <- sort(observed)
  } else {
    cell_levels <- intersect(as.character(cell_levels), observed)
    if (!length(cell_levels)) cell_levels <- sort(observed)
  }
  if (length(cell_levels) < 2L) return(NULL)

  d$Cell_type <- factor(d$Cell_type, levels = cell_levels)
  fit <- tryCatch(
    stats::aov(Score_scaled ~ Cell_type, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm) || !length(sm)) return(NULL)
  p <- tryCatch(sm[[1L]][["Pr(>F)"]][1L], error = function(e) NA_real_)
  f <- tryCatch(sm[[1L]][["F value"]][1L], error = function(e) NA_real_)
  if (is.na(p)) return(NULL)

  y_top <- max(d$Score_scaled, na.rm = TRUE)
  if (!is.finite(y_top)) y_top <- 0

  # NOTE: the no-source ANOVA case is rendered as a centred text
  # annotation (not a ggsignif bracket) -- the renderer in app.R reads
  # `attr(out, "render_kind")` and skips drawing the connector line and
  # significance stars. We only emit the bare "ANOVA (p = ...)" label
  # so the plot stays clean.
  out <- data.frame(
    xmin   = 1L,
    xmax   = length(cell_levels),
    x_mid  = (1L + length(cell_levels)) / 2,
    y_pos  = y_top * 1.08,
    p_val  = p,
    F_val  = f,
    label  = sprintf("ANOVA (p = %s)", format.pval(p, digits = 2L)),
    stringsAsFactors = FALSE
  )
  attr(out, "cell_levels") <- cell_levels
  attr(out, "y_top") <- y_top
  attr(out, "render_kind") <- "annotation"
  out
}

# ---- per-cell-type ANOVA across 3+ sources -------------------------------
#
# Used by the "Score by cell type and source" plot when a source column
# *is* mapped and contains more than 2 levels. Runs a one-way ANOVA
# (Score ~ Source) *within* each cell type and returns one bracket
# row per cell type with the F-test p-value. Cell types with fewer
# than 2 sources observed are dropped.
celltype_source_anova_pvalues <- function(df,
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

  rows <- lapply(cell_levels, function(ct) {
    sub <- d[d$Cell_type == ct, , drop = FALSE]
    srcs <- unique(sub$Source)
    if (length(srcs) < 2L) return(NULL)
    sub$Source <- factor(sub$Source, levels = srcs)
    fit <- tryCatch(stats::aov(Score_scaled ~ Source, data = sub),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    if (is.null(sm) || !length(sm)) return(NULL)
    p <- tryCatch(sm[[1L]][["Pr(>F)"]][1L], error = function(e) NA_real_)
    if (is.na(p)) return(NULL)
    lab <- if (p < 0.001) "***" else if (p < 0.01) "**" else
           if (p < 0.05)  "*"   else "ns"
    data.frame(
      Cell_type = ct,
      p_val     = p,
      label     = lab,
      n_sources = length(srcs),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) return(NULL)
  if (isTRUE(significant_only)) {
    out <- out[out$p_val < p_threshold, , drop = FALSE]
    if (!nrow(out)) return(NULL)
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

#' Column name for phenotype groups tied to a score column
#' @keywords internal
phenotype_group_column_for <- function(score_column) {
  if (is.null(score_column) || !nzchar(score_column)) return(NA_character_)
  paste0("phenotype_group_", score_column)
}

#' Resolve phenotype group column from groups table + active score column
#' @keywords internal
resolve_phenotype_group_column <- function(groups, score_column = NULL) {
  if (is.null(groups)) return(NA_character_)
  target <- phenotype_group_column_for(score_column)
  if (!is.na(target) && target %in% colnames(groups)) return(target)
  grp_cols <- grep("^phenotype_group_", colnames(groups), value = TRUE)
  if (length(grp_cols)) grp_cols[1L] else NA_character_
}

#' Active score column from Shiny score-column selector
#' @keywords internal
active_score_column <- function(scores, score_col_input = NULL) {
  sc <- score_col_input %||% ""
  if (identical(sc, "__all__") || !nzchar(sc)) {
    if (!is.null(scores) && ncol(scores) >= 1L) return(colnames(scores)[1L])
    return("")
  }
  sc
}

#' @keywords internal
.precog_pancreatic_gene_z <- function() {
  ref <- tryCatch(PhenoMapR::precog, error = function(e) NULL)
  if (is.null(ref) || !"Pancreatic" %in% colnames(ref)) return(NULL)
  z <- ref[, "Pancreatic"]
  names(z) <- rownames(ref)
  z[is.finite(z)]
}

#' Catalog metadata for the Shiny demo source dataset (CRA001160).
#' @keywords internal
shiny_demo_source_info <- function() {
  list(
    accession = "CRA001160",
    title = "Pancreatic ductal adenocarcinoma (PAAD) scRNA-seq",
    cancer_type = "Pancreatic",
    n_cells_total = 57530L,
    n_patients = 34L,
    n_tumors = 25L,
    n_controls = 9L,
    paper_label = "Peng et al. 2019",
    paper_url = "https://doi.org/10.1038/s41422-019-0195-y",
    data_source = "TISCH2",
    data_source_url = "https://tisch.compbio.cn/home/",
    reference_signature = "PRECOG \u00b7 Pancreatic"
  )
}

.shiny_demo_pool_cache <- new.env(parent = emptyenv())

.shiny_first_existing <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  for (p in paths) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NULL
}

.shiny_demo_search_roots <- function() {
  helpers <- system.file("shiny", "helpers.R", package = "PhenoMapR")
  pkg_from_helpers <- if (nzchar(helpers)) {
    normalizePath(file.path(dirname(helpers), "..", ".."), winslash = "/")
  } else {
    character(0)
  }
  roots <- c(
    getwd(),
    Sys.getenv("PHENOMAPR_PACKAGE_ROOT", unset = ""),
    system.file(package = "PhenoMapR"),
    pkg_from_helpers
  )
  unique(roots[nzchar(roots)])
}

.shiny_row_var <- function(mat) {
  if (inherits(mat, "Matrix")) {
    mu <- Matrix::rowMeans(mat)
    ex2 <- Matrix::rowMeans(mat^2)
    pmax(as.numeric(ex2 - mu^2), 0)
  } else {
    apply(mat, 1L, stats::var)
  }
}

.standardize_cra001160_metadata <- function(md) {
  n <- nrow(md)
  if (n == 0L) {
    return(data.frame(
      .cell_id = character(0),
      cell_type = character(0),
      cell_type_original = character(0),
      Source = character(0),
      Patient = character(0),
      UMAP_1 = numeric(0),
      UMAP_2 = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  .pick <- function(cands, kind = c("char", "num")) {
    kind <- match.arg(kind)
    hit <- cands[cands %in% names(md)]
    if (!length(hit)) {
      return(if (kind == "num") rep(NA_real_, n) else rep(NA_character_, n))
    }
    if (kind == "num") as.numeric(md[[hit[1L]]]) else as.character(md[[hit[1L]]])
  }
  cell_id <- if ("Cell" %in% names(md)) {
    as.character(md$Cell)
  } else {
    as.character(rownames(md))
  }
  if (length(cell_id) != n) cell_id <- as.character(rownames(md))
  major <- .pick(c(
    "Celltype (major-lineage)", "Celltype..major.lineage.",
    "cell_type_major", "Celltype", "seurat_clusters"
  ))
  minor <- .pick(c(
    "Celltype (minor-lineage)", "Celltype..minor.lineage.",
    "cell_type_minor", "Celltype (original)", "Celltype..original."
  ))
  original <- .pick(c(
    "Celltype (original)", "Celltype..original.",
    "Celltype (minor-lineage)", "Celltype..minor.lineage."
  ))
  minor <- gsub(" cell$", "", minor, ignore.case = TRUE)
  original <- gsub(" cell$", "", original, ignore.case = TRUE)
  cell_type <- minor
  na_ct <- is.na(cell_type) | !nzchar(cell_type)
  if (any(na_ct)) cell_type[na_ct] <- original[na_ct]
  na_ct <- is.na(cell_type) | !nzchar(cell_type)
  if (any(na_ct)) cell_type[na_ct] <- major[na_ct]
  data.frame(
    .cell_id = cell_id,
    cell_type = cell_type,
    cell_type_minor = minor,
    cell_type_major = major,
    cell_type_original = original,
    Source = .pick(c("Source", "Tumor_Normal", "Sample_Site")),
    Patient = .pick(c("Patient", "Sample", "orig.ident")),
    UMAP_1 = .pick(c("UMAP_1", "umap_1"), "num"),
    UMAP_2 = .pick(c("UMAP_2", "umap_2"), "num"),
    stringsAsFactors = FALSE
  )
}

#' Align metadata rows to matrix columns and attach as \code{phenomapr_coldata}.
#' @keywords internal
.attach_matrix_coldata <- function(expr, metadata, cell_id_col = ".cell_id") {
  if (is.null(expr) || is.null(metadata) || !nrow(metadata)) return(expr)
  cells <- colnames(expr)
  if (is.null(cells) || !length(cells)) return(expr)

  if (!is.null(cell_id_col) && nzchar(cell_id_col) &&
      cell_id_col %in% colnames(metadata)) {
    idx <- match(cells, metadata[[cell_id_col]])
  } else if (!is.null(rownames(metadata))) {
    idx <- match(cells, rownames(metadata))
  } else if (nrow(metadata) == length(cells)) {
    idx <- seq_along(cells)
  } else {
    return(expr)
  }
  if (anyNA(idx)) return(expr)

  coldata <- metadata[idx, , drop = FALSE]
  rownames(coldata) <- cells
  attr(expr, "phenomapr_coldata") <- coldata
  expr
}

#' Public wrapper for attaching matrix colData in the Shiny scoring path.
attach_matrix_coldata <- function(expr, metadata, cell_id_col = ".cell_id") {
  .attach_matrix_coldata(expr, metadata, cell_id_col = cell_id_col)
}

#' Attach metadata UMAP coordinates as \code{phenomapr_obsm} on a plain matrix.
#' @keywords internal
.attach_demo_matrix_obsm <- function(expr, metadata) {
  expr <- .attach_matrix_coldata(expr, metadata, cell_id_col = ".cell_id")
  if (is.null(metadata) || !all(c("UMAP_1", "UMAP_2") %in% colnames(metadata))) {
    return(expr)
  }
  cells <- colnames(expr)
  coldata <- attr(expr, "phenomapr_coldata")
  if (!is.null(coldata) && all(c("UMAP_1", "UMAP_2") %in% colnames(coldata))) {
    emb <- cbind(
      UMAP_1 = as.numeric(coldata$UMAP_1),
      UMAP_2 = as.numeric(coldata$UMAP_2)
    )
  } else {
    idx <- match(cells, metadata$.cell_id)
    if (anyNA(idx)) return(expr)
    emb <- cbind(
      UMAP_1 = as.numeric(metadata$UMAP_1[idx]),
      UMAP_2 = as.numeric(metadata$UMAP_2[idx])
    )
  }
  rownames(emb) <- cells
  attr(expr, "phenomapr_obsm") <- list(UMAP = emb)
  expr
}

.shiny_demo_bundle_path <- function() {
  candidates <- character(0)
  sf <- system.file("extdata", "shiny", "PAAD_CRA001160_demo_5000.rds", package = "PhenoMapR")
  if (nzchar(sf)) candidates <- c(candidates, sf)
  shiny_helpers <- system.file("shiny", "helpers.R", package = "PhenoMapR")
  if (nzchar(shiny_helpers)) {
    candidates <- c(
      candidates,
      normalizePath(
        file.path(dirname(shiny_helpers), "..", "extdata", "shiny", "PAAD_CRA001160_demo_5000.rds"),
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
  for (r in .shiny_demo_search_roots()) {
    candidates <- c(
      candidates,
      file.path(r, "inst", "extdata", "shiny", "PAAD_CRA001160_demo_5000.rds")
    )
  }
  candidates <- c(candidates, Sys.getenv("PHENOMAPR_SHINY_DEMO_RDS", unset = ""))
  .shiny_first_existing(unique(candidates))
}

.load_shiny_demo_pool_from_seurat <- function(path) {
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  obj <- tryCatch(Seurat::UpdateSeuratObject(obj), error = function(e) obj)
  assay <- if ("RNA" %in% names(obj@assays)) "RNA" else names(obj@assays)[1L]
  layer <- if ("data" %in% .seurat_assay_layers(obj, assay)) {
    "data"
  } else if ("counts" %in% .seurat_assay_layers(obj, assay)) {
    "counts"
  } else {
    "data"
  }
  expr <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(expr)) return(NULL)
  cells <- colnames(expr)
  md <- .standardize_cra001160_metadata(obj@meta.data)
  if (any(is.na(md$UMAP_1)) || any(is.na(md$UMAP_2))) {
    red_name <- if ("umap" %in% names(obj@reductions)) {
      "umap"
    } else if (length(obj@reductions)) {
      names(obj@reductions)[1L]
    } else {
      NA_character_
    }
    if (!is.na(red_name)) {
      emb <- tryCatch(Seurat::Embeddings(obj, red_name), error = function(e) NULL)
      if (!is.null(emb) && ncol(emb) >= 2L) {
        md$UMAP_1 <- emb[cells, 1L]
        md$UMAP_2 <- emb[cells, 2L]
      }
    }
  }
  md$.cell_id <- cells
  md <- md[match(cells, md$.cell_id), , drop = FALSE]
  if (anyNA(md$.cell_id)) return(NULL)
  list(expression = expr, metadata = md, from = "CRA001160")
}

.load_shiny_demo_pool_from_h5_tsv <- function(h5_path, tsv_path) {
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  expr <- tryCatch(Seurat::Read10X_h5(h5_path), error = function(e) NULL)
  if (is.null(expr)) return(NULL)
  md_raw <- tryCatch(
    utils::read.delim(tsv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(md_raw)) return(NULL)
  md <- .standardize_cra001160_metadata(md_raw)
  cells <- colnames(expr)
  md <- md[match(cells, md$.cell_id), , drop = FALSE]
  if (anyNA(md$.cell_id)) return(NULL)
  list(expression = expr, metadata = md, from = "CRA001160")
}

#' Load (and cache) the full CRA001160 demo expression pool.
#' @keywords internal
load_shiny_demo_pool <- function() {
  if (!is.null(.shiny_demo_pool_cache$pool)) {
    return(.shiny_demo_pool_cache$pool)
  }
  roots <- .shiny_demo_search_roots()
  pool <- NULL

  h5_paths <- c(
    Sys.getenv("PHENOMAPR_CRA001160_H5", unset = ""),
    unlist(lapply(roots, function(r) {
      file.path(r, "vignettes", "PAAD_CRA001160_expression.h5")
    }), use.names = FALSE)
  )
  tsv_paths <- c(
    Sys.getenv("PHENOMAPR_CRA001160_META", unset = ""),
    unlist(lapply(roots, function(r) {
      file.path(r, "vignettes", "PAAD_CRA001160_CellMetainfo_table.tsv")
    }), use.names = FALSE)
  )
  h5_file <- .shiny_first_existing(h5_paths)
  tsv_file <- .shiny_first_existing(tsv_paths)
  if (!is.null(h5_file) && !is.null(tsv_file)) {
    pool <- .load_shiny_demo_pool_from_h5_tsv(h5_file, tsv_file)
  }

  if (is.null(pool)) {
    rds_paths <- unlist(lapply(roots, function(r) {
      c(
        file.path(r, "inst", "extdata", "vignette_subsets", "PAAD_CRA001160_seurat_subset.rds"),
        file.path(r, "inst", "extdata", "shiny", "PAAD_CRA001160_demo_pool.rds"),
        file.path(r, "vignette_subsets", "PAAD_CRA001160_seurat_subset.rds")
      )
    }), use.names = FALSE)
    rds_paths <- c(
      rds_paths,
      Sys.getenv("PHENOMAPR_CRA001160_RDS", unset = ""),
      system.file("extdata", "vignette_subsets", "PAAD_CRA001160_seurat_subset.rds",
                  package = "PhenoMapR"),
      system.file("extdata", "shiny", "PAAD_CRA001160_demo_pool.rds", package = "PhenoMapR")
    )
    rds_file <- .shiny_first_existing(rds_paths)
    if (!is.null(rds_file)) {
      pool <- .load_shiny_demo_pool_from_seurat(rds_file)
    }
  }

  if (is.null(pool)) {
    url <- Sys.getenv("PHENOMAPR_CRA001160_RDS_URL", unset = "")
    dest <- file.path(tempdir(), "PAAD_CRA001160_seurat_subset.rds")
    dl_ok <- FALSE
    if (nzchar(url)) {
      dl_ok <- tryCatch({
        utils::download.file(url, dest, quiet = TRUE, mode = "wb")
        TRUE
      }, error = function(e) FALSE)
    }
    if (!dl_ok && requireNamespace("googledrive", quietly = TRUE)) {
      dl_ok <- tryCatch({
        googledrive::drive_deauth()
        googledrive::drive_download(
          googledrive::as_id("14p_fYIFeuuRdXBF3J-5ZsXElq_mduSzb"),
          path = dest,
          overwrite = TRUE
        )
        TRUE
      }, error = function(e) FALSE)
    }
    if (isTRUE(dl_ok) && file.exists(dest)) {
      pool <- .load_shiny_demo_pool_from_seurat(dest)
    }
  }

  if (!is.null(pool)) {
    pool$source_info <- shiny_demo_source_info()
    .shiny_demo_pool_cache$pool <- pool
    return(pool)
  }

  NULL
}

#' Clear the in-memory CRA001160 demo pool cache (mainly for tests).
#' @keywords internal
clear_shiny_demo_pool_cache <- function() {
  rm(list = ls(envir = .shiny_demo_pool_cache), envir = .shiny_demo_pool_cache)
}

#' @keywords internal
.make_shiny_demo_dataset_synthetic <- function(n_genes = 500L, n_cells = 5000L, seed = 7L) {
  set.seed(seed)
  cell_types <- c("Acinar", "Ductal", "Macrophage", "CD8T", "Fibroblast")
  sources <- c("Tumor", "Normal")
  colnames_cells <- paste0("Cell_", seq_len(n_cells))
  cell_type <- sample(cell_types, n_cells, replace = TRUE,
                      prob = c(0.28, 0.22, 0.12, 0.18, 0.20))
  source <- ifelse(cell_type %in% c("Acinar", "Ductal"),
                   sample(sources, n_cells, replace = TRUE, prob = c(0.7, 0.3)),
                   sample(sources, n_cells, replace = TRUE, prob = c(0.55, 0.45)))

  hallmark <- c(
    "TP53", "MYC", "EGFR", "BRCA1", "CDKN2A", "KRAS", "PIK3CA", "PTEN",
    "MKI67", "CD68", "CD3D", "CD8A", "COL1A1", "KRT19", "AMY2A", "S100A8",
    "FOXM1", "KLRB1", "CXCL8", "IL6", "VEGFA", "STAT1", "HLA-DRA", "GZMB"
  )
  precog_z <- .precog_pancreatic_gene_z()
  precog_genes <- if (!is.null(precog_z)) {
    names(precog_z)[abs(precog_z) >= 1]
  } else {
    character(0)
  }
  genes <- unique(c(hallmark, precog_genes))
  if (length(genes) < n_genes) {
    genes <- c(genes, paste0("GENE", seq_len(n_genes - length(genes))))
  }
  genes <- genes[seq_len(min(n_genes, length(genes)))]

  base <- matrix(stats::rnorm(length(genes) * n_cells, mean = 0.45, sd = 0.3),
                 nrow = length(genes), ncol = n_cells)
  rownames(base) <- genes
  colnames(base) <- colnames_cells

  type_idx <- split(seq_len(n_cells), cell_type)
  for (ct in names(type_idx)) {
    bump_genes <- switch(
      ct,
      Acinar = c("AMY2A", "MYC"),
      Ductal = c("KRT19", "FOXM1"),
      Macrophage = c("CD68", "S100A8", "CXCL8"),
      CD8T = c("CD3D", "CD8A", "GZMB", "KLRB1"),
      Fibroblast = c("COL1A1", "IL6"),
      character(0)
    )
    bump_genes <- intersect(bump_genes, genes)
    for (g in bump_genes) {
      base[g, type_idx[[ct]]] <- base[g, type_idx[[ct]]] +
        stats::rnorm(length(type_idx[[ct]]), mean = 1.2, sd = 0.35)
    }
  }

  expr <- pmax(base, 0)
  zero_mask <- matrix(
    stats::rbinom(length(expr), size = 1L, prob = 0.58),
    nrow = nrow(expr), ncol = ncol(expr)
  )
  expr[zero_mask == 1L] <- 0

  gene_var <- apply(expr, 1L, stats::var)
  top_h <- order(gene_var, decreasing = TRUE)[seq_len(min(50L, nrow(expr)))]
  cell_mat <- t(expr[top_h, , drop = FALSE])
  pcs <- stats::prcomp(cell_mat, center = TRUE, scale. = TRUE)
  umap1 <- pcs$x[, 1L] + stats::rnorm(n_cells, sd = 0.15)
  umap2 <- if (ncol(pcs$x) >= 2L) {
    pcs$x[, 2L] + stats::rnorm(n_cells, sd = 0.15)
  } else {
    stats::rnorm(n_cells, sd = 0.15)
  }

  metadata <- data.frame(
    .cell_id = colnames_cells,
    cell_type = cell_type,
    cell_type_minor = cell_type,
    Source = source,
    UMAP_1 = umap1,
    UMAP_2 = umap2,
    stringsAsFactors = FALSE
  )
  expr <- .attach_demo_matrix_obsm(expr, metadata)
  info <- shiny_demo_source_info()
  info$title <- paste(info$title, "(synthetic fallback)")
  info$n_cells_sampled <- n_cells
  info$n_genes_sampled <- nrow(expr)
  info$n_cells_pool <- n_cells
  list(expression = expr, metadata = metadata, source_info = info, from = "synthetic")
}

#' Load the pre-built Shiny demo bundle (5,000 pre-selected CRA001160 cells).
#'
#' Falls back to building from the full CRA001160 pool or a synthetic matrix
#' when the bundled RDS is unavailable (mainly for development).
#'
#' @return List with \code{expression}, \code{metadata}, and \code{source_info}.
#' @keywords internal
make_shiny_demo_dataset <- function(n_genes = 500L, n_cells = 5000L, seed = NULL) {
  bundle_path <- .shiny_demo_bundle_path()
  if (!is.null(bundle_path)) {
    demo <- tryCatch(readRDS(bundle_path), error = function(e) NULL)
    if (!is.null(demo) && !is.null(demo$expression) && !is.null(demo$metadata)) {
      expr <- demo$expression
      md <- demo$metadata
      if (!"cell_type" %in% names(md) && "cell_type_minor" %in% names(md)) {
        md$cell_type <- md$cell_type_minor
      }
      expr <- .attach_demo_matrix_obsm(expr, md)
      info <- demo$source_info %||% shiny_demo_source_info()
      info$is_presampled <- TRUE
      info$n_cells_sampled <- ncol(expr)
      info$n_genes_sampled <- nrow(expr)
      return(list(
        expression = expr,
        metadata = md,
        source_info = info,
        from = demo$from %||% "CRA001160"
      ))
    }
  }

  pool <- load_shiny_demo_pool()
  if (is.null(pool)) {
    if (is.null(seed)) seed <- 20250625L
    return(.make_shiny_demo_dataset_synthetic(n_genes = n_genes, n_cells = n_cells, seed = seed))
  }

  if (is.null(seed)) seed <- 20250625L
  set.seed(seed)
  cells <- colnames(pool$expression)
  n_draw <- min(as.integer(n_cells), length(cells))
  keep <- sample(cells, n_draw, replace = FALSE)

  expr_sub <- pool$expression[, keep, drop = FALSE]
  n_keep_genes <- min(as.integer(n_genes), nrow(expr_sub))
  gene_var <- .shiny_row_var(expr_sub)
  top_genes <- order(gene_var, decreasing = TRUE)[seq_len(n_keep_genes)]
  expr <- if (inherits(expr_sub, "Matrix")) {
    as.matrix(expr_sub[top_genes, , drop = FALSE])
  } else {
    expr_sub[top_genes, , drop = FALSE]
  }

  md <- pool$metadata[match(keep, pool$metadata$.cell_id), , drop = FALSE]
  rownames(md) <- NULL
  minor <- md$cell_type_minor
  if (!is.null(minor)) md$cell_type <- as.character(minor)
  expr <- .attach_demo_matrix_obsm(expr, md)

  info <- pool$source_info
  info$n_cells_pool <- length(cells)
  info$n_cells_sampled <- n_draw
  info$n_genes_sampled <- nrow(expr)
  info$bundle_seed <- seed

  list(
    expression = expr,
    metadata = md,
    source_info = info,
    from = pool$from
  )
}

build_cell_table <- function(scores,
                              groups = NULL,
                              metadata = NULL,
                              cell_id_col = NULL,
                              cell_type_col = NULL,
                              source_col = NULL,
                              score_column = NULL) {
  if (is.null(scores)) return(NULL)
  s <- as.data.frame(scores)
  s$cell_id <- rownames(scores)
  score_cols <- setdiff(colnames(s), "cell_id")
  numeric_cols <- score_cols[vapply(s[score_cols], is.numeric, logical(1))]
  if (!length(numeric_cols)) return(NULL)
  score_name <- if (!is.null(score_column) && nzchar(score_column) &&
                    score_column %in% numeric_cols) {
    score_column
  } else {
    numeric_cols[1L]
  }
  out <- data.frame(
    cell_id = s$cell_id,
    score   = as.numeric(s[[score_name]]),
    stringsAsFactors = FALSE
  )
  attr(out, "score_name") <- score_name

  if (!is.null(groups) && "cell_id" %in% colnames(groups)) {
    grp_col <- resolve_phenotype_group_column(groups, score_name)
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

    # Always store cell_type / source as character labels. Integers (e.g. FOV
    # / core / cluster IDs) would otherwise stay numeric and ggplot2 would
    # treat them as continuous when applying discrete PhenoMapR scales.
    if (!is.null(cell_type_col) && nzchar(cell_type_col) && cell_type_col != "(none)" &&
        cell_type_col %in% colnames(md)) {
      out$.cell_type <- as.character(md[match(out$cell_id, md$cell_id), cell_type_col])
      names(out)[names(out) == ".cell_type"] <- "cell_type"
    }
    if (!is.null(source_col) && nzchar(source_col) && source_col != "(none)" &&
        source_col %in% colnames(md)) {
      out$.source <- as.character(md[match(out$cell_id, md$cell_id), source_col])
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
  clear_id  <- paste0(picker_id, "_clear")

  # pick_state is the *effective* pick exposed to the rest of the app:
  # NULL means "no file selected", otherwise a list with datapath/name.
  # Two write paths into it:
  #   1. The shinyFiles picker emits a new selection -> pick_state set.
  #   2. The user clicks the small "x" remove button next to the
  #      displayed filename -> pick_state set to NULL.
  # Downstream observers depend on pick_state via the returned reactive,
  # so clearing here propagates everywhere just like a fresh pick would.
  pick_state <- shiny::reactiveVal(NULL)

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

    # shinyFiles selection -> pick_state.
    #
    # IMPORTANT (sticky picks): we do NOT clear pick_state on a
    # NULL / integer(0) / parse-failure update from the picker.
    # shinyFiles emits an empty input value any time the user opens
    # the file dialog and cancels (or sometimes when the picker is
    # re-armed mid-session); previously that would silently unload
    # the user\'s current file the moment they clicked "Browse..."
    # again -- even if they had no intention of switching files.
    # The explicit "x" remove button below is the ONLY user gesture
    # that drops a loaded file, plus selecting a new file (which
    # naturally replaces the old pick via the success branch here).
    shiny::observeEvent(input[[picker_id]], {
      sel <- input[[picker_id]]
      if (is.null(sel) || identical(sel, integer(0))) {
        # Browse dialog was opened then cancelled; keep the existing
        # pick. (No-op.)
        return()
      }
      path_df <- tryCatch(
        shinyFiles::parseFilePaths(roots, sel),
        error = function(e) NULL
      )
      if (is.null(path_df) || nrow(path_df) == 0L) {
        # Couldn\'t materialise a real file path from the selection;
        # also keep the existing pick rather than silently dropping
        # the loaded file.
        return()
      }
      pick_state(list(
        datapath = as.character(path_df$datapath[1]),
        name = as.character(path_df$name[1] %||% basename(path_df$datapath[1])),
        source = "server"
      ))
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    # Remove-file button -> clear pick_state.
    #
    # IMPORTANT: the clear button is rendered dynamically by the
    # `_chosen` renderUI below -- it only exists in the DOM after a
    # file has been picked. When Shiny first binds the button it sets
    # input[[clear_id]] from NULL to 0 (the initial actionButton
    # value). ignoreInit only suppresses the VERY FIRST evaluation
    # (when input was NULL before the button existed), so the NULL-to-0
    # transition still fires the observer. Without this guard, picking
    # a file would auto-clear it a moment later via that spurious 0
    # event, and the metadata never stays in the session.
    shiny::observeEvent(input[[clear_id]], {
      n <- input[[clear_id]] %||% 0L
      if (!is.numeric(n) || n < 1L) return()
      pick_state(NULL)
    }, ignoreInit = TRUE)

    output[[paste0(picker_id, "_chosen")]] <- shiny::renderUI({
      cur <- pick_state()
      if (is.null(cur)) {
        return(shiny::tags$span(
          class = "phenomapr-file-input-name-empty",
          "No file selected"
        ))
      }
      shiny::tags$span(
        class = "phenomapr-file-input-name-set",
        title = cur$datapath,
        shiny::tags$span(
          class = "phenomapr-file-input-name-text",
          basename(cur$datapath)
        ),
        # The "remove" affordance: a tiny x button that resets the
        # picker back to "No file selected" and propagates NULL through
        # the pick reactive (which downstream observers depend on).
        # We use a hand-rolled <button> instead of actionButton() so it
        # inherits the surrounding flex layout cleanly. The trailing
        # ".btn" class lets Bootstrap style it consistently if present.
        shiny::tags$button(
          id = clear_id,
          type = "button",
          class = "btn action-button phenomapr-file-input-clear",
          title = "Remove loaded file",
          `aria-label` = "Remove loaded file",
          shiny::HTML("&times;")
        )
      )
    })

    return(shiny::reactive({
      pick_state()
    }))
  }

  # Fallback: shinyFiles not installed -> a plain fileInput was rendered
  # under `id`. Still surface a remove affordance and clear semantics.
  shiny::observeEvent(input[[id]], {
    up <- input[[id]]
    if (is.null(up) || nrow(up) == 0L) {
      pick_state(NULL)
      return()
    }
    pick_state(list(
      datapath = as.character(up$datapath[1]),
      name = as.character(up$name[1]),
      source = "upload"
    ))
  }, ignoreInit = TRUE, ignoreNULL = FALSE)

  shiny::observeEvent(input[[clear_id]], {
    pick_state(NULL)
  }, ignoreInit = TRUE)

  shiny::reactive({ pick_state() })
}

# ============================================================================
# Panel download affordance (plots + tables)
# ============================================================================
#
# Every plot and table panel in the app surfaces a small download button
# in the right-hand side of its card header (or a slim banner row when
# the panel sits inside a nav_panel that has no card chrome). The
# helpers are:
#
#   phenomapr_dl_btn(download_id, tooltip)
#     UI helper. A small icon-only downloadButton with the
#     `.phenomapr-panel-download-btn` class. Used directly inside
#     card headers and banner rows.
#
#   phenomapr_card_header_dl(..., download_id, tooltip)
#     UI helper. Returns a bslib::card_header() with the title
#     content on the left and the download button pinned to the
#     right via flex layout (see .phenomapr-card-header-dl in
#     styles.css).
#
#   phenomapr_panel_banner_dl(download_id, tooltip)
#     UI helper. Returns a slim, right-aligned banner row for cases
#     where the panel has no card header (e.g., nav_panel children).
#
#   phenomapr_register_plot_download(output, id, plot_fn,
#                                    width, height, dpi)
#     Server helper. Registers `output[[paste0(id, "_download")]]` as a
#     downloadHandler that calls `plot_fn()` for the current ggplot
#     object and saves it via ggplot2::ggsave(). plot_fn is typically a
#     reactiveVal getter that the renderPlot body has been instrumented
#     to populate (see "capture pattern" below).
#
#   phenomapr_register_table_download(output, id, data_fn)
#     Server helper. Registers `output[[paste0(id, "_download")]]` as a
#     downloadHandler that calls `data_fn()` for the current data.frame
#     and writes it as TSV.
#
# Capture pattern in renderPlot / renderDT bodies:
#   render*({ ...; p <- ggplot(...) + ...; rv(p); p })
# i.e., the body builds the plot/table, *captures* it via the
# reactiveVal `rv`, and returns it as usual. Downstream
# downloadHandler reads rv() at click time. Capture-on-render works
# because the user can only click the download button after the panel
# has rendered.

phenomapr_dl_btn <- function(download_id, tooltip = "Download",
                             label = NULL) {
  # When `label` is provided we render an icon + visible text button
  # (still styled by .phenomapr-panel-download-btn but with the
  # `--labeled` modifier so it widens to fit the label and the icon
  # gets a small right gap). Callers that want the compact 22x22
  # icon-only button just omit `label` -- this is what every other
  # panel header has been doing forever, so no migration is needed.
  classes <- if (!is.null(label) && nzchar(label)) {
    "phenomapr-panel-download-btn phenomapr-panel-download-btn-labeled"
  } else {
    "phenomapr-panel-download-btn"
  }
  shiny::downloadButton(
    download_id,
    label = label,
    icon  = shiny::icon("download"),
    class = classes,
    title = tooltip
  )
}

phenomapr_card_header_dl <- function(..., download_id, tooltip = "Download") {
  bslib::card_header(
    class = "phenomapr-card-header-dl",
    shiny::tags$div(class = "phenomapr-card-header-dl-title", ...),
    phenomapr_dl_btn(download_id, tooltip)
  )
}

phenomapr_panel_banner_dl <- function(download_id, tooltip = "Download") {
  shiny::tags$div(
    class = "phenomapr-panel-banner",
    phenomapr_dl_btn(download_id, tooltip)
  )
}

# ---- Plot-only modal-trigger download buttons -----------------------------
#
# Plots get a download icon that *opens a modal* instead of starting
# the download immediately. The modal exposes width / height / DPI /
# file format / base font size controls and contains the actual
# `downloadButton("plot_dl_action", ...)` whose downloadHandler is
# registered once at server start. Tables keep the direct
# downloadButton flow above -- one click -> file -- since they have
# fewer customisation knobs.
#
# Wiring:
#   - phenomapr_modal_dl_btn(panel_id) renders an actionButton with
#     inputId = paste0(panel_id, "_modal_btn") and shares the
#     `.phenomapr-panel-download-btn` styling so it visually matches
#     the table download buttons.
#   - app.R server registers one observeEvent per plot panel that
#     captures the click, stashes `panel_id` in a reactiveVal, and
#     calls showModal(phenomapr_plot_download_modal(...)).
#   - The shared `output$plot_dl_action` downloadHandler reads the
#     stashed panel id and the modal inputs (width/height/dpi/format/
#     base_size) to produce the file via .phenomapr_save_plot().

phenomapr_modal_dl_btn <- function(panel_id,
                                   tooltip = "Download (with options)",
                                   label = NULL) {
  # See phenomapr_dl_btn for the labeled/icon-only convention; the
  # modal trigger uses the same `--labeled` modifier so the two
  # buttons in a `data + plot` header line up visually.
  classes <- if (!is.null(label) && nzchar(label)) {
    "phenomapr-panel-download-btn phenomapr-panel-download-btn-labeled"
  } else {
    "phenomapr-panel-download-btn"
  }
  shiny::actionButton(
    inputId = paste0(panel_id, "_modal_btn"),
    label   = label,
    icon    = shiny::icon("download"),
    class   = classes,
    title   = tooltip
  )
}

phenomapr_card_header_modal_dl <- function(..., panel_id,
                                           tooltip = "Download (with options)",
                                           plot_label = NULL,
                                           data_download_id = NULL,
                                           data_tooltip = "Download plot data (TSV)",
                                           data_label = NULL) {
  # If a panel exposes both its rendered plot AND its underlying
  # tabular data for download (e.g. "Score by cell type and source",
  # whose data is the per-cell score table), we render TWO download
  # affordances clustered at the right end of the header: a plain
  # downloadButton for the data (`data_download_id`) followed by the
  # modal-trigger button for the plot (`panel_id`). They live inside
  # a `.phenomapr-card-header-dl-actions` flex wrapper so the
  # existing `margin-left: auto` styling that pins a single button to
  # the right still works without leaving a wide gap between the two
  # buttons (which the default `justify-content: space-between` would
  # otherwise insert).
  #
  # `plot_label` / `data_label` (both default NULL) opt the buttons
  # into the labeled style: when set, the button shows icon + text
  # so users with two adjacent download buttons can tell them apart
  # at a glance ("Download plot" vs "Download plot data") instead of
  # having to hover for the tooltip.
  modal_btn <- phenomapr_modal_dl_btn(panel_id, tooltip, label = plot_label)
  title_div <- shiny::tags$div(class = "phenomapr-card-header-dl-title", ...)
  if (is.null(data_download_id)) {
    bslib::card_header(
      class = "phenomapr-card-header-dl",
      title_div,
      modal_btn
    )
  } else {
    bslib::card_header(
      class = "phenomapr-card-header-dl",
      title_div,
      shiny::tags$div(
        class = "phenomapr-card-header-dl-actions",
        phenomapr_dl_btn(data_download_id, data_tooltip,
                         label = data_label),
        modal_btn
      )
    )
  }
}

phenomapr_panel_banner_modal_dl <- function(panel_id,
                                            tooltip = "Download (with options)") {
  shiny::tags$div(
    class = "phenomapr-panel-banner",
    phenomapr_modal_dl_btn(panel_id, tooltip)
  )
}

# Builds the modal dialog. `panel_label` is shown in the title; the
# `defaults` list seeds the width / height / dpi / format / base_size
# inputs so each panel can suggest sensible export dimensions.
phenomapr_plot_download_modal <- function(panel_label = "plot",
                                          defaults = list()) {
  d <- utils::modifyList(
    list(width = 8, height = 6, dpi = 300,
         format = "png", base_size = 14, radius_pt = 3,
         units = "in"),
    defaults
  )
  # Build a flex header that puts the dialog title on the left and a
  # right-aligned Cancel / Download actions cluster. We pass this in
  # via `title` (because Bootstrap's modal-header already lays its
  # children out as a flex row) and clear out `footer` so the buttons
  # only live in the header. The implicit `.close` (x) glyph in the
  # top-right corner stays as a third dismiss affordance.
  title_row <- shiny::tags$div(
    class = "phenomapr-plot-dl-modal-title-row",
    shiny::tags$span(
      class = "phenomapr-plot-dl-modal-title-text",
      shiny::icon("download"),
      " ",
      paste0("Download plot: ", panel_label)
    ),
    shiny::tags$div(
      class = "phenomapr-plot-dl-modal-actions",
      shiny::modalButton("Cancel"),
      shiny::downloadButton("plot_dl_action", "Download",
                            class = "btn-primary")
    )
  )
  ## Convenience labels for Width / Height that carry the active unit.
  unit_lbl <- function(side) sprintf("%s (%s)", side,
                                     if (identical(d$units, "cm")) "cm" else "inches")
  shiny::modalDialog(
    title = title_row,
    size  = "l",
    easyClose = TRUE,
    shiny::div(
      class = "phenomapr-plot-dl-modal",
      # ---- Live preview ---------------------------------------------------
      # Mirrors what the saved file will look like at the chosen width /
      # height / base font size. DPI and file format do not affect the
      # on-screen preview (vector / raster + DPI are file-only) -- that
      # caveat is called out in the helpText below. The renderPlot in
      # app.R reads `active_plot_dl()` + the modal inputs and re-runs
      # whenever they change, so users can dial in their layout before
      # ever clicking Download.
      shiny::div(
        class = "phenomapr-plot-dl-preview",
        shiny::tags$div(
          class = "phenomapr-plot-dl-preview-title",
          shiny::icon("eye"), " Preview"
        ),
        shiny::div(
          class = "phenomapr-plot-dl-preview-frame",
          shiny::plotOutput("plot_dl_preview",
                            width = "100%", height = "320px")
        ),
        shiny::tags$small(
          class = "text-muted",
          "Updates live as you change width / height / base font size. ",
          "DPI and file format affect the saved file only -- not the preview."
        )
      ),
      shiny::tags$hr(class = "phenomapr-plot-dl-divider"),
      # ---- Sizing & format controls ---------------------------------------
      # Width / height live in a row WITH a unit selector (inches vs cm).
      # The selector is wired in app.R: when the user toggles units, the
      # width / height fields are rescaled (x 2.54 to cm, /2.54 to in) so
      # the physical size stays constant and the labels update to match.
      shiny::fluidRow(
        shiny::column(4,
          shiny::radioButtons(
            "plot_dl_units", "Units",
            choices  = c("inches" = "in", "cm" = "cm"),
            selected = d$units,
            inline   = TRUE
          )
        ),
        shiny::column(4,
          shiny::numericInput("plot_dl_width", unit_lbl("Width"),
                              value = d$width, min = 1,
                              max = if (identical(d$units, "cm")) 80 else 30,
                              step = if (identical(d$units, "cm")) 1 else 0.5)),
        shiny::column(4,
          shiny::numericInput("plot_dl_height", unit_lbl("Height"),
                              value = d$height, min = 1,
                              max = if (identical(d$units, "cm")) 80 else 30,
                              step = if (identical(d$units, "cm")) 1 else 0.5))
      ),
      shiny::fluidRow(
        shiny::column(6,
          shiny::selectInput(
            "plot_dl_format", "File format",
            choices  = c("PNG (raster)"  = "png",
                         "JPEG (raster)" = "jpeg",
                         "TIFF (raster)" = "tiff",
                         "PDF (vector)"  = "pdf",
                         "SVG (vector)"  = "svg"),
            selected = d$format
          )),
        shiny::column(6,
          shiny::numericInput("plot_dl_dpi", "DPI (raster only)",
                              value = d$dpi, min = 72, max = 600, step = 10))
      ),
      shiny::fluidRow(
        shiny::column(6,
          shiny::sliderInput("plot_dl_base_size",
                             "Plot base font size",
                             min = 8, max = 24,
                             value = d$base_size, step = 1)),
        shiny::column(6,
          shiny::sliderInput("plot_dl_radius_pt",
                             "Bar / box corner sharpness (pt)",
                             min = 0, max = 5,
                             value = min(d$radius_pt, 5),
                             step = 0.5))
      ),
      shiny::helpText(
        "Choose your preferred unit (inches or cm) -- width / height are ",
        "rescaled automatically when you switch so the physical size stays ",
        "constant. DPI only affects raster formats (PNG / JPEG / TIFF); ",
        "PDF and SVG are vector formats and render at any zoom. Base font ",
        "size and corner sharpness affect this export only (the on-screen ",
        "plot is unchanged). Corner sharpness applies to ggchicklet2-rounded ",
        "bars, stacks, histograms and boxplots; non-chicklet plots ",
        "ignore it gracefully."
      )
    ),
    # Footer-less: Cancel + Download live in the title row above.
    footer = NULL
  )
}

# Export-dispatch: write a captured plot object to file in the user-
# chosen format. ggplot / patchwork objects go through ggplot2::ggsave
# (which handles png / jpeg / tiff / pdf / svg via the `device` arg);
# anything else (ComplexHeatmap, recordedplot, plain grobs) falls back
# to grDevices::<device>() + print(). When `base_size` is set we apply
# a theme override at export time only so the on-screen plot is
# unchanged.
.phenomapr_save_plot <- function(file, plot_obj, format,
                                 width, height, dpi,
                                 base_size = NULL,
                                 radius_pt = NULL) {
  format <- tolower(format)
  is_gg  <- inherits(plot_obj, c("ggplot", "patchwork"))

  if (is_gg && !is.null(base_size) && is.numeric(base_size) &&
      is.finite(base_size) && base_size > 0) {
    plot_obj <- plot_obj +
      ggplot2::theme(text = ggplot2::element_text(size = base_size))
  }
  if (is_gg && !is.null(radius_pt)) {
    plot_obj <- .apply_chicklet_radius(plot_obj, radius_pt)
  }

  if (is_gg) {
    ggplot2::ggsave(
      file, plot = plot_obj, device = format,
      width = width, height = height,
      units = "in", dpi = dpi,
      limitsize = FALSE
    )
    return(invisible(NULL))
  }

  switch(format,
    png  = grDevices::png(file, width = width * dpi, height = height * dpi, res = dpi),
    jpeg = grDevices::jpeg(file, width = width * dpi, height = height * dpi, res = dpi, quality = 95),
    tiff = grDevices::tiff(file, width = width * dpi, height = height * dpi, res = dpi),
    pdf  = grDevices::pdf(file, width = width, height = height),
    svg  = grDevices::svg(file, width = width, height = height),
    stop("Unsupported format: ", format)
  )
  on.exit(if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(),
                                              silent = TRUE),
          add = TRUE)
  print(plot_obj)
  invisible(NULL)
}

phenomapr_dl_filename <- function(stem, ext) {
  sprintf("phenomapr_%s_%s.%s",
          stem,
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ext)
}

# Write a small placeholder PNG when there is no plot to download.
.phenomapr_write_placeholder_png <- function(file, msg = "No plot to download") {
  grDevices::png(file, width = 480, height = 240, res = 110)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::text(0.5, 0.5, msg, cex = 1.2)
}

phenomapr_register_plot_download <- function(output, id, plot_fn,
                                             width = 8, height = 6, dpi = 150) {
  output[[paste0(id, "_download")]] <- shiny::downloadHandler(
    filename = function() phenomapr_dl_filename(id, "png"),
    content = function(file) {
      p <- tryCatch(plot_fn(), error = function(e) NULL)
      if (is.null(p)) {
        .phenomapr_write_placeholder_png(file)
        return(invisible(NULL))
      }
      # ggplot2::ggsave handles ggplot objects natively. For other
      # plot-like objects (e.g. lattice, base R), fall back to a
      # direct PNG render via grDevices::png().
      if (inherits(p, c("ggplot", "ggplot2::ggplot", "patchwork"))) {
        ggplot2::ggsave(file, plot = p,
                        width = width, height = height, dpi = dpi,
                        units = "in", limitsize = FALSE)
      } else {
        grDevices::png(file,
                       width = width * dpi, height = height * dpi,
                       res = dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(p)
      }
      invisible(NULL)
    }
  )
}

phenomapr_register_table_download <- function(output, id, data_fn,
                                              filename = NULL) {
  output[[paste0(id, "_download")]] <- shiny::downloadHandler(
    filename = function() {
      if (!is.null(filename)) return(filename())
      phenomapr_dl_filename(id, "tsv")
    },
    content = function(file) {
      df <- tryCatch(data_fn(), error = function(e) NULL)
      if (is.null(df) || (is.data.frame(df) && nrow(df) == 0L)) {
        writeLines("# No rows to download", file)
        return(invisible(NULL))
      }
      utils::write.table(df, file, sep = "\t",
                         row.names = FALSE, quote = FALSE,
                         na = "")
      invisible(NULL)
    }
  )
}

# ============================================================================
# Centered "busy" overlay (R-driven custom message + client debounce)
# ============================================================================
#
# Timing model
# ------------
# The popup is GATED by an explicit phenomapr_busy_show() call from R
# plus a SHOW_DELAY_MS debounce on the JS side; it is hidden by an
# explicit phenomapr_busy_hide() call (or any of several safety nets
# described below).
#
# Why not gate purely on Shiny's busy/idle state? Because R is single
# threaded and a long synchronous observer body (e.g. reading a
# multi-GB .h5 file, computing PhenoMap scores) BLOCKS libuv from
# delivering the websocket "busy: true" / "busy: false" frames to the
# browser until *after* the work completes. The browser would then
# receive the busy and idle messages back-to-back and the user would
# never see a popup, no matter how long the operation took. That is
# exactly the bug the previous busy-class-poll model could not
# escape: there's no point polling .shiny-busy if the server can't
# actually update it during the slow work.
#
# Instead, the public R API is:
#
#   phenomapr_busy_show("Computing PhenoMap scores...", "Cancer: LUAD")
#   ...heavy work scheduled into later::later(...)...
#   phenomapr_busy_hide()
#
# The CRITICAL discipline call-sites must follow is to schedule the
# heavy synchronous work into a `later::later(fn, delay = 0)`
# callback. This lets flushReact for the OBSERVER complete (which
# is what actually sends the busy_show custom message to the
# browser via the websocket) BEFORE libuv is blocked by the heavy
# work. Without that deferral, busy_show and busy_hide queue
# together inside the same blocking observer and arrive at the
# client back-to-back -- popup never appears. See
# observeEvent(input$run_score, ...) etc. in app.R for the canonical
# pattern.
#
# JS-side behaviour:
#   * handleShow(payload):  store the text, schedule a SHOW_DELAY_MS
#                            timer. If still no hide arrived when the
#                            timer fires, render the popup. (Short
#                            ops where show + hide arrive within
#                            SHOW_DELAY_MS therefore never flash.)
#   * handleHide():          cancel any pending show timer; hide
#                            popup if visible; clear stored text.
#   * shiny:idle:            advisory hide (defence-in-depth) --
#                            even if R forgot to call busy_hide(),
#                            the popup goes away when shiny goes
#                            idle for at least IDLE_SETTLE_MS.
#   * shiny:disconnected:    force hide.
#   * Absolute 3-min cap:    last-resort hide if a popup somehow
#                            survives all of the above.
#
# The `delay_seconds` parameter on phenomapr_busy_show() is retained
# for back-compat but IGNORED. The client-side SHOW_DELAY_MS debounce
# is the single source of truth for show-side timing.

phenomapr_busy_show <- function(message,
                                detail = NULL,
                                hint = "This may take a moment...",
                                delay_seconds = 0,
                                session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) return(invisible(NULL))
  session$sendCustomMessage("phenomapr-busy-show", list(
    message = as.character(message %||% ""),
    detail  = if (is.null(detail) || !nzchar(detail)) "" else as.character(detail),
    hint    = if (is.null(hint)   || !nzchar(hint))   "" else as.character(hint)
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
#   - Visibility is driven by R-side custom messages (phenomapr-busy-show
#     and phenomapr-busy-hide) with a SHOW_DELAY_MS debounce on the
#     show side. Fast ops where show + hide arrive within
#     SHOW_DELAY_MS never flash the popup. Long ops show the popup
#     after the debounce and hide it the moment R signals end.
#   - For the popup to actually appear during a slow synchronous
#     observer, R MUST schedule the heavy work via
#     later::later(fn, delay = 0) AFTER calling phenomapr_busy_show().
#     This lets flushReact deliver the show message to the browser
#     before libuv gets blocked by the heavy work; without it,
#     show + hide queue together server-side and arrive back-to-back.
#     See the canonical pattern in inst/shiny/app.R observers.
#   - shiny:idle is wired as an ADVISORY safety hide (defence in
#     depth, in case R forgot busy_hide()). It does NOT decide show
#     timing; only the explicit phenomapr-busy-show message does.
#   - A 3-minute absolute hard cap force-hides any popup that
#     somehow survives all of the above.
#   - Append `?busyDebug=1` to the URL to enable console.debug logging
#     of every visibility transition. Useful when diagnosing future
#     timing reports without touching the code.
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
  // Visibility model
  // ----------------
  // 1. R calls phenomapr_busy_show("...", "...") which sends a
  //    `phenomapr-busy-show` custom message. handleShow() stores the
  //    text and schedules a SHOW_DELAY_MS timer.
  // 2. If R schedules its heavy work into later::later(fn, delay=0),
  //    flushReact for the calling observer completes immediately,
  //    libuv flushes the websocket queue, and the show message
  //    reaches the browser. The browser now has SHOW_DELAY_MS to
  //    decide whether the op is slow enough to visualize.
  // 3. When the timer fires, if R has not already sent a
  //    `phenomapr-busy-hide` (which would have cancelled the timer),
  //    the popup is rendered. R\'s heavy synchronous work may still
  //    be in progress at this point; that is exactly when we WANT
  //    the popup visible.
  // 4. When the heavy work completes, R calls phenomapr_busy_hide().
  //    handleHide() cancels any pending timer and hides the popup
  //    if visible.
  //
  // Defence-in-depth: shiny:idle also triggers a (delayed) hide,
  // shiny:disconnected forces an immediate hide, and an absolute
  // 3-minute cap force-hides the popup if everything else fails.
  var SHOW_DELAY_MS     = 600;    // debounce before showing after busy_show
  var IDLE_SETTLE_MS    = 200;    // delay before advisory shiny:idle hide
  var ABSOLUTE_TIMEOUT  = 180000; // last-resort hard cap (3 min)
  var FILE_INPUT_SUFFIX = "_server"; // shinyFiles inputs created by phenomapr_file_input

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
      // Clear any inline `display: none` that forceHide() may have
      // pinned, so the .is-visible class can take effect again.
      ov.style.display = "";
      ov.classList.add("is-visible");
      isVisible = true;
      shownAt = nowMs();
      dbg("overlay shown", currentMessage);
    }
  }
  function renderHide() {
    // Forwarder kept for code paths that still call renderHide()
    // directly. The actual logic lives in forceHide() so every hide
    // path in this file goes through the same DOM probe + class-strip.
    forceHide("renderHide");
  }
  function nowMs() {
    return (typeof performance !== "undefined" && performance.now)
      ? performance.now()
      : Date.now();
  }
  // ---- State ---------------------------------------------------------
  var defaultMessage    = { message: "Working...", detail: "", hint: "" };
  var currentMessage    = { message: "Working...", detail: "", hint: "" };
  var isVisible         = false;
  var shownAt           = 0;     // ms timestamp when popup last became visible
  var dismissedThisRun  = false; // user-clicked dismiss; cleared by handleHide / shiny:idle
  var pendingShowTimer  = null;  // setTimeout id for the SHOW_DELAY_MS debounce
  var pendingShow       = false; // R has called busy_show but the timer has not yet fired
  var idleHideTimer     = null;  // setTimeout id for the advisory shiny:idle hide
  var lastShowAt        = 0;     // ms timestamp of the most recent R-side busy_show
  var lastHideAt        = 0;     // ms timestamp of the most recent R-side busy_hide

  function resetMessage() {
    currentMessage = {
      message: defaultMessage.message,
      detail:  defaultMessage.detail,
      hint:    defaultMessage.hint
    };
  }

  function clearPendingShow() {
    if (pendingShowTimer !== null) {
      clearTimeout(pendingShowTimer);
      pendingShowTimer = null;
    }
    pendingShow = false;
  }
  function clearIdleHide() {
    if (idleHideTimer !== null) {
      clearTimeout(idleHideTimer);
      idleHideTimer = null;
    }
  }

  // ---- Shiny advisory listeners --------------------------------------
  // shiny:idle is an ADVISORY safety hide -- if R signals busy_show
  // and then forgets to call busy_hide (e.g. an observer body
  // throws before reaching its on.exit), the popup still goes away
  // when shiny becomes idle. We delay by IDLE_SETTLE_MS to absorb
  // the tiny gap between an observer returning and a deferred
  // later::later() callback firing the next reactive cascade.
  // shiny:busy is intentionally NOT used to drive show timing --
  // see the file comment above for why.
  function onShinyBusy() {
    dbg("shiny:busy event");
    // We deliberately do NOT clear the in-flight idle hide here.
    // pendingShow gates the actual hide -- if the user is still
    // inside the SHOW_DELAY_MS debounce window, the timer will be
    // a no-op. After SHOW_DELAY_MS the popup is rendered and
    // pendingShow is false, so any subsequent shiny:idle (e.g. the
    // tail of a heavy cascade after busy_hide arrived) hides the
    // popup. Cancelling on every shiny:busy used to defeat this:
    // a cascade with multiple back-to-back flushReacts could keep
    // re-cancelling the hide and leave the popup stuck on screen.
  }
  function onShinyIdle() {
    dbg("shiny:idle event");
    clearAllFileInputLoading();
    // dismissedThisRun is per-busy-cycle; clear once the cycle ends.
    dismissedThisRun = false;

    // FAST PATH: popup is visible AND pendingShow is false. This
    // means the SHOW_DELAY_MS debounce has already elapsed (popup
    // is fully rendered) and R is now idle -- which can only happen
    // AFTER our heavy-work later() callback completed and any
    // resulting reactive cascade finished. Hide immediately, no
    // need to wait IDLE_SETTLE_MS more.
    if (overlayDomVisible() && !pendingShow) {
      console.log("[phenomapr-busy] shiny:idle -> immediate hide (cascade complete)");
      clearIdleHide();
      clearPendingShow();
      forceHide("shiny:idle (popup up, cascade complete)");
      resetMessage();
      return;
    }

    // SLOW PATH: popup is NOT yet visible (pendingShow may be true
    // because we are between the busy_show message and the show
    // timer firing) OR popup is in some inconsistent state. Schedule
    // an advisory hide on a IDLE_SETTLE_MS deadline so a "fast op"
    // that completes inside the show-debounce window still cleans
    // up if R never sends an explicit busy_hide. Does NOT reschedule
    // if a timer is already armed -- otherwise rapid busy/idle
    // toggling could push the deadline indefinitely.
    if ((overlayDomVisible() || isVisible || pendingShow) &&
        idleHideTimer === null) {
      idleHideTimer = setTimeout(function () {
        idleHideTimer = null;
        if (!pendingShow) {
          dbg("idle settle: advisory hide");
          clearPendingShow();
          forceHide("shiny:idle settle");
          resetMessage();
        }
      }, IDLE_SETTLE_MS);
    }
  }
  function overlayDomVisible() {
    var ov = document.getElementById("phenomapr-busy-overlay");
    return !!(ov && ov.classList.contains("is-visible"));
  }

  function onShinyDisconnected() {
    clearPendingShow();
    clearIdleHide();
    if (isVisible) renderHide();
    resetMessage();
    clearAllFileInputLoading();
    dismissedThisRun = false;
    dbg("disconnected: forced hide");
  }

  // ---- rAF watchdog --------------------------------------------------
  // Handles two cases that can otherwise leave the popup stuck on
  // screen:
  //
  //   1. STALE-HIDE: R sent busy_hide AFTER the most recent busy_show,
  //      but the popup is still visible (e.g. handleHide ran before
  //      renderShow could set isVisible, or a race re-added the
  //      visible class). If lastHideAt > lastShowAt and the popup is
  //      DOM-visible, force-hide on every tick.
  //
  //   2. ABSOLUTE: popup has been visible for >= ABSOLUTE_TIMEOUT
  //      (3 min) with no R-side intervention. Force-hide as a hard
  //      cap. Should essentially never fire.
  function watchdogTick() {
    var now = nowMs();
    var ov = document.getElementById("phenomapr-busy-overlay");
    var domVisible = !!(ov && ov.classList.contains("is-visible"));

    // Case 1: R has signalled hide more recently than show, but the
    // popup is still up. Honour the most recent R signal.
    if (domVisible && lastHideAt > 0 && lastHideAt >= lastShowAt) {
      console.warn("[phenomapr-busy] watchdog: stale-hide -> force-hiding");
      clearPendingShow();
      clearIdleHide();
      forceHide("watchdog stale-hide");
      resetMessage();
    }

    // Case 2: absolute hard cap.
    if (isVisible && shownAt > 0 && (now - shownAt) >= ABSOLUTE_TIMEOUT) {
      console.warn("[phenomapr-busy] watchdog: absolute timeout -> force-hiding");
      clearPendingShow();
      clearIdleHide();
      forceHide("watchdog absolute timeout");
      resetMessage();
    }

    if (typeof requestAnimationFrame === "function") {
      requestAnimationFrame(watchdogTick);
    } else {
      setTimeout(watchdogTick, 250);
    }
  }

  // ---- User-driven dismiss -------------------------------------------
  // The user can click the X (or backdrop) to dismiss the popup
  // mid-op. We mark `dismissedThisRun` so the next handleShow()
  // does not re-open it inside the same logical cycle. The flag is
  // cleared on the next shiny:idle or handleHide().
  function userDismiss(reason) {
    dbg("user dismissed", reason);
    clearPendingShow();
    clearIdleHide();
    if (isVisible) renderHide();
    dismissedThisRun = true;
  }

  // ---- R-side message API (drives visibility) -----------------------
  function handleShow(p) {
    console.log("[phenomapr-busy] handleShow received:",
                (p && p.message) || "(no message)");
    clearIdleHide();
    currentMessage = {
      message: (p && p.message) || "Working...",
      detail:  (p && p.detail)  || "",
      hint:    (p && p.hint)    || "This may take a moment..."
    };
    lastShowAt = nowMs();
    if (isVisible) {
      applyOverlayText();
      dbg("R-side show: live text update", currentMessage);
      return;
    }
    if (dismissedThisRun) {
      dbg("R-side show: suppressed (user dismissed earlier this cycle)");
      return;
    }
    // Schedule a debounced show. If R sends busy_hide before the
    // timer fires (i.e. the op turned out to be fast), we cancel.
    clearPendingShow();
    pendingShow = true;
    dbg("R-side show: queued, debouncing " + SHOW_DELAY_MS + "ms");
    pendingShowTimer = setTimeout(function () {
      pendingShowTimer = null;
      if (!pendingShow) return;        // already cancelled by hide
      pendingShow = false;
      if (dismissedThisRun) return;    // user dismissed during debounce
      dbg("R-side show: debounce elapsed, rendering popup");
      renderShow();
    }, SHOW_DELAY_MS);
  }
  // CRITICAL: handleHide MUST declare exactly one parameter, even
  // though we don\'t use it. Shiny.addCustomMessageHandler() does
  // `handler.length !== 1` and silently refuses to register a
  // zero-arg or multi-arg function with an "Uncaught handler must be
  // a function that takes one argument." error in the console.
  // Without that one parameter the entire phenomapr-busy-hide path
  // is dead -- popup never disappears AND dismissedThisRun stays
  // pinned to true after a manual dismiss, suppressing every
  // subsequent show. Do NOT remove the unused argument.
  function handleHide(_msg) {
    console.log("[phenomapr-busy] handleHide received");
    clearPendingShow();
    clearIdleHide();
    lastHideAt = nowMs();
    forceHide("R-side busy_hide");
    resetMessage();
    // R explicitly signalled end-of-op; allow next op to show again.
    dismissedThisRun = false;
  }

  // Robust hide: hides the popup unconditionally even if the JS
  // `isVisible` flag has somehow drifted out of sync with the DOM
  // (which can happen across hot reloads, partial renders, etc.).
  // Idempotent. Belt-and-suspenders -- both removes the
  // `.is-visible` class AND forces an inline `display: none` style,
  // so a CSS specificity issue (e.g. a rule with !important from
  // some downstream stylesheet) cannot leave the popup stuck on
  // screen. The inline style is cleared by renderShow() before
  // showing the popup again.
  function forceHide(reason) {
    var ov = document.getElementById("phenomapr-busy-overlay");
    if (ov) {
      var hadClass = ov.classList.contains("is-visible");
      ov.classList.remove("is-visible");
      ov.style.display = "none";
      if (hadClass) {
        dbg("forceHide(" + (reason || "") + "): removed .is-visible + inline display:none");
      }
    }
    isVisible = false;
    shownAt = 0;
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
  //   2. shiny:idle -> swap `.is-loading` for `.is-loading-complete`
  //      on every still-loading file-input.
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

  // ---- Wire up event listeners ---------------------------------------
  // We listen to BOTH native and jQuery flavours of shiny:busy /
  // shiny:idle because Shiny dispatches via jQuery internally, but
  // many Shiny+browser combos also bubble the events to document.
  // Both firing is harmless because the handlers are idempotent
  // (they short-circuit if the popup is already in the requested
  // state).
  if (document && document.addEventListener) {
    document.addEventListener("shiny:busy",          onShinyBusy,          false);
    document.addEventListener("shiny:idle",          onShinyIdle,          false);
    document.addEventListener("shiny:disconnected",  onShinyDisconnected,  false);
    document.addEventListener("shiny:inputchanged",  onInputChanged,       false);
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

    // Toggle the `open` attribute and an optional `mu-prompt`
    // class on a static HTML5 <details> element by id. Used to
    // auto-open the metadata-upload panel (and apply the salmon
    // prompt look) when the freshly-loaded expression file did
    // NOT carry its own cell-level metadata, and to auto-close +
    // un-prompt it once metadata is available -- without ever
    // re-rendering the <details> wrapper itself (which would tear
    // down the file picker inside it).
    //
    // msg = { id: "<element id>",
    //         open: true|false,
    //         prompt: true|false }
    Shiny.addCustomMessageHandler(
      "phenomapr-set-details-state",
      function(msg) {
        if (!msg || !msg.id) return;
        var el = document.getElementById(msg.id);
        if (!el) return;
        if (msg.open) {
          el.setAttribute("open", "open");
        } else {
          el.removeAttribute("open");
        }
        if (typeof msg.prompt !== "undefined") {
          if (msg.prompt) {
            el.classList.add("mu-prompt");
          } else {
            el.classList.remove("mu-prompt");
          }
        }
      }
    );

    var $j = window.jQuery || window.$;
    if ($j) {
      $j(document).on("shiny:busy",          onShinyBusy);
      $j(document).on("shiny:idle",          onShinyIdle);
      $j(document).on("shiny:disconnected",  onShinyDisconnected);
      $j(document).on("shiny:inputchanged",  onInputChanged);
    } else {
      dbg("jQuery not available at bind(); relying on native listeners");
    }

    // Start the absolute-timeout watchdog. Idempotent: bind() is
    // itself idempotent (called once at DOMContentLoaded), so this
    // starts exactly one rAF chain per page load.
    if (typeof requestAnimationFrame === "function") {
      requestAnimationFrame(watchdogTick);
    } else {
      setTimeout(watchdogTick, 1000);
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
