# AnnData / h5ad integration: CosMx demo (python-free hdf5r path)

root <- .integration_state$config$root
h5ad_path <- integration_ensure_path("cosmx_h5ad", download = FALSE, root = root)
helpers <- integration_source_shiny_helpers()

integration_run_test(
  "h5ad:cosmx:path",
  "CosMx h5ad file available",
  skip_if = if (is.na(h5ad_path)) {
    "CosMx h5ad not found (set PHENOMAPR_COSMX_H5AD or place in ~/Downloads)"
  } else NULL,
  expr = list(path = h5ad_path)
)

if (!is.na(h5ad_path)) {
  integration_run_test(
    "h5ad:cosmx:parse",
    "parse_expression_upload preserves matrix and phenomapr_obs",
    skip_if = if (!requireNamespace("hdf5r", quietly = TRUE)) "hdf5r not installed" else NULL,
    expr = {
      parsed <- helpers$parse_expression_upload(h5ad_path, basename(h5ad_path))
      obj <- parsed$object
      integration_assert(!is.null(obj), "parse returned NULL object")
      integration_assert(is.matrix(obj) || inherits(obj, "Matrix"), "not a matrix")
      obs <- attr(obj, "phenomapr_obs")
      integration_assert(is.data.frame(obs), "phenomapr_obs missing")
      integration_assert(ncol(obj) == nrow(obs), "cell count mismatch")
      integration_assert(any(c("fov", "x_FOV_px", "y_FOV_px") %in% names(obs)),
                         "expected CosMx obs columns")
      list(n_cells = ncol(obj), n_genes = nrow(obj), obs_cols = ncol(obs))
    }
  )

  integration_run_test(
    "h5ad:cosmx:phenomap",
    "PhenoMap on h5ad-derived matrix",
    expr = {
      parsed <- helpers$parse_expression_upload(h5ad_path, basename(h5ad_path))
      obj <- parsed$object
      scores <- PhenoMap(
        expression = obj,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      sc <- scores[[1]]
      integration_assert(all(is.finite(sc)), "non-finite scores")
      list(n = length(sc))
    }
  )

  integration_run_test(
    "h5ad:cosmx:markers",
    "find_phenotype_markers on h5ad-derived matrix",
    expr = {
      parsed <- helpers$parse_expression_upload(h5ad_path, basename(h5ad_path))
      obj <- parsed$object
      scores <- PhenoMap(expression = obj, reference = "precog",
                         cancer_type = "Pancreatic", verbose = FALSE)
      groups <- define_phenotype_groups(scores, percentile = 0.10)
      gcol <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
      mk <- find_phenotype_markers(
        expression = obj,
        group_labels = groups,
        group_column = gcol,
        cell_id_column = "cell_id",
        marker_scope = "phenotype_groups",
        verbose = FALSE,
        max_cells_per_ident = 1500L
      )
      integration_assert(nrow(mk$adverse_markers) > 0, "no adverse markers")
      list(n_adverse = nrow(mk$adverse_markers))
    }
  )

  integration_run_test(
    "h5ad:cosmx:polygons",
    "No false polygon path when boundaries absent",
    expr = {
      parsed <- helpers$parse_expression_upload(h5ad_path, basename(h5ad_path))
      obj <- parsed$object
      has_poly <- helpers$spatial_polygons_available(obj)
      integration_assert(!isTRUE(has_poly), "polygons should not be available for CosMx demo")
      emb <- helpers$list_available_embeddings(obj)
      list(embeddings = paste(emb, collapse = ", "), polygons = has_poly)
    }
  )

  integration_run_test(
    "h5ad:cosmx:combined_coords",
    "Global spatial coords for combined sections view",
    expr = {
      parsed <- helpers$parse_expression_upload(h5ad_path, basename(h5ad_path))
      obj <- parsed$object
      seg_emb <- helpers$extract_embedding(obj, "segmentation")
      integration_assert(!is.null(seg_emb) && nrow(seg_emb) > 0, "no segmentation embedding")
      combined_emb <- helpers$apply_global_spatial_coords_for_combined(obj, seg_emb)
      integration_assert(nrow(combined_emb) == ncol(obj), "coord row count mismatch")
      list(n = nrow(combined_emb), fov_local = helpers$.spatial_coords_are_fov_local(seg_emb))
    }
  )
}
