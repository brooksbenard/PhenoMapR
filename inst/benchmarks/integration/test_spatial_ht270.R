# Spatial integration: HT270 CytoSPACE pancreatic Visium deconvolution

root <- .integration_state$config$root
dl <- integration_get_option("download")
cyto_path <- integration_ensure_path("spatial_cyto", download = dl, root = root)

integration_run_test(
  "spatial:ht270:path",
  "HT270 CytoSPACE processed RDS available",
  skip_if = if (is.na(cyto_path)) "HT270 processed RDS not found" else NULL,
  expr = list(path = cyto_path)
)

if (!is.na(cyto_path)) {
  helpers <- integration_source_shiny_helpers()

  integration_run_test(
    "spatial:ht270:phenomap",
    "PhenoMap on spatial/cell Seurat object",
    expr = {
      obj <- readRDS(cyto_path)
      integration_assert(!is.null(obj), "failed to read RDS")
      scores <- PhenoMap(
        expression = obj,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      sc <- scores[[grep("Pancreatic", colnames(scores), value = TRUE)[1]]]
      integration_assert(all(is.finite(sc)), "non-finite spatial scores")
      groups <- define_phenotype_groups(scores, percentile = 0.05)
      list(n_cells = nrow(scores), mean_score = mean(sc, na.rm = TRUE))
    }
  )

  integration_run_test(
    "spatial:ht270:embeddings",
    "Spatial and segmentation reductions detected",
    expr = {
      obj <- readRDS(cyto_path)
      emb <- helpers$list_available_embeddings(obj)
      integration_assert(length(emb) >= 1L, "no embeddings detected")
      has_spatial <- any(tolower(emb) %in% c("spatial", "umap", "segmentation"))
      integration_assert(has_spatial, "expected spatial/UMAP-like embedding")
      list(embeddings = paste(emb, collapse = ", "))
    }
  )

  integration_run_test(
    "spatial:ht270:coloc_tables",
    "Precomputed spatial colocalization tables load",
    expr = {
      data_dir <- file.path(root, "inst", "data")
      files <- c(
        "spatial_coloc_cellchat_pairs.rds",
        "spatial_coloc_cellchat_lr_pairs.rds",
        "spatial_coloc_cellchat_assoc.rds",
        "spatial_pair_cell_metrics.rds"
      )
      dims <- vapply(files, function(f) {
        fp <- file.path(data_dir, f)
        integration_assert(file.exists(fp), paste("missing", f))
        obj <- readRDS(fp)
        n <- if (is.data.frame(obj)) nrow(obj) else length(obj)
        if (f == "spatial_pair_cell_metrics.rds" && n == 0L) return(0)
        integration_assert(n > 0L, paste(f, "is empty"))
        n
      }, numeric(1))
      list(rows = paste(names(dims), dims, sep = "=", collapse = "; "))
    }
  )
}
