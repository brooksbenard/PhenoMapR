# Single-cell integration: CRA001160 (PAAD scRNA-seq)

root <- .integration_state$config$root
helpers <- integration_source_shiny_helpers()
demo_path <- integration_ensure_path("shiny_demo", download = FALSE, root = root)
seurat_path <- integration_resolve_path("cra001160_rds", root = root)

load_cra_subset <- function(n_cells = 800L, seed = 42L) {
  if (!is.na(seurat_path) && isTRUE(integration_get_option("full"))) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat required for full CRA001160 RDS test")
    }
    obj <- readRDS(seurat_path)
    obj <- tryCatch(Seurat::UpdateSeuratObject(obj), error = function(e) obj)
    cells <- colnames(obj)
    set.seed(seed)
    keep <- sample(cells, min(n_cells, length(cells)))
    return(list(
      expression = obj[, keep],
      kind = "seurat",
      metadata = obj@meta.data[keep, , drop = FALSE]
    ))
  }
  if (!is.na(demo_path) && file.exists(demo_path)) {
    bundle <- readRDS(demo_path)
    expr <- bundle$expression
    md <- bundle$metadata
    if (!"cell_type" %in% names(md) && "cell_type_minor" %in% names(md)) {
      md$cell_type <- md$cell_type_minor
    }
    set.seed(seed)
    cells <- colnames(expr)
    keep <- sample(cells, min(as.integer(n_cells), length(cells)))
    expr <- expr[, keep, drop = FALSE]
    md <- md[match(keep, md$.cell_id), , drop = FALSE]
    if (exists(".attach_demo_matrix_obsm", envir = helpers, inherits = FALSE)) {
      expr <- helpers$.attach_demo_matrix_obsm(expr, md)
    }
    return(list(expression = expr, metadata = md, kind = "matrix"))
  }
  Sys.setenv(PHENOMAPR_SHINY_DEMO_RDS = demo_path)
  if (exists("clear_shiny_demo_pool_cache", envir = helpers, inherits = FALSE)) {
    helpers$clear_shiny_demo_pool_cache()
  }
  demo <- helpers$make_shiny_demo_dataset(
    n_genes = 300L,
    n_cells = min(n_cells, 5000L),
    seed = seed
  )
  list(expression = demo$expression, metadata = demo$metadata, kind = "matrix")
}

run_marker_scope <- function(expr, groups, scope_name) {
  scopes <- list(
    cohort = list(marker_scope = "phenotype_groups", celltype_contrast = "within_cell_type"),
    vs_opposite = list(marker_scope = "cell_type_specific", celltype_contrast = "vs_opposite_tail"),
    within_ct = list(marker_scope = "cell_type_specific", celltype_contrast = "within_cell_type")
  )
  cfg <- scopes[[scope_name]]
  ct_col <- if ("cell_type" %in% colnames(groups)) "cell_type" else NULL
  score_col <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
  find_phenotype_markers(
    expression = expr,
    group_labels = groups,
    group_column = score_col,
    cell_id_column = "cell_id",
    cell_type_column = ct_col,
    marker_scope = cfg$marker_scope,
    celltype_contrast = cfg$celltype_contrast,
    verbose = FALSE,
    max_cells_per_ident = 2000L
  )
}

integration_run_test(
  "sc:cra001160:demo_available",
  "Shiny demo bundle or CRA001160 Seurat RDS available",
  skip_if = if (is.na(demo_path) && is.na(seurat_path)) {
    "No CRA001160 demo bundle or Seurat RDS found"
  } else NULL,
  expr = list(demo = demo_path, seurat = seurat_path)
)

if (!is.na(demo_path) || !is.na(seurat_path)) {
  dat <- load_cra_subset(n_cells = if (integration_get_option("full")) 2000L else 1200L)

  integration_run_test(
    "sc:cra001160:phenomap",
    "PhenoMap scores and phenotype tail directionality",
    expr = {
      scores <- PhenoMap(
        expression = dat$expression,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      integration_assert(is.data.frame(scores), "scores not data.frame")
      sc <- scores[[grep("Pancreatic", colnames(scores), value = TRUE)[1]]]
      names(sc) <- rownames(scores)
      groups <- define_phenotype_groups(scores, percentile = 0.10)
      gcol <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
      adverse <- groups$cell_id[groups[[gcol]] == "Most Adverse"]
      favorable <- groups$cell_id[groups[[gcol]] == "Most Favorable"]
      integration_assert(length(adverse) > 0 && length(favorable) > 0, "empty tails")
      mean_adv <- mean(sc[adverse], na.rm = TRUE)
      mean_fav <- mean(sc[favorable], na.rm = TRUE)
      integration_assert(
        mean_adv > mean_fav,
        sprintf("adverse mean %.4f should exceed favorable %.4f", mean_adv, mean_fav)
      )
      n <- nrow(groups)
      expect_adv <- ceiling(0.10 * n)
      integration_assert(
        abs(length(adverse) - expect_adv) <= 2L,
        sprintf("adverse tail size %d vs expected ~%d", length(adverse), expect_adv)
      )
      if ("cell_type" %in% colnames(dat$metadata)) {
        groups$cell_type <- dat$metadata$cell_type[match(groups$cell_id, dat$metadata$.cell_id %||% dat$metadata$cell_id)]
      }
      list(
        n_cells = n,
        mean_adverse = mean_adv,
        mean_favorable = mean_fav,
        groups = groups,
        scores = scores
      )
    }
  )

  # Re-load for marker tests (avoid depending on previous test return value)
  dat <- load_cra_subset(n_cells = if (integration_get_option("full")) 2000L else 1200L)
  scores <- PhenoMap(
    expression = dat$expression,
    reference = "precog",
    cancer_type = "Pancreatic",
    verbose = FALSE
  )
  groups <- define_phenotype_groups(scores, percentile = 0.10)
  if ("cell_type" %in% colnames(dat$metadata)) {
    id_col <- if (".cell_id" %in% colnames(dat$metadata)) ".cell_id" else "cell_id"
    groups$cell_type <- dat$metadata$cell_type[match(groups$cell_id, dat$metadata[[id_col]])]
  }

  for (scope in c("cohort", "vs_opposite", "within_ct")) {
    integration_run_test(
      paste0("sc:cra001160:markers:", scope),
      paste("find_phenotype_markers scope:", scope),
      expr = {
        mk <- run_marker_scope(dat$expression, groups, scope)
        integration_assert(is.list(mk), "markers not a list")
        if (scope == "within_ct" && nrow(mk$adverse_markers) == 0L) {
          integration_assert(
            any(!is.na(groups$cell_type)),
            "within_ct empty: no cell types in groups"
          )
          return(list(skipped = "no within-cell-type markers on this subset"))
        }
        integration_assert(nrow(mk$adverse_markers) > 0, "no adverse markers")
        integration_assert(nrow(mk$favorable_markers) > 0, "no favorable markers")
        if ("avg_log2FC" %in% names(mk$adverse_markers)) {
          integration_assert(
            mean(mk$adverse_markers$avg_log2FC > 0, na.rm = TRUE) >= 0.5,
            "majority of adverse markers should have positive avg_log2FC"
          )
        }
        top_adv <- integration_top_genes(mk$adverse_markers, 20L)
        golden <- integration_check_golden(
          paste0("cra001160_markers_", scope, ".rds"),
          top_adv,
          threshold = 0.8
        )
        list(
          n_adverse = nrow(mk$adverse_markers),
          n_favorable = nrow(mk$favorable_markers),
          golden_overlap = golden$overlap
        )
      }
    )
  }

  integration_run_test(
    "sc:cra001160:seurat_crosscheck",
    "Top markers overlap Seurat::FindMarkers on fixed subset",
    skip_if = if (!requireNamespace("Seurat", quietly = TRUE)) "Seurat not installed" else NULL,
    expr = {
      set.seed(99)
      expr_use <- dat$expression
      if (inherits(expr_use, "Seurat")) {
        cells <- colnames(expr_use)
        keep <- sample(cells, min(500L, length(cells)))
        expr_use <- expr_use[, keep]
        gsub <- groups[groups$cell_id %in% keep, , drop = FALSE]
        obj <- expr_use
        Seurat::Idents(obj) <- gsub[[grep("^phenotype_group_", colnames(gsub), value = TRUE)[1]]][
          match(colnames(obj), gsub$cell_id)
        ]
      } else {
        mat <- expr_use
        cells <- colnames(mat)
        keep <- sample(cells, min(500L, length(cells)))
        mat <- mat[, keep, drop = FALSE]
        gsub <- groups[groups$cell_id %in% keep, , drop = FALSE]
        obj <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA")
        obj <- Seurat::NormalizeData(obj, verbose = FALSE)
        gcol <- grep("^phenotype_group_", colnames(gsub), value = TRUE)[1]
        Seurat::Idents(obj) <- gsub[[gcol]][match(colnames(obj), gsub$cell_id)]
        expr_use <- mat
        groups_use <- gsub
      }
      seu_mk <- Seurat::FindMarkers(
        obj,
        ident.1 = "Most Adverse",
        ident.2 = "Most Favorable",
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        verbose = FALSE
      )
      seu_top <- rownames(seu_mk[order(seu_mk$p_val), , drop = FALSE])[seq_len(min(20L, nrow(seu_mk)))]
      pm_mk <- find_phenotype_markers(
        expression = expr_use,
        group_labels = if (exists("groups_use")) groups_use else groups,
        group_column = grep("^phenotype_group_", colnames(groups), value = TRUE)[1],
        cell_id_column = "cell_id",
        marker_scope = "phenotype_groups",
        verbose = FALSE,
        max_cells_per_ident = 2000L
      )
      pm_top <- integration_top_genes(pm_mk$adverse_markers, 20L)
      overlap <- integration_jaccard(seu_top, pm_top)
      integration_assert(overlap >= 0.25, sprintf("Seurat overlap %.3f < 0.25", overlap))
      list(overlap = overlap)
    }
  )

  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {
    integration_run_test(
      "sc:cra001160:heatmap",
      "plot_phenotype_markers produces PNG output",
      expr = {
        mk <- run_marker_scope(dat$expression, groups, "cohort")
        score_col <- grep("Pancreatic", colnames(scores), value = TRUE)[1]
        gcol <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
        meta <- data.frame(
          cell_id = groups$cell_id,
          stringsAsFactors = FALSE
        )
        meta[[score_col]] <- scores[[score_col]][match(groups$cell_id, rownames(scores))]
        meta[[gcol]] <- groups[[gcol]]
        if ("cell_type" %in% colnames(groups)) meta$cell_type <- groups$cell_type
        tmp <- tempfile(fileext = ".png")
        grDevices::png(tmp, width = 800, height = 600)
        plot_phenotype_markers(
          markers = mk,
          expr_mat = dat$expression,
          meta = meta,
          cell_id_col = "cell_id",
          group_col = gcol,
          score_col = score_col,
          celltype_col = if ("cell_type" %in% names(meta)) "cell_type" else NULL,
          top_n_markers = 10L,
          draw = TRUE
        )
        if (grDevices::dev.cur() > 1L) grDevices::dev.off()
        integration_assert(file.exists(tmp), "heatmap PNG not created")
        integration_assert(file.info(tmp)$size > 5000L, "heatmap PNG too small")
        list(bytes = file.info(tmp)$size)
      }
    )
  }
}
