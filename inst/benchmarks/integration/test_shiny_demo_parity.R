# Shiny demo bundle API parity (mirrors server workflow without browser)

root <- .integration_state$config$root
helpers <- integration_source_shiny_helpers()
demo_path <- integration_ensure_path("shiny_demo", download = FALSE, root = root)

integration_run_test(
  "shiny:demo:bundle",
  "make_shiny_demo_dataset matches bundled RDS dimensions",
  skip_if = if (is.na(demo_path)) "Shiny demo RDS not found" else NULL,
  expr = {
    bundle <- readRDS(demo_path)
    integration_assert(ncol(bundle$expression) == 5000L, "expected 5000-cell demo bundle")
    set.seed(11L)
    keep <- sample(colnames(bundle$expression), 500L)
    expr <- bundle$expression[, keep, drop = FALSE]
    md <- bundle$metadata[match(keep, bundle$metadata$.cell_id), , drop = FALSE]
    demo <- list(expression = expr, metadata = md)
    integration_assert(ncol(demo$expression) == 500L, "subset cell count")
    integration_assert(nrow(demo$metadata) == ncol(demo$expression), "metadata row mismatch")
    tbl <- helpers$build_cell_table(
      scores = NULL,
      groups = NULL,
      metadata = demo$metadata
    )
    integration_assert(is.null(tbl) || nrow(tbl) == ncol(demo$expression),
                       "build_cell_table nrow mismatch without scores")
    list(n_cells = ncol(demo$expression), n_genes = nrow(demo$expression))
  }
)

if (!is.na(demo_path)) {
  integration_run_test(
    "shiny:demo:full_chain",
    "Demo expression -> PhenoMap -> groups -> markers (Shiny server parity)",
    expr = {
      Sys.setenv(PHENOMAPR_SHINY_DEMO_RDS = demo_path)
      helpers$clear_shiny_demo_pool_cache()
      demo <- helpers$make_shiny_demo_dataset(n_genes = 250L, n_cells = 400L, seed = 21L)
      scores <- PhenoMap(
        expression = demo$expression,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      groups <- define_phenotype_groups(scores, percentile = 0.10)
      score_col <- helpers$active_score_column(scores, score_col_input = NULL)
      grp_col <- helpers$resolve_phenotype_group_column(groups, score_col)
      integration_assert(!is.na(grp_col), "could not resolve phenotype group column")
      groups$cell_type <- demo$metadata$cell_type[match(groups$cell_id, demo$metadata$.cell_id)]
      mk <- find_phenotype_markers(
        expression = demo$expression,
        group_labels = groups,
        group_column = grp_col,
        cell_id_column = "cell_id",
        cell_type_column = "cell_type",
        marker_scope = "phenotype_groups",
        verbose = FALSE,
        max_cells_per_ident = 2000L
      )
      integration_assert(nrow(mk$adverse_markers) > 0, "no adverse markers in demo chain")
      cell_tbl <- helpers$build_cell_table(scores, groups, demo$metadata, score_column = score_col)
      integration_assert("group" %in% colnames(cell_tbl), "cell table missing group column")
      list(
        n_markers = nrow(mk$adverse_markers),
        grp_col = grp_col,
        score_col = score_col
      )
    }
  )
}
