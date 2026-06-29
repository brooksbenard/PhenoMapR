source_shiny_helpers <- function() {
  src_helpers <- tryCatch(
    testthat::test_path("..", "..", "inst", "shiny", "helpers.R"),
    error = function(e) NULL
  )
  helpers <- if (!is.null(src_helpers) && file.exists(src_helpers)) {
    src_helpers
  } else {
    system.file("shiny", "helpers.R", package = "PhenoMapR")
  }
  testthat::skip_if(!nzchar(helpers) || !file.exists(helpers),
                    "shiny/helpers.R not found")
  env <- new.env(parent = globalenv())
  sys.source(helpers, envir = env)
  if (exists("clear_shiny_demo_pool_cache", envir = env, inherits = FALSE)) {
    env$clear_shiny_demo_pool_cache()
  }
  pkg_root <- normalizePath(file.path(dirname(helpers), "..", ".."), winslash = "/")
  bundle <- file.path(pkg_root, "inst", "extdata", "shiny", "PAAD_CRA001160_demo_5000.rds")
  if (file.exists(bundle)) {
    Sys.setenv(PHENOMAPR_SHINY_DEMO_RDS = bundle)
    env$.shiny_demo_bundle_override <- bundle
  }
  env
}

test_that("make_shiny_demo_dataset returns expression with UMAP metadata", {
  e <- source_shiny_helpers()
  bundle <- e$.shiny_demo_bundle_override %||% e$.shiny_demo_bundle_path()
  has_bundle <- !is.null(bundle) && nzchar(bundle) && file.exists(bundle)
  demo <- e$make_shiny_demo_dataset(n_genes = 200L, n_cells = 80L, seed = 42L)
  expect_equal(ncol(demo$expression), if (has_bundle) 5000L else 80L)
  expect_lte(nrow(demo$expression), if (has_bundle) 500L else 200L)
  expect_equal(nrow(demo$metadata), ncol(demo$expression))
  expect_equal(length(demo$metadata$UMAP_1), ncol(demo$expression))
  expect_equal(length(demo$metadata$UMAP_2), ncol(demo$expression))
  expect_true(all(c("UMAP_1", "UMAP_2", "cell_type", "Source") %in% colnames(demo$metadata)))
  expect_true(!is.null(demo$source_info))
  expect_true(!is.null(demo$source_info$accession))
})

test_that("make_shiny_demo_dataset draws different cells without a fixed seed", {
  e <- source_shiny_helpers()
  bundle <- e$.shiny_demo_bundle_override %||% e$.shiny_demo_bundle_path()
  has_bundle <- !is.null(bundle) && nzchar(bundle) && file.exists(bundle)
  d1 <- e$make_shiny_demo_dataset(n_genes = 100L, n_cells = 50L, seed = 1L)
  d2 <- e$make_shiny_demo_dataset(n_genes = 100L, n_cells = 50L, seed = 2L)
  if (!has_bundle) {
    expect_false(identical(colnames(d1$expression), colnames(d2$expression)))
  } else {
    expect_equal(colnames(d1$expression), colnames(d2$expression))
  }
})

test_that("make_shiny_demo_dataset exposes UMAP on the expression matrix", {
  e <- source_shiny_helpers()
  demo <- e$make_shiny_demo_dataset(n_genes = 100L, n_cells = 40L, seed = 42L)
  emb <- e$list_available_embeddings(demo$expression)
  expect_true(length(emb) >= 1L)
  expect_true("UMAP" %in% emb)
})

test_that("make_shiny_demo_dataset attaches colData for matrix pseudobulk", {
  e <- source_shiny_helpers()
  demo <- e$make_shiny_demo_dataset(n_genes = 100L, n_cells = 40L, seed = 42L)
  coldata <- attr(demo$expression, "phenomapr_coldata")
  expect_true(is.data.frame(coldata))
  expect_equal(nrow(coldata), ncol(demo$expression))
  expect_true(all(c("Patient", "Source") %in% colnames(coldata)))
  expect_equal(rownames(coldata), colnames(demo$expression))

  genes <- rownames(demo$expression)
  ref <- data.frame(
    row.names = genes[seq_len(min(5L, length(genes)))],
    s = rep(1, min(5L, length(genes)))
  )
  scores <- PhenoMapR::PhenoMap(
    expression = demo$expression,
    reference = ref,
    pseudobulk = TRUE,
    group_by = "Patient",
    verbose = FALSE
  )
  n_patients <- length(unique(coldata$Patient[!is.na(coldata$Patient)]))
  expect_equal(nrow(scores), n_patients)
  expect_lt(nrow(scores), ncol(demo$expression))
})

test_that("build_cell_table honors score_column selection", {
  e <- source_shiny_helpers()
  scores <- data.frame(
    score_a = c(1, 2, 3),
    score_b = c(10, 20, 30),
    row.names = c("C1", "C2", "C3")
  )
  tbl <- e$build_cell_table(scores, score_column = "score_b")
  expect_equal(tbl$score, c(10, 20, 30))
  expect_equal(attr(tbl, "score_name"), "score_b")
})

test_that("resolve_phenotype_group_column matches score-specific group column", {
  e <- source_shiny_helpers()
  groups <- data.frame(
    cell_id = c("C1", "C2"),
    phenotype_group_score_a = c("Most Adverse", "Other"),
    phenotype_group_score_b = c("Other", "Most Favorable"),
    stringsAsFactors = FALSE
  )
  expect_equal(
    e$resolve_phenotype_group_column(groups, "score_b"),
    "phenotype_group_score_b"
  )
})
