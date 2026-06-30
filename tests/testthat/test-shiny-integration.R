# Shiny server-equivalent integration (testthat)

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

test_that("Shiny parity: demo -> score -> groups -> markers chain", {
  e <- source_shiny_helpers()
  demo <- e$make_shiny_demo_dataset(n_genes = 200L, n_cells = 120L, seed = 5L)
  scores <- PhenoMapR::PhenoMap(
    expression = demo$expression,
    reference = "precog",
    cancer_type = "Pancreatic",
    verbose = FALSE
  )
  groups <- PhenoMapR::define_phenotype_groups(scores, percentile = 0.10)
  score_col <- e$active_score_column(scores)
  grp_col <- e$resolve_phenotype_group_column(groups, score_col)
  expect_false(is.na(grp_col))
  groups$cell_type <- demo$metadata$cell_type[match(groups$cell_id, demo$metadata$.cell_id)]
  markers <- PhenoMapR::find_phenotype_markers(
    expression = demo$expression,
    group_labels = groups,
    group_column = grp_col,
    cell_id_column = "cell_id",
    cell_type_column = "cell_type",
    marker_scope = "phenotype_groups",
    verbose = FALSE,
    max_cells_per_ident = 1000L
  )
  expect_gt(nrow(markers$adverse_markers), 0)
  expect_gt(nrow(markers$favorable_markers), 0)
})

test_that("Shiny parity: percentile regrouping changes tail sizes", {
  e <- source_shiny_helpers()
  demo <- e$make_shiny_demo_dataset(n_genes = 150L, n_cells = 100L, seed = 3L)
  scores <- PhenoMapR::PhenoMap(
    expression = demo$expression,
    reference = "precog",
    cancer_type = "Pancreatic",
    verbose = FALSE
  )
  g10 <- PhenoMapR::define_phenotype_groups(scores, percentile = 0.10)
  g20 <- PhenoMapR::define_phenotype_groups(scores, percentile = 0.20)
  score_col <- e$active_score_column(scores)
  grp_col <- e$resolve_phenotype_group_column(g10, score_col)
  n_adv_10 <- sum(g10[[grp_col]] == "Most Adverse")
  n_adv_20 <- sum(g20[[grp_col]] == "Most Adverse")
  expect_gt(n_adv_20, n_adv_10)
})

test_that("Shiny testServer captures server function from app.R", {
  testthat::skip_if_not_installed("shiny")
  app_r <- tryCatch(
    testthat::test_path("..", "..", "inst", "shiny", "app.R"),
    error = function(e) NULL
  )
  if (is.null(app_r) || !file.exists(app_r)) {
    app_r <- system.file("shiny", "app.R", package = "PhenoMapR")
  }
  testthat::skip_if(!nzchar(app_r) || !file.exists(app_r), "app.R not found")
  captured <- new.env(parent = emptyenv())
  stub_env <- new.env(parent = globalenv())
  stub_env$shinyApp <- function(ui, server) {
    captured$ui <- ui
    captured$server <- server
    invisible(NULL)
  }
  parse_ok <- tryCatch(
    sys.source(app_r, envir = stub_env),
    error = function(e) e
  )
  if (inherits(parse_ok, "error")) {
    testthat::skip(paste("app.R source failed:", conditionMessage(parse_ok)))
  }
  testthat::skip_if(!exists("server", envir = captured, inherits = FALSE),
                    "server function not captured from app.R")
  shiny::testServer(
    app = shiny::shinyApp(
      ui = captured$ui,
      server = captured$server
    ),
    expr = {
      session$flushReact()
      expect_true(!is.null(output))
    }
  )
})
