# test-find-prognostic-markers.R
# Tests for find_phenotype_markers (matrix path and group_label resolution)

test_that("find_phenotype_markers with matrix and group vector returns list", {
  set.seed(1)
  n_genes <- 50
  n_cells <- 60
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  group_vec <- c(
    rep("Most Adverse", 10),
    rep("Most Favorable", 10),
    rep("Other", 40)
  )
  out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_s3_class(out$adverse_markers, "data.frame")
  expect_s3_class(out$favorable_markers, "data.frame")
  expect_true("gene" %in% names(out$adverse_markers))
  expect_true("p_val" %in% names(out$adverse_markers))
})

test_that("find_phenotype_markers with verbose=TRUE runs (matrix path)", {
  set.seed(99)
  expr <- matrix(pmax(0, rnorm(40 * 45)), 40, 45, dimnames = list(paste0("G", 1:40), paste0("C", 1:45)))
  group_vec <- c(rep("Most Adverse", 8), rep("Most Favorable", 8), rep("Other", 29))
  msg <- capture.output(
    out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = TRUE),
    type = "message"
  )
  expect_type(out, "list")
  expect_true(length(msg) >= 0)
})

test_that("find_phenotype_markers with group_labels data.frame and group_column", {
  set.seed(2)
  n_genes <- 40
  n_cells <- 50
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  groups_df <- data.frame(
    cell_id = paste0("C", seq_len(n_cells)),
    pg = c(rep("Most Adverse", 8), rep("Most Favorable", 8), rep("Other", 34)),
    stringsAsFactors = FALSE
  )
  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers supports cell_type_specific marker detection", {
  set.seed(6)
  n_genes <- 20
  n_cells <- 12
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  # Two cell types with exactly 2 adverse, 2 favorable, 2 other per type
  cell_type <- c(rep("T", 6), rep("B", 6))
  pg <- c(
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2),
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2)
  )
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = pg,
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )

  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
  expect_true(all(unique(out$adverse_markers$cell_type) %in% c("T", "B")))
})

test_that("find_phenotype_markers cell_type_specific handles empty per-cell-type results", {
  set.seed(66)
  n_genes <- 30
  n_cells <- 24
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  # T cells have adverse/favorable/other; B cells have only "Other"
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = c(
      rep("Most Adverse", 4), rep("Most Favorable", 4), rep("Other", 4),
      rep("Other", 12)
    ),
    cell_type = c(rep("T", 12), rep("B", 12)),
    stringsAsFactors = FALSE
  )

  expect_no_error({
    out <- find_phenotype_markers(
      expr,
      group_labels = groups_df,
      group_column = "pg",
      cell_id_column = "cell_id",
      marker_scope = "cell_type_specific",
      cell_type_column = "cell_type",
      min.pct = 0,
      logfc.threshold = 0,
      pval_threshold = 1,
      max_cells_per_ident = Inf,
      verbose = FALSE
    )
  })

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
  expect_true("cell_type" %in% names(out$favorable_markers))
})

test_that("find_phenotype_markers errors when group_labels length mismatch", {
  expr <- matrix(1, nrow = 5, ncol = 10, dimnames = list(paste0("G", 1:5), paste0("C", 1:10)))
  group_vec <- rep("Other", 5)
  expect_error(
    find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE),
    "Length of .* must match"
  )
})

test_that("find_phenotype_markers warns when few adverse or favorable cells", {
  set.seed(3)
  expr <- matrix(pmax(0, rnorm(30 * 20)), 30, 20, dimnames = list(paste0("G", 1:30), paste0("C", 1:20)))
  group_vec <- c(rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 16))
  expect_warning(
    out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE),
    "Fewer than 5"
  )
  expect_type(out, "list")
})

test_that("find_phenotype_markers accepts dgCMatrix", {
  skip_if_not_installed("Matrix")
  set.seed(4)
  n_genes <- 30
  n_cells <- 40
  expr_dense <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  expr_sparse <- Matrix::Matrix(expr_dense, sparse = TRUE)
  group_vec <- c(rep("Most Adverse", 6), rep("Most Favorable", 6), rep("Other", 28))
  out <- find_phenotype_markers(expr_sparse, group_labels = group_vec, verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers with max_cells_per_ident subsamples", {
  set.seed(5)
  n_genes <- 25
  n_cells <- 200
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  group_vec <- c(
    rep("Most Adverse", 80),
    rep("Most Favorable", 80),
    rep("Other", 40)
  )
  out <- find_phenotype_markers(
    expr,
    group_labels = group_vec,
    max_cells_per_ident = 20L,
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers with invalid group labels returns empty markers", {
  expr <- matrix(1, nrow = 5, ncol = 5, dimnames = list(paste0("G", 1:5), paste0("C", 1:5)))
  bad_df <- data.frame(id = paste0("C", 1:5), grp = c("A", "B", "C", "D", "E"))
  out <- find_phenotype_markers(expr, group_labels = bad_df, group_column = "grp", verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_equal(nrow(out$adverse_markers), 0)
  expect_equal(nrow(out$favorable_markers), 0)
})

test_that("find_phenotype_markers with Seurat object uses FindMarkers path", {
  skip_if_not_installed("Seurat")
  set.seed(42)
  n_genes <- 40
  n_cells <- 45
  counts <- matrix(
    pmax(0, as.integer(rnorm(n_genes * n_cells, 10, 3))),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  groups <- c(
    rep("Most Adverse", 8),
    rep("Most Favorable", 8),
    rep("Other", 29)
  )
  out <- find_phenotype_markers(obj, group_labels = groups, assay = "RNA", slot = "data", verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_s3_class(out$adverse_markers, "data.frame")
  expect_s3_class(out$favorable_markers, "data.frame")
})
