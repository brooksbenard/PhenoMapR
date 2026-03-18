# test-derive-reference.R
# Minimal tests for derive_reference_from_bulk (binary and continuous)

test_that("derive_reference_from_bulk returns data.frame for binary phenotype", {
  # Genes (rows) x samples (columns)
  set.seed(1)
  n_samp <- 20
  n_genes <- 10
  expr <- matrix(rnorm(n_genes * n_samp), nrow = n_genes, ncol = n_samp)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samp))
  pheno <- data.frame(
    sample_id = colnames(expr),
    response = rep(c("R", "NR"), each = 10),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = "sample_id",
    phenotype_column = "response",
    phenotype_type = "binary",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
  expect_true(nrow(ref) >= 1)
  expect_true(ncol(ref) >= 1)
})

test_that("derive_reference_from_bulk errors when no sample overlap", {
  # Genes x samples; colnames = sample IDs. No overlap with pheno.
  expr <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("G1", "G2"), c("A", "B", "C")))
  pheno <- data.frame(sample_id = c("X", "Y"), y = c(1, 0))
  expect_error(
    derive_reference_from_bulk(expr, pheno, sample_id_column = "sample_id",
                               phenotype_column = "y", phenotype_type = "binary", verbose = FALSE),
    "No overlapping sample IDs"
  )
})

test_that("derive_reference_from_bulk errors on non-matrix bulk_expression", {
  pheno <- data.frame(id = "S1", y = 1)
  expect_error(
    derive_reference_from_bulk("not_a_matrix", pheno, phenotype_type = "binary", verbose = FALSE),
    "must be a matrix"
  )
})

test_that("derive_reference_from_bulk continuous phenotype returns data.frame", {
  set.seed(2)
  n_samp <- 25
  n_genes <- 15
  expr <- matrix(rnorm(n_genes * n_samp), nrow = n_genes, ncol = n_samp)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samp))
  pheno <- data.frame(
    sample_id = colnames(expr),
    continuous_y = rnorm(n_samp),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = "sample_id",
    phenotype_column = "continuous_y",
    phenotype_type = "continuous",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
  expect_true(nrow(ref) >= 1)
  expect_true(ncol(ref) >= 1)
})

test_that("derive_reference_from_bulk accepts genes as rows and samples as columns", {
  set.seed(3)
  n_genes <- 8
  n_samp <- 12
  expr <- matrix(rnorm(n_genes * n_samp), nrow = n_genes, ncol = n_samp)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samp))
  pheno <- data.frame(
    id = colnames(expr),
    y = rep(c(0, 1), length.out = n_samp),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = "id",
    phenotype_column = "y",
    phenotype_type = "binary",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
  expect_true(nrow(ref) >= 1)
})

test_that("derive_reference_from_bulk uses first phenotype column when sample_id_column NULL", {
  set.seed(4)
  # Genes x samples (enough genes so heuristic does not transpose)
  expr <- matrix(rnorm(100 * 10), 100, 10)
  rownames(expr) <- paste0("G", seq_len(100))
  colnames(expr) <- paste0("S", 1:10)
  pheno <- data.frame(
    sample_id = colnames(expr),
    response = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = NULL,
    phenotype_column = "response",
    phenotype_type = "binary",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
})

test_that("derive_reference_from_bulk survival phenotype returns data.frame", {
  skip_if_not_installed("survival")
  set.seed(5)
  n_samp <- 30
  n_genes <- 20
  expr <- matrix(rnorm(n_genes * n_samp), nrow = n_genes, ncol = n_samp)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samp))
  pheno <- data.frame(
    sample_id = colnames(expr),
    time = runif(n_samp, 1, 10),
    event = sample(c(0, 1), n_samp, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = "sample_id",
    phenotype_type = "survival",
    survival_time = "time",
    survival_event = "event",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
  expect_true(nrow(ref) >= 1)
})

test_that("derive_reference_from_bulk phenotype_type auto infers binary", {
  set.seed(6)
  # Genes (rows) x samples (columns)
  expr <- matrix(rnorm(100 * 15), 100, 15)
  rownames(expr) <- paste0("G", seq_len(100))
  colnames(expr) <- paste0("S", seq_len(15))
  pheno <- data.frame(
    sample_id = colnames(expr),
    y = rep(c(0, 1), length.out = 15),
    stringsAsFactors = FALSE
  )
  ref <- derive_reference_from_bulk(
    bulk_expression = expr,
    phenotype = pheno,
    sample_id_column = "sample_id",
    phenotype_column = "y",
    phenotype_type = "auto",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
})

test_that("derive_reference_from_bulk errors when phenotype has no rows", {
  # Genes x samples
  expr <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("G1", "G2"), c("S1", "S2", "S3")))
  pheno <- data.frame(sample_id = character(0), y = numeric(0))
  expect_error(
    derive_reference_from_bulk(expr, pheno, phenotype_column = "y", phenotype_type = "binary", verbose = FALSE),
    "no rows"
  )
})

test_that("derive_reference_from_bulk errors when sample_id_column not found", {
  # Genes x samples (enough genes so no transpose)
  expr <- matrix(rnorm(20 * 2), 20, 2, dimnames = list(paste0("G", 1:20), c("S1", "S2")))
  pheno <- data.frame(id = c("S1", "S2"), y = c(0, 1))
  expect_error(
    derive_reference_from_bulk(expr, pheno, sample_id_column = "wrong_col", phenotype_column = "y", phenotype_type = "binary", verbose = FALSE),
    "Sample ID column|not found"
  )
})

test_that("plot_reference_signature runs with derive_reference output", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  ref <- data.frame(z = rnorm(60), row.names = paste0("G", 1:60))
  colnames(ref)[1] <- "Survival z-score"
  pdf(NULL)
  out <- plot_reference_signature(ref, n_label = 5L, row_title = "Survival z-score")
  dev.off()
  expect_true(inherits(out, "HeatmapList"))
})

