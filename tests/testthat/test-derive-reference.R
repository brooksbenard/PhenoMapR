# test-derive-reference.R
# Minimal tests for derive_reference_from_bulk (binary and continuous)

test_that("derive_reference_from_bulk returns data.frame for binary phenotype", {
  # Samples in rows, genes in columns
  set.seed(1)
  n_samp <- 20
  n_genes <- 10
  expr <- matrix(rnorm(n_samp * n_genes), nrow = n_samp, ncol = n_genes)
  rownames(expr) <- paste0("S", seq_len(n_samp))
  colnames(expr) <- paste0("G", seq_len(n_genes))
  pheno <- data.frame(
    sample_id = rownames(expr),
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
  expr <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("A", "B"), c("G1", "G2", "G3")))
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
  expr <- matrix(rnorm(n_samp * n_genes), nrow = n_samp, ncol = n_genes)
  rownames(expr) <- paste0("S", seq_len(n_samp))
  colnames(expr) <- paste0("G", seq_len(n_genes))
  pheno <- data.frame(
    sample_id = rownames(expr),
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

test_that("derive_reference_from_bulk with gene_axis rows (genes as rows) works", {
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
    gene_axis = "rows",
    verbose = FALSE
  )
  expect_s3_class(ref, "data.frame")
  expect_true(nrow(ref) >= 1)
})

test_that("derive_reference_from_bulk uses first phenotype column when sample_id_column NULL", {
  set.seed(4)
  expr <- matrix(rnorm(10 * 5), 10, 5)
  rownames(expr) <- paste0("S", 1:10)
  colnames(expr) <- paste0("G", 1:5)
  pheno <- data.frame(
    sample_id = rownames(expr),
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
