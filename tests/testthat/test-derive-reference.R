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
