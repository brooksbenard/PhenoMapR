# test-phenomap.R
# Tests for PhenoMap() with matrix and custom reference

test_that("PhenoMap with matrix and precog returns data.frame of scores", {
  # Use genes that exist in precog; precog uses "Breast" not "BRCA"
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(30, nrow(precog)))]
  n_samp <- 5
  expr <- matrix(
    pmax(0, rnorm(length(genes) * n_samp)),  # non-negative to avoid warning
    nrow = length(genes),
    ncol = n_samp,
    dimnames = list(genes, paste0("S", seq_len(n_samp)))
  )
  scores <- PhenoMap(expression = expr, reference = "precog", cancer_type = "Breast", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_samp)
  expect_true(ncol(scores) >= 1)
  expect_equal(rownames(scores), paste0("S", seq_len(n_samp)))
})

test_that("PhenoMap with custom reference data.frame works", {
  custom_ref <- data.frame(
    row.names = c("TP53", "MYC", "EGFR"),
    my_sig = c(2.5, -1.8, 2.0)
  )
  expr <- matrix(
    pmax(0, rnorm(3 * 4)),
    nrow = 3,
    ncol = 4,
    dimnames = list(c("TP53", "MYC", "EGFR"), paste0("C", 1:4))
  )
  scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 4)
  expect_true(any(grepl("my_sig|weighted_sum", colnames(scores))))
})

test_that("PhenoMap errors on invalid reference name", {
  expr <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("C1", "C2")))
  expect_error(
    PhenoMap(expression = expr, reference = "not_a_ref", cancer_type = "Breast"),
    "Unknown reference|must be one of"
  )
})

test_that("PhenoMap errors when cancer_type missing for built-in reference", {
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[1:5]
  expr <- matrix(pmax(0, rnorm(5 * 2)), 5, 2, dimnames = list(genes, c("A", "B")))
  expect_error(
    PhenoMap(expression = expr, reference = "precog", cancer_type = NULL),
    "cancer_type"
  )
})

test_that("PhenoMap errors on non-matrix expression", {
  custom_ref <- data.frame(row.names = "G1", sig = 1)
  expect_error(
    PhenoMap(expression = "not_a_matrix", reference = custom_ref),
    "Unable to detect input type|Supported"
  )
})

test_that("PhenoMap accepts data.frame expression (process_matrix path)", {
  custom_ref <- data.frame(row.names = c("G1", "G2"), sig = c(1, -1))
  expr_df <- as.data.frame(matrix(pmax(0, rnorm(2 * 4)), 2, 4))
  rownames(expr_df) <- c("G1", "G2")
  colnames(expr_df) <- paste0("C", 1:4)
  scores <- PhenoMap(expression = expr_df, reference = custom_ref, verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 4)
})

test_that("PhenoMap with matrix lacking colnames generates cell names and returns scores", {
  custom_ref <- data.frame(row.names = c("A", "B"), s = c(3, -2))
  expr <- matrix(1, nrow = 2, ncol = 3)
  rownames(expr) <- c("A", "B")
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE),
    "no column names"
  )
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 3)
  expect_true(ncol(scores) >= 1)
})

