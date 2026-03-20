# test-weighted-sum.R
# Tests that hit calculate_weighted_scores / compute_scores paths (via PhenoMap)
# and normalize_scores edge cases

test_that("PhenoMap with few overlapping genes warns", {
  custom_ref <- data.frame(
    row.names = paste0("G", 1:15),
    s = c(rep(3, 5), rep(-3, 5), rep(1, 5))
  )
  expr <- matrix(
    pmax(0, rnorm(15 * 4)),
    nrow = 15,
    ncol = 4,
    dimnames = list(paste0("G", 1:15), paste0("C", 1:4))
  )
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE),
    "Fewer than 50 overlapping"
  )
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 4)
})

test_that("PhenoMap reference_sign = -1 negates scores relative to reference_sign = 1", {
  custom_ref <- data.frame(
    row.names = paste0("G", 1:20),
    s = c(rep(3, 10), rep(-3, 10))
  )
  expr <- matrix(
    pmax(0, rnorm(20 * 5)),
    nrow = 20,
    ncol = 5,
    dimnames = list(paste0("G", 1:20), paste0("C", 1:5))
  )
  expect_warning(
    s1 <- PhenoMap(
      expression = expr,
      reference = custom_ref,
      z_score_cutoff = 2,
      verbose = FALSE,
      reference_sign = 1L
    ),
    "Fewer than 50 overlapping"
  )
  expect_warning(
    s2 <- PhenoMap(
      expression = expr,
      reference = custom_ref,
      z_score_cutoff = 2,
      verbose = FALSE,
      reference_sign = -1L
    ),
    "Fewer than 50 overlapping"
  )
  expect_equal(s1[[1]], -s2[[1]])
})

test_that("normalize_scores with NA in input", {
  x <- c(1, 2, NA, 4, 5)
  z <- normalize_scores(x)
  expect_length(z, 5)
  expect_true(is.na(z[3]))
  expect_equal(mean(z[-3], na.rm = TRUE), 0, tolerance = 1e-10)
  expect_equal(sd(z[-3], na.rm = TRUE), 1, tolerance = 1e-10)
})

test_that("normalize_scores with constant vector warns", {
  expect_warning(z <- normalize_scores(rep(5, 3)), "Standard deviation is 0")
  expect_equal(z, rep(5, 3))
})

test_that("calculate_weighted_scores errors when reference_data is not data.frame or matrix", {
  expr <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("C1", "C2")))
  expect_error(
    PhenoMapR:::calculate_weighted_scores(expr, reference_data = list(), verbose = FALSE),
    "reference_data must be a data.frame or matrix"
  )
})

test_that("compute_scores errors when expression_data has no rownames", {
  expr <- matrix(1:4, 2, 2)
  scores <- setNames(c(1, -1), c("A", "B"))
  expect_error(
    PhenoMapR:::compute_scores(expr, scores, pseudobulk = FALSE, verbose = FALSE),
    "rownames"
  )
})

test_that("compute_scores errors when prognostic_scores has no names", {
  expr <- matrix(1:4, 2, 2, dimnames = list(c("A", "B"), c("C1", "C2")))
  scores <- c(1, -1)
  expect_error(
    PhenoMapR:::compute_scores(expr, scores, pseudobulk = FALSE, verbose = FALSE),
    "named"
  )
})
