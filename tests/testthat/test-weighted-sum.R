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
