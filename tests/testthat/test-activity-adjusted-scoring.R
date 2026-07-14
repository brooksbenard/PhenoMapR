test_that("PhenoMap activity_adjusted mode runs on matrix input", {
  skip_if_not_installed("Seurat")
  custom_ref <- data.frame(
    row.names = paste0("G", 1:60),
    s = c(rep(3, 30), rep(-3, 30))
  )
  expr <- matrix(
    pmax(0, rpois(60 * 40, lambda = 5)),
    nrow = 60,
    ncol = 40,
    dimnames = list(paste0("G", 1:60), paste0("C", 1:40))
  )
  md <- data.frame(
    cell_id = colnames(expr),
    nCount_RNA = colSums(expr),
    stringsAsFactors = FALSE
  )
  expect_warning(
    scores <- PhenoMap(
      expression = expr,
      reference = custom_ref,
      z_score_cutoff = 2,
      score_mode = "activity_adjusted",
      vars_to_regress = "nCount_RNA",
      sample_metadata = md,
      verbose = FALSE
    ),
    "negative values"
  )
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 40)
})

test_that("PhenoMap permutation_n adds empirical p column", {
  custom_ref <- data.frame(
    row.names = paste0("G", 1:60),
    s = c(rep(3, 30), rep(-3, 30))
  )
  expr <- matrix(
    pmax(0, rpois(60 * 30, lambda = 5)),
    nrow = 60,
    ncol = 30,
    dimnames = list(paste0("G", 1:60), paste0("C", 1:30))
  )
  scores <- suppressWarnings(PhenoMap(
      expression = expr,
      reference = custom_ref,
      z_score_cutoff = 2,
      permutation_n = 10L,
      permutation_seed = 1L,
      verbose = FALSE
    ))
  expect_true(any(grepl("^empirical_p_", colnames(scores))))
  expect_true(all(scores[[grep("^empirical_p_", colnames(scores), value = TRUE)[1]]] > 0))
})
