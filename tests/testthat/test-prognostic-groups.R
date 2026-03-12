# test-prognostic-groups.R
# Tests for define_prognostic_groups

test_that("define_prognostic_groups returns data.frame with group column", {
  scores_df <- data.frame(
    sample = paste0("S", 1:20),
    score = c(rep(1, 5), rep(2, 10), rep(3, 5)),
    stringsAsFactors = FALSE
  )
  rownames(scores_df) <- scores_df$sample
  scores_df$sample <- NULL
  groups <- define_prognostic_groups(scores_df, percentile = 0.2)
  expect_s3_class(groups, "data.frame")
  expect_true("cell_id" %in% names(groups))
  expect_true(any(grepl("prognostic_group", names(groups))))
  expect_equal(sort(unique(groups[[2]])), c("Most Adverse", "Most Favorable", "Other"))
})

test_that("define_prognostic_groups errors on non-data.frame", {
  expect_error(define_prognostic_groups(c(1, 2, 3)), "must be a data.frame")
})

test_that("define_prognostic_groups errors when percentile out of range", {
  scores_df <- data.frame(x = 1:5)
  expect_error(define_prognostic_groups(scores_df, percentile = 0), "between 0 and 0.5")
  expect_error(define_prognostic_groups(scores_df, percentile = 0.5), "between 0 and 0.5")
})

test_that("define_prognostic_groups uses score_columns when provided", {
  scores_df <- data.frame(
    id = paste0("C", 1:10),
    sc1 = rnorm(10),
    sc2 = rnorm(10),
    stringsAsFactors = FALSE
  )
  rownames(scores_df) <- scores_df$id
  groups <- define_prognostic_groups(scores_df, percentile = 0.1, score_columns = "sc1")
  expect_true("prognostic_group_sc1" %in% names(groups))
  expect_false("prognostic_group_sc2" %in% names(groups))
})

test_that("define_prognostic_groups errors when no numeric score columns", {
  scores_df <- data.frame(id = letters[1:5], label = c("A", "B", "A", "B", "A"))
  expect_error(
    define_prognostic_groups(scores_df, percentile = 0.1),
    "No numeric score columns"
  )
})

test_that("define_prognostic_groups generates cell_id when rownames NULL", {
  scores_df <- data.frame(score = rnorm(6))
  rownames(scores_df) <- NULL
  groups <- define_prognostic_groups(scores_df, percentile = 0.2)
  expect_true("cell_id" %in% names(groups))
  expect_equal(groups$cell_id, as.character(1:6))
})

test_that("define_prognostic_groups skips missing score_columns with warning", {
  scores_df <- data.frame(sc1 = rnorm(8))
  rownames(scores_df) <- paste0("C", 1:8)
  expect_warning(
    groups <- define_prognostic_groups(scores_df, percentile = 0.1, score_columns = c("sc1", "nonexistent")),
    "not found"
  )
  expect_true("prognostic_group_sc1" %in% names(groups))
})

test_that("define_prognostic_groups handles all-NA score column", {
  scores_df <- data.frame(score = rep(NA_real_, 5))
  rownames(scores_df) <- paste0("C", 1:5)
  groups <- define_prognostic_groups(scores_df, percentile = 0.2)
  expect_true("prognostic_group_score" %in% names(groups))
  expect_true(all(is.na(groups$prognostic_group_score)))
})
