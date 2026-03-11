# test-plot-score-distribution.R
# Tests for plot_score_distribution

test_that("plot_score_distribution returns a ggplot", {
  scores_vec <- c(1, 2, 2, 3, 3, 3, 4, 5)
  p <- plot_score_distribution(scores_vec)
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_distribution accepts data.frame and score_column", {
  scores_df <- data.frame(
    cell_id = paste0("C", 1:10),
    my_score = rnorm(10),
    stringsAsFactors = FALSE
  )
  p <- plot_score_distribution(scores_df, score_column = "my_score")
  expect_s3_class(p, "ggplot")
})
