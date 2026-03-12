# test-utils.R
# Tests for list_cancer_types, get_celltype_palette, normalize_scores, get_gene_coverage

test_that("list_cancer_types returns a named list when reference is NULL", {
  out <- list_cancer_types()
  expect_type(out, "list")
  expect_named(out)
  expect_true(all(c("precog", "tcga", "pediatric_precog", "ici_precog") %in% names(out)))
})

test_that("list_cancer_types returns character vector for valid reference", {
  for (ref in c("precog", "tcga", "pediatric_precog", "ici_precog")) {
    out <- list_cancer_types(ref)
    expect_type(out, "character")
    expect_true(length(out) >= 1)
  }
})

test_that("list_cancer_types errors on invalid reference", {
  expect_error(list_cancer_types("invalid_ref"), "must be one of|precog|tcga")
})

test_that("get_celltype_palette returns named character vector", {
  types <- c("Acinar", "Ductal", "UnknownType")
  pal <- get_celltype_palette(types)
  expect_type(pal, "character")
  expect_named(pal)
  expect_equal(names(pal), types)
  expect_equal(pal["Acinar"], cell_type_colors["Acinar"])
  expect_equal(pal["Ductal"], cell_type_colors["Ductal"])
})

test_that("get_celltype_palette handles empty and NA", {
  pal_empty <- get_celltype_palette(character(0))
  expect_equal(length(pal_empty), 0)
  pal_na <- get_celltype_palette(c(NA, "Acinar"))
  expect_equal(names(pal_na), "Acinar")
})

test_that("normalize_scores returns z-scores", {
  x <- c(1, 2, 3, 4, 5)
  z <- normalize_scores(x)
  expect_equal(mean(z), 0, tolerance = 1e-10)
  expect_equal(sd(z), 1, tolerance = 1e-10)
  expect_length(z, length(x))
})

test_that("normalize_scores errors on non-numeric", {
  expect_error(normalize_scores(c("a", "b")), "scores must be numeric")
})

test_that("normalize_scores handles constant vector", {
  expect_warning(z <- normalize_scores(rep(5, 3)), "Standard deviation is 0")
  expect_equal(z, rep(5, 3))
})

test_that("get_gene_coverage returns data.frame for valid reference", {
  genes <- c("TP53", "MYC", "BRCA1")
  out <- get_gene_coverage(genes, "precog")
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) >= 1)
})

test_that("get_top_prognostic_genes returns data.frame", {
  out <- get_top_prognostic_genes("precog", "Breast", n = 5)
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) <= 5)
})

test_that("get_top_prognostic_genes direction positive and negative", {
  out_pos <- get_top_prognostic_genes("precog", "Breast", n = 3, direction = "positive")
  expect_s3_class(out_pos, "data.frame")
  expect_true(all(out_pos$z_score > 0))
  out_neg <- get_top_prognostic_genes("precog", "Breast", n = 3, direction = "negative")
  expect_true(all(out_neg$z_score < 0))
})

test_that("get_gene_coverage with reference NULL checks all references", {
  genes <- c("TP53", "MYC", "BRCA1")
  out <- get_gene_coverage(genes, reference = NULL)
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) >= 4)  # precog, tcga, pediatric_precog, ici_precog
  expect_true("reference" %in% names(out))
})

test_that("get_gene_coverage with empty genes returns coverage_pct NA", {
  out <- get_gene_coverage(character(0), reference = "precog")
  expect_s3_class(out, "data.frame")
  expect_equal(out$total_genes, 0)
  expect_equal(out$covered_genes, 0)
  expect_true(is.na(out$coverage_pct))
})

test_that("list_cancer_types accepts pediatric and ici aliases", {
  peds <- list_cancer_types("pediatric")
  expect_type(peds, "character")
  expect_true(length(peds) >= 1)
  ici <- list_cancer_types("ici")
  expect_type(ici, "character")
})

test_that("load_rds_fast errors on missing file", {
  expect_error(load_rds_fast("/nonexistent/path/to/file.rds"), "File not found")
})

test_that("plot_score_distribution uses first column when score_column NULL and data.frame", {
  df <- data.frame(a = 1:5, b = 10:14)
  p <- plot_score_distribution(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_distribution handles all NA scores", {
  df <- data.frame(score = rep(NA_real_, 5))
  p <- plot_score_distribution(df)
  expect_s3_class(p, "ggplot")
})
