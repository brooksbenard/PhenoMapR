# test-add-scores.R
# Tests for add_scores_to_seurat and add_scores_to_sce (conditional on packages)

test_that("add_scores_to_seurat adds score columns and returns Seurat object", {
  skip_if_not_installed("Seurat")
  scores <- data.frame(
    row.names = c("c1", "c2", "c3"),
    weighted_sum_score_sig = c(0.5, -0.2, 0.8)
  )
  counts <- matrix(1, nrow = 10, ncol = 3, dimnames = list(paste0("G", 1:10), c("c1", "c2", "c3")))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  out <- add_scores_to_seurat(obj, scores)
  expect_true(inherits(out, "Seurat"))
  expect_true("weighted_sum_score_sig" %in% names(out@meta.data))
  expect_equal(out@meta.data["c1", "weighted_sum_score_sig"], 0.5)
})

test_that("add_scores_to_seurat errors when no common cells", {
  skip_if_not_installed("Seurat")
  scores <- data.frame(row.names = c("x", "y"), s = c(1, 2))
  counts <- matrix(1, nrow = 5, ncol = 2, dimnames = list(paste0("G", 1:5), c("c1", "c2")))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  expect_error(add_scores_to_seurat(obj, scores), "No common cells")
})

test_that("add_scores_to_seurat warns when some cells in scores not in object", {
  skip_if_not_installed("Seurat")
  scores <- data.frame(
    row.names = c("c1", "c2", "c3", "c4"),
    s = c(1, 2, 3, 4)
  )
  counts <- matrix(1, nrow = 5, ncol = 2, dimnames = list(paste0("G", 1:5), c("c1", "c2")))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  expect_warning(out <- add_scores_to_seurat(obj, scores), "not found")
  expect_equal(ncol(out), 2)
})

test_that("add_scores_to_seurat errors on non-Seurat object", {
  expect_error(add_scores_to_seurat(list(), data.frame(row.names = "a", x = 1)), "Seurat object")
})

test_that("add_scores_to_sce adds score columns to colData", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  scores <- data.frame(
    row.names = c("c1", "c2"),
    score_sig = c(0.1, -0.3)
  )
  mat <- matrix(1, nrow = 10, ncol = 2, dimnames = list(paste0("G", 1:10), c("c1", "c2")))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat),
    colData = S4Vectors::DataFrame(cell = c("c1", "c2"))
  )
  out <- add_scores_to_sce(sce, scores)
  expect_true(inherits(out, "SingleCellExperiment"))
  expect_true("score_sig" %in% names(SummarizedExperiment::colData(out)))
})

test_that("add_scores_to_sce errors when no common cells", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  scores <- data.frame(row.names = c("x", "y"), s = c(1, 2))
  mat <- matrix(1, nrow = 5, ncol = 2, dimnames = list(paste0("G", 1:5), c("c1", "c2")))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
  expect_error(add_scores_to_sce(sce, scores), "No common cells")
})

test_that("add_scores_to_sce errors on non-SCE object", {
  expect_error(add_scores_to_sce(list(), data.frame(row.names = "a", x = 1)), "SingleCellExperiment")
})

test_that("add_scores_to_seurat with prefix adds prefixed column names", {
  skip_if_not_installed("Seurat")
  scores <- data.frame(row.names = c("c1", "c2"), s = c(1, 2))
  counts <- matrix(1, nrow = 5, ncol = 2, dimnames = list(paste0("G", 1:5), c("c1", "c2")))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  out <- add_scores_to_seurat(obj, scores, prefix = "pref_")
  expect_true("pref_s" %in% names(out@meta.data))
})
