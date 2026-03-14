# test-phenomap-seurat-sce.R
# PhenoMap with Seurat and SingleCellExperiment (data_handlers + scoring paths)

test_that("PhenoMap with Seurat object returns scores", {
  skip_if_not_installed("Seurat")
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(25, nrow(precog)))]
  n_cells <- 6
  counts <- matrix(
    pmax(0, rnorm(length(genes) * n_cells)),
    nrow = length(genes),
    ncol = n_cells,
    dimnames = list(genes, paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  scores <- PhenoMap(expression = obj, reference = "precog", cancer_type = "Breast", assay = "RNA", slot = "data", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_cells)
  expect_true(any(grepl("weighted_sum", colnames(scores))))
})

test_that("PhenoMap with SingleCellExperiment returns scores", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(25, nrow(precog)))]
  n_cells <- 5
  mat <- matrix(
    pmax(0, rnorm(length(genes) * n_cells)),
    nrow = length(genes),
    ncol = n_cells,
    dimnames = list(genes, paste0("C", seq_len(n_cells)))
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat, logcounts = log1p(mat))
  )
  scores <- PhenoMap(expression = sce, reference = "precog", cancer_type = "Breast", assay = "logcounts", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_cells)
})
