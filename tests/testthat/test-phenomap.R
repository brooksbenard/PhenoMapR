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

test_that("PhenoMap accepts dgCMatrix (e.g. from Read10X_h5)", {
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(20, nrow(precog)))]
  n_samp <- 4
  # Sparse matrix with same structure as Read10X_h5 output
  expr_dense <- matrix(
    pmax(0, rnorm(length(genes) * n_samp)),
    nrow = length(genes),
    ncol = n_samp,
    dimnames = list(genes, paste0("C", seq_len(n_samp)))
  )
  expr_sparse <- Matrix::Matrix(expr_dense, sparse = TRUE)
  expect_s4_class(expr_sparse, "dgCMatrix")
  scores <- PhenoMap(expression = expr_sparse, reference = "precog", cancer_type = "Breast", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_samp)
  expect_equal(rownames(scores), paste0("C", seq_len(n_samp)))
})

test_that("PhenoMap with TCGA reference returns scores", {
  data(tcga, package = "PhenoMapR", envir = environment())
  genes <- rownames(tcga)[seq_len(min(25, nrow(tcga)))]
  if (!"PAAD" %in% colnames(tcga)) skip("TCGA PAAD not in tcga data")
  expr <- matrix(
    pmax(0, rnorm(length(genes) * 4)),
    nrow = length(genes),
    ncol = 4,
    dimnames = list(genes, paste0("S", 1:4))
  )
  scores <- PhenoMap(expression = expr, reference = "tcga", cancer_type = "PAAD", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 4)
  expect_true(any(grepl("PAAD|weighted_sum", colnames(scores))))
})

test_that("PhenoMap verbose=TRUE runs without error", {
  custom_ref <- data.frame(row.names = c("G1", "G2", "G3"), s = c(2, -1, 2))
  expr <- matrix(pmax(0, rnorm(3 * 3)), 3, 3, dimnames = list(c("G1", "G2", "G3"), c("C1", "C2", "C3")))
  capture.output(scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = TRUE))
  expect_s3_class(scores, "data.frame")
})

test_that("PhenoMap with matrix containing NA warns", {
  custom_ref <- data.frame(row.names = c("A", "B"), s = c(1, -1))
  expr <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("A", "B"), c("C1", "C2", "C3")))
  expr[1, 1] <- NA_real_
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE),
    "NA"
  )
  expect_s3_class(scores, "data.frame")
})

test_that("PhenoMap with matrix containing negative values warns", {
  custom_ref <- data.frame(row.names = c("A", "B"), s = c(1, -1))
  expr <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("A", "B"), c("C1", "C2", "C3")))
  expr[1, 1] <- -0.5
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE),
    "negative"
  )
  expect_s3_class(scores, "data.frame")
})

test_that("PhenoMap with custom ref and no genes passing z cutoff warns", {
  custom_ref <- data.frame(row.names = c("G1", "G2"), s = c(0.5, -0.5))
  expr <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("G1", "G2"), c("C1", "C2")))
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, z_score_cutoff = 10, verbose = FALSE),
    "No genes pass|z-score cutoff"
  )
  expect_s3_class(scores, "data.frame")
})

test_that("PhenoMap with custom ref and no common genes warns", {
  custom_ref <- data.frame(row.names = c("X", "Y", "Z"), s = c(3, -2, 2))
  expr <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("C1", "C2")))
  expect_warning(
    scores <- PhenoMap(expression = expr, reference = custom_ref, verbose = FALSE),
    "No common genes"
  )
  expect_s3_class(scores, "data.frame")
})

test_that("PhenoMap errors when pseudobulk TRUE but group_by NULL", {
  custom_ref <- data.frame(row.names = c("G1"), s = 1)
  expr <- matrix(1, nrow = 1, ncol = 2, dimnames = list("G1", c("C1", "C2")))
  expect_error(
    PhenoMap(expression = expr, reference = custom_ref, pseudobulk = TRUE, verbose = FALSE),
    "group_by"
  )
})

test_that("PhenoMap with Seurat slot counts uses counts layer", {
  skip_if_not_installed("Seurat")
  data(precog, package = "PhenoMapR", envir = environment())
  genes <- rownames(precog)[seq_len(min(30, nrow(precog)))]
  n_cells <- 5
  counts <- matrix(
    pmax(0, as.integer(rnorm(length(genes) * n_cells, 50, 20))),
    nrow = length(genes),
    ncol = n_cells,
    dimnames = list(genes, paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  scores <- PhenoMap(expression = obj, reference = "precog", cancer_type = "Breast", assay = "RNA", slot = "counts", verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), n_cells)
})

test_that("PhenoMap with ici_precog returns scores when column matches", {
  data(ici, package = "PhenoMapR", envir = environment())
  if (ncol(ici) == 0) skip("ICI data has no columns")
  cols <- colnames(ici)
  metastatic_col <- cols[grepl("_Metastatic_", cols)][1]
  if (is.na(metastatic_col)) skip("No ICI metastatic column found")
  cancer_type <- sub("_.+", "", metastatic_col)
  cancer_type <- paste0(cancer_type, "_Metastatic")
  genes <- rownames(ici)[seq_len(min(25, nrow(ici)))]
  expr <- matrix(
    pmax(0, rnorm(length(genes) * 4)),
    nrow = length(genes),
    ncol = 4,
    dimnames = list(genes, paste0("S", 1:4))
  )
  scores <- PhenoMap(expression = expr, reference = "ici_precog", cancer_type = cancer_type, verbose = FALSE)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 4)
})

