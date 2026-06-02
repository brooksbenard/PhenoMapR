test_that("detect_expression_format flags bulk raw counts with HUGO IDs", {
  set.seed(42)
  m <- matrix(rpois(200 * 30, 5), nrow = 200,
              dimnames = list(c("TP53", "MYC", "EGFR", "BRCA1", "CDKN2A",
                                "KRAS", "PIK3CA", "PTEN", "MKI67", "CD68",
                                paste0("GENE", 1:190)),
                              paste0("S", 1:30)))
  d <- detect_expression_format(m)
  expect_equal(d$format, "raw_counts")
  expect_equal(d$gene_id_kind, "hugo")
  expect_equal(d$sc_or_bulk, "bulk")
  expect_true(any(grepl("Raw counts detected", d$recommendations)))
})


test_that("detect_expression_format flags single-cell raw counts via sparsity", {
  set.seed(7)
  sc <- matrix(0, nrow = 200, ncol = 400,
               dimnames = list(c("TP53", "MYC", paste0("Gene", 1:198)),
                               paste0("Cell_", 1:400)))
  nz <- sample(length(sc), 0.2 * length(sc))
  sc[nz] <- rpois(length(nz), 2)
  d <- detect_expression_format(sc)
  expect_equal(d$format, "raw_counts")
  expect_equal(d$sc_or_bulk, "single_cell")
  expect_true(d$sparsity >= 0.5)
  # Single-cell modality => label is denominated in CELLS, not samples.
  expect_match(d$sc_or_bulk_label, "Single-cell-like\\b")
  expect_match(d$sc_or_bulk_label, "\\bcells\\b")
  expect_false(grepl("\\bsamples\\b", d$sc_or_bulk_label))
})


test_that("detect_expression_format keeps 'samples' wording for bulk-like data", {
  set.seed(13)
  ## Dense, low-sparsity matrix with few columns -> bulk.
  m <- matrix(runif(200 * 30, min = 0.5, max = 12), nrow = 200,
              dimnames = list(paste0("G", 1:200), paste0("S", 1:30)))
  d <- detect_expression_format(m)
  expect_equal(d$sc_or_bulk, "bulk")
  expect_match(d$sc_or_bulk_label, "Bulk-like\\b")
  expect_match(d$sc_or_bulk_label, "\\bsamples\\b")
  expect_false(grepl("\\bcells\\b", d$sc_or_bulk_label))
})


test_that("detect_expression_format flags CPM and log-normalized matrices", {
  set.seed(3)
  m <- matrix(rpois(200 * 30, 5), nrow = 200,
              dimnames = list(paste0("G", 1:200), paste0("S", 1:30)))
  cpm <- sweep(m, 2, pmax(colSums(m), 1), "/") * 1e6
  d_cpm <- detect_expression_format(cpm)
  expect_equal(d_cpm$format, "cpm_or_tpm")

  logcpm <- log2(cpm + 1)
  d_log <- detect_expression_format(logcpm)
  expect_equal(d_log$format, "log_normalized")
})


test_that("detect_expression_format flags z-scaled matrices", {
  set.seed(11)
  m <- matrix(rnorm(50 * 20), nrow = 50,
              dimnames = list(paste0("G", 1:50), paste0("S", 1:20)))
  d <- detect_expression_format(m)
  expect_equal(d$format, "z_scaled")
})


test_that("detect_expression_format flags ENSG IDs and mixed gene-ID styles", {
  set.seed(2)
  m <- matrix(rpois(200 * 30, 5), nrow = 200,
              dimnames = list(sprintf("ENSG%011d", 1:200),
                              paste0("S", 1:30)))
  d_e <- detect_expression_format(m)
  expect_equal(d_e$gene_id_kind, "ensembl")
  expect_true(any(grepl("Ensembl", d_e$recommendations)))

  # 100 ENSG + 100 HUGO -> neither group exceeds 50%, so classified mixed
  rownames(m)[1:100] <- paste0("HUGO", 1:100)
  d_mix <- detect_expression_format(m)
  expect_equal(d_mix$gene_id_kind, "mixed")
})


test_that("detect_expression_format finds duplicate gene IDs", {
  m <- matrix(rpois(100 * 10, 3), nrow = 100,
              dimnames = list(c("TP53", "TP53", "MYC", "MYC", "MYC",
                                paste0("G", 1:95)),
                              paste0("S", 1:10)))
  d <- detect_expression_format(m)
  expect_equal(d$n_dup, 2L)
  expect_true("TP53" %in% d$dup_examples || "MYC" %in% d$dup_examples)
})


test_that("clean_matrix_input log2(CPM+1) normalizes bulk raw counts", {
  set.seed(4)
  m <- matrix(rpois(200 * 30, 5), nrow = 200,
              dimnames = list(paste0("G", 1:200), paste0("S", 1:30)))
  res <- clean_matrix_input(m, do_hugo = FALSE, mode = "bulk",
                            verbose = FALSE)
  expect_lt(max(res$matrix), 25)
  expect_gte(min(res$matrix), 0)
  expect_true(any(grepl("log2", res$steps)))
})


test_that("clean_matrix_input applies Seurat-style normalization for sc data", {
  set.seed(5)
  sc <- matrix(0, nrow = 200, ncol = 400,
               dimnames = list(paste0("G", 1:200), paste0("C", 1:400)))
  nz <- sample(length(sc), 0.2 * length(sc))
  sc[nz] <- rpois(length(nz), 2)
  res <- clean_matrix_input(sc, do_hugo = FALSE, mode = "single_cell",
                            verbose = FALSE)
  expect_true(any(grepl("Seurat", res$steps)))
  expect_lt(max(res$matrix), 15)
  expect_gte(min(res$matrix), 0)
})


test_that("clean_matrix_input collapses duplicate gene rows by mean", {
  set.seed(6)
  m <- matrix(rpois(100 * 10, 3), nrow = 100,
              dimnames = list(c("TP53", "TP53", "TP53", "MYC", "MYC",
                                paste0("G", 1:95)),
                              paste0("S", 1:10)))
  res <- clean_matrix_input(m, do_hugo = FALSE, do_normalize = FALSE,
                            mode = "bulk", verbose = FALSE)
  expect_equal(res$n_collapsed, 3L)
  expect_false(anyDuplicated(rownames(res$matrix)) > 0)
})


test_that("clean_matrix_input skips normalization for already-log data", {
  set.seed(8)
  m <- matrix(runif(200 * 30, 0, 8), nrow = 200,
              dimnames = list(paste0("G", 1:200), paste0("S", 1:30)))
  res <- clean_matrix_input(m, do_hugo = FALSE, mode = "bulk",
                            verbose = FALSE)
  expect_true(any(grepl("already log-normalized", res$steps)))
  expect_equal(dim(res$matrix), dim(m))
})


test_that("clean_matrix_input errors on rownameless input", {
  m <- matrix(rpois(100, 3), 10, 10)
  expect_error(clean_matrix_input(m), "rownames")
})


test_that("detect_expression_format handles 10X-scale dgCMatrix without warning or slowness", {
  skip_if_not_installed("Matrix")
  set.seed(7)
  n_genes <- 5000L; n_cells <- 8000L  # smaller than 30k x 50k for speed
  target_nnz <- as.integer(0.08 * n_genes * n_cells)
  i <- sample.int(n_genes, target_nnz, replace = TRUE)
  j <- sample.int(n_cells, target_nnz, replace = TRUE)
  vals <- as.integer(rpois(target_nnz, lambda = 3) + 1L)
  m <- Matrix::sparseMatrix(i = i, j = j, x = vals,
    dims = c(n_genes, n_cells),
    dimnames = list(c("TP53", "MYC", paste0("Gene", 1:(n_genes - 2))),
                    paste0("Cell_", 1:n_cells)),
    check = FALSE)

  warns <- character(0)
  t0 <- Sys.time()
  d <- withCallingHandlers(
    detect_expression_format(m),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  expect_lt(elapsed, 5)                # must be sub-second-ish even on CI
  expect_length(warns, 0)               # no NAs-introduced-by-coercion noise
  expect_equal(d$format, "raw_counts")
  expect_equal(d$sc_or_bulk, "single_cell")
  exact <- 1 - length(m@x) / (as.numeric(nrow(m)) * as.numeric(ncol(m)))
  expect_lt(abs(d$sparsity - exact), 0.001)
})
