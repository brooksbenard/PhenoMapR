test_that("collapse_probes_to_genes averages probes mapping to the same symbol", {
  expr <- matrix(c(2, 4, 6, 8), nrow = 2L, ncol = 2L, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("S1", "S2")))
  annot <- data.frame(
    ID = c("p1", "p2"),
    `Gene Symbol` = c("AAA", "AAA"),
    check.names = FALSE
  )
  out <- collapse_probes_to_genes(expr, annot)
  expect_equal(nrow(out), 1L)
  expect_equal(rownames(out), "AAA")
  expect_equal(as.numeric(out["AAA", "S1"]), 4)
  expect_equal(as.numeric(out["AAA", "S2"]), 6)
})

test_that("preprocess_microarray_expression quantile-normalizes and scales genes", {
  set.seed(1)
  x <- matrix(rnorm(40 * 6), nrow = 40, ncol = 6,
              dimnames = list(paste0("G", 1:40), paste0("S", 1:6)))
  out <- preprocess_microarray_expression(x, verbose = FALSE)
  expect_equal(dim(out), dim(x))
  gene_sd <- apply(out, 1, sd, na.rm = TRUE)
  expect_true(all(abs(gene_sd - 1) < 0.05 | gene_sd == 0))
})

test_that("derive_reference_from_bulk microarray path uses probe annotation", {
  skip_if_not_installed("HGNChelper")
  set.seed(2)
  bulk <- matrix(runif(6 * 12, min = 50, max = 5000), 6, 12,
                 dimnames = list(paste0("ILMN_", 1:6), paste0("S", 1:12)))
  pheno <- data.frame(
    sample_id = colnames(bulk),
    status = rep(c("A", "B"), each = 6),
    stringsAsFactors = FALSE
  )
  annot <- data.frame(
    ID = rownames(bulk),
    Symbol = c("TP53", "MYC", "EGFR", "TP53", "KRAS", "PIK3CA"),
    stringsAsFactors = FALSE
  )
  ref <- suppressWarnings(derive_reference_from_bulk(
    bulk_expression = bulk,
    phenotype = pheno,
    sample_id_column = "sample_id",
    phenotype_column = "status",
    phenotype_type = "binary",
    platform = "microarray",
    probe_annotation = annot,
    verbose = FALSE
  ))
  expect_s3_class(ref, "data.frame")
  expect_gt(nrow(ref), 3L)
  expect_equal(attr(ref, "platform"), "microarray")
})

test_that("derive_meta_z_from_bulk_studies combines per-study z-scores", {
  set.seed(3)
  make_study <- function(prefix, n = 20L) {
    genes <- paste0("G", 1:30)
    bulk <- matrix(rnorm(30 * n), 30, n,
                   dimnames = list(genes, paste0(prefix, "_S", seq_len(n))))
    pheno <- data.frame(
      sample_id = colnames(bulk),
      y = rnorm(n),
      stringsAsFactors = FALSE
    )
    list(
      bulk_expression = bulk,
      phenotype = pheno,
      sample_id_column = "sample_id",
      phenotype_column = "y",
      phenotype_type = "continuous"
    )
  }
  studies <- list(a = make_study("A"), b = make_study("B"))
  ref <- suppressWarnings(derive_meta_z_from_bulk_studies(
    studies, meta_z_label = "test_meta_z", verbose = FALSE
  ))
  expect_s3_class(ref, "data.frame")
  expect_equal(colnames(ref), "test_meta_z")
  expect_gt(nrow(ref), 10L)
  expect_equal(attr(ref, "n_studies"), 2L)
})
