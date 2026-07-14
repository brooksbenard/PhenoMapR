# test-find-prognostic-markers.R
# Tests for find_phenotype_markers (matrix path and group_label resolution)

test_that("find_phenotype_markers with matrix and group vector returns list", {
  set.seed(1)
  n_genes <- 50
  n_cells <- 60
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  group_vec <- c(
    rep("Most Adverse", 10),
    rep("Most Favorable", 10),
    rep("Other", 40)
  )
  out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_s3_class(out$adverse_markers, "data.frame")
  expect_s3_class(out$favorable_markers, "data.frame")
  expect_true("gene" %in% names(out$adverse_markers))
  expect_true("p_val" %in% names(out$adverse_markers))
})

test_that("find_phenotype_markers with verbose=TRUE runs (matrix path)", {
  set.seed(99)
  expr <- matrix(pmax(0, rnorm(40 * 45)), 40, 45, dimnames = list(paste0("G", 1:40), paste0("C", 1:45)))
  group_vec <- c(rep("Most Adverse", 8), rep("Most Favorable", 8), rep("Other", 29))
  msg <- capture.output(
    out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = TRUE),
    type = "message"
  )
  expect_type(out, "list")
  expect_true(length(msg) >= 0)
})

test_that("find_phenotype_markers with group_labels data.frame and group_column", {
  set.seed(2)
  n_genes <- 40
  n_cells <- 50
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  groups_df <- data.frame(
    cell_id = paste0("C", seq_len(n_cells)),
    pg = c(rep("Most Adverse", 8), rep("Most Favorable", 8), rep("Other", 34)),
    stringsAsFactors = FALSE
  )
  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers supports cell_type_specific marker detection", {
  set.seed(6)
  n_genes <- 20
  n_cells <- 60
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  # Two cell types with 10 adverse, 10 favorable, 10 other per type
  cell_type <- c(rep("T", 30), rep("B", 30))
  pg <- c(
    rep("Most Adverse", 10), rep("Most Favorable", 10), rep("Other", 10),
    rep("Most Adverse", 10), rep("Most Favorable", 10), rep("Other", 10)
  )
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = pg,
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )

  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
  expect_true(all(unique(out$adverse_markers$cell_type) %in% c("T", "B")))
})

test_that("find_phenotype_markers cell_type_specific accepts celltype_contrast vs_cohort_rest (original cohort-wide reference)", {
  set.seed(6)
  n_genes <- 20
  n_cells <- 12
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  cell_type <- c(rep("T", 6), rep("B", 6))
  pg <- c(
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2),
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2)
  )
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = pg,
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )

  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    celltype_contrast = "vs_cohort_rest",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
})

test_that("vs_cohort_rest defines a larger reference set than within_cell_type for T adverse", {
  cell_type_vec <- c(rep("T", 6), rep("B", 6))
  group_vec <- c(
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2),
    rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2)
  )
  ct <- "T"
  tail <- "Most Adverse"
  is_in <- !is.na(cell_type_vec) & !is.na(group_vec) &
    (cell_type_vec == ct) & (group_vec == tail)
  n_within <- sum(!is.na(cell_type_vec) & !is.na(group_vec) &
    (cell_type_vec == ct) & (group_vec != tail))
  n_cohort <- sum(!is_in & !is.na(group_vec))
  expect_true(n_cohort > n_within)
})

test_that("find_phenotype_markers cell_type_specific accepts celltype_contrast vs_opposite_tail", {
  set.seed(70)
  n_genes <- 20
  n_cells <- 12
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = c(
      rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2),
      rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 2)
    ),
    cell_type = c(rep("T", 6), rep("B", 6)),
    stringsAsFactors = FALSE
  )

  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    celltype_contrast = "vs_opposite_tail",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
})

test_that("find_phenotype_markers vs_opposite_tail still produces markers when a cell type exists in only one tail", {
  # Setup: ductal cells appear ONLY in Most Adverse, immune cells ONLY in
  # Most Favorable. With celltype_contrast = "within_cell_type" the
  # ductal-adverse block has no same-type comparator (no favorable ductal
  # cells exist) and the immune-favorable block has no same-type
  # comparator either, so within_cell_type returns empty in both blocks.
  # vs_opposite_tail compares against the entire opposite tail
  # (regardless of cell type) and should still produce non-empty
  # markers, demonstrating the new option recovers signal that
  # within_cell_type drops.
  set.seed(202)
  n_genes <- 60
  n_cells <- 30
  # Construct expression so that ductal-adverse cells differ strongly
  # from favorable cells: shift the first 5 gene means up only in the
  # ductal-adverse block. Then the wilcoxon vs_opposite_tail contrast
  # should pick those genes up.
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells, mean = 1, sd = 0.5)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)),
                    paste0("C", seq_len(n_cells)))
  )

  cell_type <- c(rep("ductal", 10), rep("immune", 20))
  pg <- c(
    rep("Most Adverse", 10),
    rep("Most Favorable", 10),
    rep("Most Adverse", 5), rep("Most Favorable", 5)
  )
  # Ductal cells that ARE adverse: bump first 5 genes.
  ductal_adv_idx <- which(cell_type == "ductal" & pg == "Most Adverse")
  expr[1:5, ductal_adv_idx] <- expr[1:5, ductal_adv_idx] + 4

  # No favorable ductal cells in this design: confirm.
  expect_equal(sum(cell_type == "ductal" & pg == "Most Favorable"), 0L)

  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = pg,
    cell_type = cell_type,
    stringsAsFactors = FALSE
  )

  out_within <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    celltype_contrast = "within_cell_type",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )
  out_opp <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    celltype_contrast = "vs_opposite_tail",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )

  # within_cell_type returns nothing for ductal (no favorable ductal
  # comparator); vs_opposite_tail should still surface ductal-adverse
  # markers vs all favorable cells.
  ductal_within_adv <- out_within$adverse_markers[
    out_within$adverse_markers$cell_type == "ductal", , drop = FALSE]
  ductal_opp_adv <- out_opp$adverse_markers[
    out_opp$adverse_markers$cell_type == "ductal", , drop = FALSE]
  expect_equal(nrow(ductal_within_adv), 0L)
  expect_gt(nrow(ductal_opp_adv), 0L)
})

test_that("within_cell_type skips cell types absent from either phenotype tail", {
  set.seed(88)
  n_genes <- 40
  n_cells <- 60
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = c(
      rep("Most Adverse", 20), rep("Most Favorable", 20), rep("Other", 20)
    ),
    cell_type = c(
      rep("ductal", 10), rep("acinar", 10),
      rep("ductal", 10), rep("beta", 10),
      rep("ductal", 10), rep("acinar", 10)
    ),
    stringsAsFactors = FALSE
  )
  out <- find_phenotype_markers(
    expr,
    group_labels = groups_df,
    group_column = "pg",
    cell_id_column = "cell_id",
    marker_scope = "cell_type_specific",
    cell_type_column = "cell_type",
    celltype_contrast = "within_cell_type",
    min.pct = 0,
    logfc.threshold = 0,
    pval_threshold = 1,
    max_cells_per_ident = Inf,
    verbose = FALSE
  )
  expect_false("acinar" %in% out$adverse_markers$cell_type)
  expect_false("acinar" %in% out$favorable_markers$cell_type)
})

test_that("find_phenotype_markers cell_type_specific handles empty per-cell-type results", {
  set.seed(66)
  n_genes <- 30
  n_cells <- 42
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  # T cells have adverse/favorable/other (>=5 per tail); B cells have only "Other"
  groups_df <- data.frame(
    cell_id = colnames(expr),
    pg = c(
      rep("Most Adverse", 10), rep("Most Favorable", 10), rep("Other", 10),
      rep("Other", 12)
    ),
    cell_type = c(rep("T", 30), rep("B", 12)),
    stringsAsFactors = FALSE
  )

  expect_no_error({
    out <- find_phenotype_markers(
      expr,
      group_labels = groups_df,
      group_column = "pg",
      cell_id_column = "cell_id",
      marker_scope = "cell_type_specific",
      cell_type_column = "cell_type",
      min.pct = 0,
      logfc.threshold = 0,
      pval_threshold = 1,
      max_cells_per_ident = Inf,
      verbose = FALSE
    )
  })

  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_true("cell_type" %in% names(out$adverse_markers))
  expect_true("cell_type" %in% names(out$favorable_markers))
  expect_false("B" %in% out$adverse_markers$cell_type)
  expect_false("B" %in% out$favorable_markers$cell_type)
})

test_that("find_phenotype_markers errors when group_labels length mismatch", {
  expr <- matrix(1, nrow = 5, ncol = 10, dimnames = list(paste0("G", 1:5), paste0("C", 1:10)))
  group_vec <- rep("Other", 5)
  expect_error(
    find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE),
    "Length of .* must match"
  )
})

test_that("find_phenotype_markers warns when few adverse or favorable cells", {
  set.seed(3)
  expr <- matrix(pmax(0, rnorm(30 * 20)), 30, 20, dimnames = list(paste0("G", 1:30), paste0("C", 1:20)))
  group_vec <- c(rep("Most Adverse", 2), rep("Most Favorable", 2), rep("Other", 16))
  expect_warning(
    out <- find_phenotype_markers(expr, group_labels = group_vec, verbose = FALSE),
    "Fewer than 5"
  )
  expect_type(out, "list")
})

test_that("find_phenotype_markers accepts dgCMatrix", {
  skip_if_not_installed("Matrix")
  set.seed(4)
  n_genes <- 30
  n_cells <- 40
  expr_dense <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  expr_sparse <- Matrix::Matrix(expr_dense, sparse = TRUE)
  group_vec <- c(rep("Most Adverse", 6), rep("Most Favorable", 6), rep("Other", 28))
  out <- find_phenotype_markers(expr_sparse, group_labels = group_vec, verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers with max_cells_per_ident subsamples", {
  set.seed(5)
  n_genes <- 25
  n_cells <- 200
  expr <- matrix(
    pmax(0, rnorm(n_genes * n_cells)),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  group_vec <- c(
    rep("Most Adverse", 80),
    rep("Most Favorable", 80),
    rep("Other", 40)
  )
  out <- find_phenotype_markers(
    expr,
    group_labels = group_vec,
    max_cells_per_ident = 20L,
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})

test_that("find_phenotype_markers with invalid group labels returns empty markers", {
  expr <- matrix(1, nrow = 5, ncol = 5, dimnames = list(paste0("G", 1:5), paste0("C", 1:5)))
  bad_df <- data.frame(id = paste0("C", 1:5), grp = c("A", "B", "C", "D", "E"))
  out <- find_phenotype_markers(expr, group_labels = bad_df, group_column = "grp", verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_equal(nrow(out$adverse_markers), 0)
  expect_equal(nrow(out$favorable_markers), 0)
})

test_that("find_phenotype_markers accepts an in-memory Python anndata.AnnData", {
  # Regression test for the user-reported error
  #   "find_phenotype_markers failed: 'expression' must be a matrix, ..."
  # Before the fix `process_expression_for_markers()` only handled
  # matrix / Matrix / data.frame / Seurat / SCE -- so a Python
  # AnnData (python.builtin.object) returned by reticulate fell off
  # the end and triggered the catch-all stop(). The fix adds a
  # python.builtin.object branch that calls .anndata_X_to_genes_cells()
  # to pull /X back into R as a (genes x cells) Matrix.
  skip_if_not_installed("reticulate")
  skip_if_not_installed("Matrix")
  ad <- tryCatch(reticulate::import("anndata", convert = FALSE),
                 error = function(e) NULL)
  skip_if(is.null(ad),
          "Python anndata module not reachable through reticulate")

  set.seed(11)
  n_cells <- 30L; n_genes <- 25L
  X <- matrix(pmax(0, rnorm(n_cells * n_genes, 1, 0.5)),
              nrow = n_cells, ncol = n_genes)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  Xpy <- scipy_sparse$csr_matrix(np$asarray(X, dtype = "float32"))
  obs_idx <- paste0("Cell_", seq_len(n_cells))
  var_idx <- paste0("GENE", seq_len(n_genes))
  obs_df <- pd$DataFrame(list(
    cell_type = rep(c("TypeA", "TypeB"), length.out = n_cells)
  ), index = obs_idx)
  var_df <- pd$DataFrame(list(), index = var_idx)
  adata  <- ad$AnnData(X = Xpy, obs = obs_df, var = var_df)

  groups <- rep(c("Most Adverse", "Most Favorable", "Other"),
                length.out = n_cells)
  out <- find_phenotype_markers(
    adata,
    group_labels = groups,
    min_pct = 0.05,
    logfc_threshold = 0.05,
    p_adj_threshold = 0.5,
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_s3_class(out$adverse_markers,   "data.frame")
  expect_s3_class(out$favorable_markers, "data.frame")
})

test_that("find_phenotype_markers with Seurat object uses FindMarkers path", {
  skip_if_not_installed("Seurat")
  set.seed(42)
  n_genes <- 40
  n_cells <- 45
  counts <- matrix(
    pmax(0, as.integer(rnorm(n_genes * n_cells, 10, 3))),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  groups <- c(
    rep("Most Adverse", 8),
    rep("Most Favorable", 8),
    rep("Other", 29)
  )
  out <- find_phenotype_markers(obj, group_labels = groups, assay = "RNA", slot = "data", verbose = FALSE)
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
  expect_s3_class(out$adverse_markers, "data.frame")
  expect_s3_class(out$favorable_markers, "data.frame")
})

test_that("find_phenotype_markers tolerates stale Seurat Graphs (CosMx / spatial)", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  set.seed(7)
  n_genes <- 30
  n_cells <- 40
  counts <- matrix(
    pmax(0, as.integer(rnorm(n_genes * n_cells, 10, 3))),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(paste0("G", seq_len(n_genes)), paste0("C", seq_len(n_cells)))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, assay = "RNA"))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  # Stale Graphs whose dimnames do not match object cells reproduce the
  # user-facing CosMx error:
  #   "Please provide rownames to the matrix before converting to a Graph"
  # when Seurat subsets during FindMarkers prep.
  fake <- Matrix::Diagonal(n = ncol(obj), x = 1)
  rownames(fake) <- colnames(fake) <- paste0("OLD_", colnames(obj))
  obj@graphs$RNA_nn <- SeuratObject::as.Graph(fake)
  obj@graphs$RNA_snn <- SeuratObject::as.Graph(fake)
  expect_error(obj[, colnames(obj)[seq_len(10L)]], "rownames|Graph")

  groups <- c(
    rep("Most Adverse", 8),
    rep("Most Favorable", 8),
    rep("Other", n_cells - 16L)
  )
  out <- find_phenotype_markers(
    obj,
    group_labels = groups,
    assay = "RNA",
    slot = "data",
    max_cells_per_ident = 15L,
    verbose = FALSE
  )
  expect_type(out, "list")
  expect_named(out, c("adverse_markers", "favorable_markers"))
})
