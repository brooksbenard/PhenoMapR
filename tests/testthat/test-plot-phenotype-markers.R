# plot_phenotype_markers — requires ComplexHeatmap + circlize (Suggests)

test_that(".global_marker_heatmap_cell_order sorts columns by score low to high", {
  meta <- data.frame(
    id = paste0("C", 1:5),
    s = c(3, 1, 2, NA, 0),
    stringsAsFactors = FALSE
  )
  expr <- matrix(1, nrow = 2, ncol = 5, dimnames = list(c("G1", "G2"), paste0("C", 1:5)))
  o <- PhenoMapR:::.global_marker_heatmap_cell_order(
    meta = meta,
    expr_mat = expr,
    cell_id_col = "id",
    score_col = "s"
  )
  expect_equal(o$cell_order, c("C5", "C2", "C3", "C1", "C4"))
})

test_that("plot_phenotype_markers returns Heatmap object with fake marker tables", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(3)
  genes <- paste0("G", 1:20)
  cells <- paste0("C", 1:30)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 10),
    score = rnorm(30),
    celltype_original = rep(c("T1", "T2", "T3"), 10),
    stringsAsFactors = FALSE
  )

  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2", "G3"),
      avg_log2FC = c(1.2, 1.0, 0.8),
      p_adj = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G4", "G5", "G6"),
      avg_log2FC = c(1.1, 0.9, 0.7),
      p_adj = c(0.01, 0.02, 0.04),
      stringsAsFactors = FALSE
    )
  )

  ht <- plot_phenotype_markers(
    markers = markers,
    expr_mat = expr,
    meta = meta,
    group_col = "phenotype_group",
    score_col = "score",
    heatmap_type = "global",
    top_n_markers = 5L,
    n_mark_labels = 2L,
    p_adj_threshold = 0.05,
    heatmap_width = grid::unit(90, "mm"),
    heatmap_height = grid::unit(55, "mm"),
    draw = FALSE
  )
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phenotype_markers cell_type_specific returns Heatmap with fake tables", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  genes <- paste0("G", 1:15)
  cells <- paste0("C", 1:24)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 8),
    score = rnorm(24),
    celltype_original = rep(c("T1", "T2"), 12),
    stringsAsFactors = FALSE
  )

  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2"),
      cell_type = c("T1", "T1"),
      avg_log2FC = c(1.2, 1.0),
      p_adj = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4"),
      cell_type = c("T2", "T2"),
      avg_log2FC = c(1.1, 0.9),
      p_adj = c(0.01, 0.03),
      stringsAsFactors = FALSE
    )
  )

  ht <- plot_phenotype_markers(
    markers = markers,
    expr_mat = expr,
    meta = meta,
    group_col = "phenotype_group",
    score_col = "score",
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L,
    n_mark_labels = 2L,
    p_adj_threshold = 0.05,
    heatmap_width = 120,  # mm
    draw = FALSE
  )
  expect_s4_class(ht, "Heatmap")
})

test_that(".pick_marker_genes can rank by lfc vs p_adj with different filters", {
  df <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    avg_log2FC = c(2.0, 1.5, 0.8, 1.2),
    p_adj = c(0.04, 0.001, 1e-6, 0.03),
    stringsAsFactors = FALSE
  )

  # LFC ranking: significant (p_adj < 0.05) and avg_log2FC > 0; pick highest LFC.
  g_lfc <- PhenoMapR:::.pick_marker_genes(df, n_keep = 3, p_adj_threshold = 0.05, rank_by = "lfc")
  expect_equal(g_lfc[1], "G1")

  # p_adj ranking: still significant, but enforce avg_log2FC > 1; then rank by p_adj.
  g_p <- PhenoMapR:::.pick_marker_genes(df, n_keep = 3, p_adj_threshold = 0.05, rank_by = "p_adj")
  expect_true(all(g_p %in% c("G1", "G2", "G4")))
  expect_false("G3" %in% g_p) # fails avg_log2FC > 1 filter
})

test_that(".meta_score_numeric_vector uses first column of matrix score column", {
  meta <- data.frame(
    id = 1:3,
    score = I(cbind(c(10, 20, 30))),
    stringsAsFactors = FALSE
  )
  v <- PhenoMapR:::.meta_score_numeric_vector(meta, "score")
  expect_equal(v, c(10, 20, 30))
})

test_that(".score_ann_for_heatmap indexes meta rows safely", {
  meta <- data.frame(score = c(1, 99, 3), stringsAsFactors = FALSE)
  idx <- c(1L, 3L, NA_integer_)
  expect_equal(PhenoMapR:::.score_ann_for_heatmap(meta, "score", idx), c(1, 3, NA_real_))
})

test_that(".rect_native_ct_marker_blocks returns jmin/jmax per row_split block", {
  row_split <- factor(
    c("Most Favorable||T1", "Most Favorable||T1", "Most Adverse||T2"),
    levels = c("Most Favorable||T1", "Most Adverse||T2")
  )
  col_block_key <- c(
    rep("Most Favorable||T1", 3L),
    rep("Most Adverse||T2", 2L)
  )
  r <- PhenoMapR:::.rect_native_ct_marker_blocks(col_block_key, row_split)
  expect_equal(nrow(r), 2L)
  expect_equal(r$block, c("Most Favorable||T1", "Most Adverse||T2"))
  expect_equal(r$jmin, c(1L, 4L))
  expect_equal(r$jmax, c(3L, 5L))
  expect_equal(r$r1, c(1L, 3L))
  expect_equal(r$r2, c(2L, 3L))
})
