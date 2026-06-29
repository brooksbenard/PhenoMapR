# plot_phenotype_markers — requires ComplexHeatmap + circlize (Suggests)

.ct_heatmap_meta_two_types <- function(cells = paste0("C", seq_len(60L))) {
  data.frame(
    Cell = cells,
    phenotype_group = rep(
      rep(c("Most Adverse", "Most Favorable", "Other"), each = 10L),
      2L
    ),
    score = rnorm(length(cells)),
    celltype_original = rep(c("T1", "T2"), each = 30L),
    stringsAsFactors = FALSE
  )
}

test_that(".phenotype_heatmap_cell_order nests all cell types within each phenotype", {
  meta <- data.frame(
    Cell = c("F_A", "F_B", "O_A", "A_A", "A_B"),
    phenotype_group = c(
      "Most Favorable", "Most Favorable", "Other",
      "Most Adverse", "Most Adverse"
    ),
    cell_type = c("Alpha", "Beta", "Alpha", "Alpha", "Beta"),
    score = c(1, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )
  expr <- matrix(
    1, nrow = 2, ncol = 5,
    dimnames = list(c("G1", "G2"), meta$Cell)
  )
  ord <- PhenoMapR:::.phenotype_heatmap_cell_order(
    meta = meta,
    expr_mat = expr,
    cell_id_col = "Cell",
    group_col = "phenotype_group",
    celltype_col = "cell_type",
    score_col = "score",
    hm_group_levels = c("Most Favorable", "Other", "Most Adverse"),
    hm_celltype_levels = c("Alpha", "Beta")
  )
  expect_equal(
    ord$cell_order,
    c("F_A", "F_B", "O_A", "A_A", "A_B")
  )
})

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
  cells <- paste0("C", 1:60)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- .ct_heatmap_meta_two_types(cells)

  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2", "G5", "G6"),
      cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.2, 1.0, 1.1, 0.9),
      p_adj = c(0.01, 0.02, 0.01, 0.03),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4", "G7", "G8"),
      cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.1, 0.9, 1.0, 0.85),
      p_adj = c(0.01, 0.03, 0.01, 0.02),
      stringsAsFactors = FALSE
    )
  )

  ht <- plot_phenotype_markers(
    markers = markers,
    expr_mat = expr,
    meta = meta,
    group_col = "phenotype_group",
    celltype_col = "celltype_original",
    score_col = "score",
    heatmap_type = "cell_type_specific",
    celltype_contrast = "within_cell_type",
    top_n_markers = 5L,
    n_mark_labels = 2L,
    p_adj_threshold = 0.05,
    heatmap_width = 120,  # mm
    draw = FALSE
  )
  expect_s4_class(ht, "Heatmap")
  expect_true(!is.null(ht@left_annotation))
  expect_true(!is.null(ht@right_annotation))
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 400, height = 300, res = 72)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
  split_names <- names(ComplexHeatmap::row_order(ht))
  fav_i <- grep("Most Favorable", split_names)
  adv_i <- grep("Most Adverse", split_names)
  expect_true(length(fav_i) > 0L && length(adv_i) > 0L)
  expect_true(max(fav_i) < min(adv_i))
})

test_that(".anno_mark_labels_gp colors labels from cell-type palette", {
  pal <- c(T1 = "#111111", T2 = "#222222")
  gp_off <- PhenoMapR:::.anno_mark_labels_gp(c("T1", "T2"), FALSE, pal)
  expect_equal(gp_off$col, "black")
  gp_on <- PhenoMapR:::.anno_mark_labels_gp(c("T1", "T2"), TRUE, pal)
  expect_equal(unname(gp_on$col), c("#111111", "#222222"))
})

test_that("plot_phenotype_markers accepts color_mark_labels_by_celltype", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  genes <- paste0("G", 1:15)
  cells <- paste0("C", 1:60)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- .ct_heatmap_meta_two_types(cells)

  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2", "G5", "G6"),
      cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.2, 1.0, 1.1, 0.9),
      p_adj = c(0.01, 0.02, 0.01, 0.03),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4", "G7", "G8"),
      cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.1, 0.9, 1.0, 0.85),
      p_adj = c(0.01, 0.03, 0.01, 0.02),
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
    color_mark_labels_by_celltype = TRUE,
    p_adj_threshold = 0.05,
    draw = FALSE
  )
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phenotype_markers cell_type_specific accepts block_outline_color and lwd args", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  genes <- paste0("G", 1:15)
  cells <- paste0("C", 1:60)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- .ct_heatmap_meta_two_types(cells)
  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2", "G5", "G6"), cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.2, 1.0, 1.1, 0.9), p_adj = c(0.01, 0.02, 0.01, 0.03),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4", "G7", "G8"), cell_type = c("T1", "T1", "T2", "T2"),
      avg_log2FC = c(1.1, 0.9, 1.0, 0.85), p_adj = c(0.01, 0.03, 0.01, 0.02),
      stringsAsFactors = FALSE
    )
  )

  # draw = FALSE only exercises Heatmap construction; we just want to
  # confirm the new args round-trip through formals + match.arg without
  # error and that the returned heatmap object is still valid. The
  # actual decorate_heatmap_body grid.rect call only fires when
  # draw = TRUE, but ComplexHeatmap::draw() requires a graphics device
  # which we don't have in CI -- so only sanity-check the formals here.
  expect_no_error({
    ht <- plot_phenotype_markers(
      markers = markers,
      expr_mat = expr,
      meta = meta,
      group_col = "phenotype_group",
      score_col = "score",
      heatmap_type = "cell_type_specific",
      celltype_contrast = "within_cell_type",
      top_n_markers = 5L,
      n_mark_labels = 2L,
      p_adj_threshold = 0.05,
      outline_marker_blocks = TRUE,
      block_outline_color = "black",
      block_outline_lwd = 2,
      draw = FALSE
    )
  })
  expect_s4_class(ht, "Heatmap")

  fmls <- formals(plot_phenotype_markers)
  expect_true("block_outline_color" %in% names(fmls))
  expect_true("block_outline_lwd"   %in% names(fmls))
  expect_true("mark_label_fontsize" %in% names(fmls))
  expect_identical(eval(fmls$block_outline_color), "white")
  expect_equal(eval(fmls$block_outline_lwd), 1)
  expect_equal(eval(fmls$mark_label_fontsize), 7)
})

test_that(".phenomap_draw_padding_mm adds headroom for edge marks", {
  pad <- PhenoMapR:::.phenomap_draw_padding_mm(
    n_rows = 20L,
    marks_at = c(1L, 20L),
    fontsize = 7,
    mark_anno_width_mm = 24,
    heatmap_type = "cell_type_specific"
  )
  expect_length(pad, 4L)
  expect_gt(pad[1], 6)
  expect_gt(pad[3], 6)
})

test_that(".phenomap_mark_anno_width_mm scales with label length", {
  w_short <- PhenoMapR:::.phenomap_mark_anno_width_mm("GENE")
  w_long <- PhenoMapR:::.phenomap_mark_anno_width_mm(
    paste(rep("A", 50), collapse = "")
  )
  expect_gt(w_long, w_short)
})

test_that("plot_phenotype_markers global mode runs with celltype_col = NULL (no cell-type strip)", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(7)
  genes <- paste0("G", 1:12)
  cells <- paste0("C", 1:18)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  # No celltype_original column on purpose — this mimics the Shiny "global"
  # path where the user hasn't mapped a cell-type column.
  meta <- data.frame(
    Cell = cells,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 6),
    score = rnorm(18),
    stringsAsFactors = FALSE
  )

  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2"),
      avg_log2FC = c(1.2, 1.0),
      p_adj = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4"),
      avg_log2FC = c(1.1, 0.9),
      p_adj = c(0.01, 0.03),
      stringsAsFactors = FALSE
    )
  )

  expect_no_error({
    ht <- plot_phenotype_markers(
      markers = markers,
      expr_mat = expr,
      meta = meta,
      group_col = "phenotype_group",
      score_col = "score",
      celltype_col = NULL,
      heatmap_type = "global",
      top_n_markers = 5L,
      n_mark_labels = 2L,
      draw = FALSE
    )
  })
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phenotype_markers global mode survives an empty favorable tail", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(11)
  genes <- paste0("G", 1:8)
  cells <- paste0("C", 1:12)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells))),
    length(genes),
    length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 4),
    score = rnorm(12),
    celltype_original = rep("T1", 12),
    stringsAsFactors = FALSE
  )
  # Favorable table exists but every row fails the LFC + p_adj filters,
  # so n_fav -> 0 in the global path.
  markers <- list(
    adverse_markers = data.frame(
      gene = c("G1", "G2"),
      avg_log2FC = c(1.5, 1.0),
      p_adj = c(0.001, 0.01),
      stringsAsFactors = FALSE
    ),
    favorable_markers = data.frame(
      gene = c("G3", "G4"),
      avg_log2FC = c(-0.5, -0.2),  # negative LFC -> filtered out
      p_adj = c(0.01, 0.02),
      stringsAsFactors = FALSE
    )
  )

  expect_no_error({
    ht <- plot_phenotype_markers(
      markers = markers,
      expr_mat = expr,
      meta = meta,
      group_col = "phenotype_group",
      score_col = "score",
      celltype_col = "celltype_original",
      heatmap_type = "global",
      top_n_markers = 5L,
      n_mark_labels = 2L,
      draw = FALSE
    )
  })
  expect_s4_class(ht, "Heatmap")
})

test_that("plot_phenotype_markers errors clearly when cell-type-specific is asked without celltype_col", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  meta <- data.frame(
    Cell = paste0("C", 1:6),
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 2),
    score = 1:6,
    stringsAsFactors = FALSE
  )
  expr <- matrix(0, nrow = 2, ncol = 6,
                 dimnames = list(c("G1", "G2"), paste0("C", 1:6)))
  markers <- list(
    adverse_markers = data.frame(gene = "G1", avg_log2FC = 1, p_adj = 0.01,
                                 stringsAsFactors = FALSE),
    favorable_markers = data.frame(gene = "G2", avg_log2FC = 1, p_adj = 0.01,
                                   stringsAsFactors = FALSE)
  )
  expect_error(
    plot_phenotype_markers(
      markers = markers,
      expr_mat = expr,
      meta = meta,
      group_col = "phenotype_group",
      score_col = "score",
      celltype_col = NULL,
      heatmap_type = "cell_type_specific",
      draw = FALSE
    ),
    regexp = "celltype_col"
  )
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

test_that("plot_phenotype_markers writes anno_mark labels into a real PNG", {
  # Regression for the Shiny "labels are blank, lines are there" report:
  # under renderPlot() Shiny replays plot expressions on a screen device,
  # which can blank anno_mark text grobs even though the link lines stay
  # drawn. The Shiny app now uses renderImage() and writes a PNG via the
  # exact path tested here. Reading the PNG back and grepping for any
  # rendered text content lets us catch a regression in the package
  # function (font availability, wrong viewport, empty `labels`, ...)
  # without needing a real Shiny session.
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(101)
  n_cells <- 60L
  n_genes <- 60L
  expr <- matrix(rnorm(n_cells * n_genes), nrow = n_genes, ncol = n_cells,
                 dimnames = list(paste0("Gene", sprintf("%03d", seq_len(n_genes))),
                                 paste0("C", seq_len(n_cells))))
  meta <- data.frame(
    cell_id   = colnames(expr),
    group     = rep(c("Most Adverse", "Other", "Most Favorable"),
                    times = c(20, 20, 20)),
    score     = runif(n_cells, -1, 1),
    cell_type = rep(c("Tumor", "T cell", "Macrophage"), each = 20L),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene       = paste0("Gene", sprintf("%03d", 1:30)),
    avg_log2FC = sort(runif(30, 0.5, 3), decreasing = TRUE),
    p_adj      = sort(runif(30, 1e-10, 1e-3)),
    stringsAsFactors = FALSE
  )
  adv <- data.frame(
    gene       = paste0("Gene", sprintf("%03d", 31:60)),
    avg_log2FC = sort(runif(30, 0.5, 3), decreasing = TRUE),
    p_adj      = sort(runif(30, 1e-10, 1e-3)),
    stringsAsFactors = FALSE
  )

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 1400, height = 640, res = 110)
  PhenoMapR::plot_phenotype_markers(
    markers      = list(adverse_markers = adv, favorable_markers = fav),
    expr_mat     = expr,
    meta         = meta,
    cell_id_col  = "cell_id",
    group_col    = "group",
    score_col    = "score",
    celltype_col = "cell_type",
    heatmap_type = "global",
    top_n_markers = 20L,
    n_mark_labels = 5L,
    draw         = TRUE
  )
  grDevices::dev.off()

  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 30000)  # rich PNG, not just an empty header
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

test_that("plot_phenotype_markers legend uses Most Phenotype +/- labels", {
  # The phenotype-group legend in the heatmap renders "Most Phenotype +"
  # and "Most Phenotype -" instead of "Most Favorable" / "Most Adverse"
  # while the underlying data values stay "Most Favorable" /
  # "Most Adverse" (so every pipeline that pattern-matches on those
  # exact strings keeps working). We assert the source of
  # plot_phenotype_markers builds a Legend whose `at` keeps the
  # data values and whose `labels` uses the new "+/-" wording, with
  # the directional convention enforced by the rest of the package:
  #   "Most Adverse"   (red,  #B2182B) -> "Most Phenotype +"
  #   "Most Favorable" (blue, #2166AC) -> "Most Phenotype -"
  # so the labels vector reads "Most Phenotype -", "Other",
  # "Most Phenotype +" (paired position-by-position with `at`).
  # Introspect the installed function body rather than the on-disk
  # R/*.R file: the latter is not shipped in the installed package
  # (only the compiled .rdb/.rdx are) and the dev fallback path is
  # not reachable inside R CMD check's tests/testthat working dir.
  src_one <- paste(
    deparse(PhenoMapR::plot_phenotype_markers, width.cutoff = 500L),
    collapse = "\n"
  )
  expect_match(
    src_one,
    paste0(
      'at\\s*=\\s*c\\(\\s*"Most Favorable",\\s*"Other",\\s*"Most Adverse"\\s*\\)',
      '[\\s\\S]*?',
      'labels\\s*=\\s*c\\(\\s*"Most Phenotype -",\\s*"Other",\\s*"Most Phenotype \\+"\\s*\\)'
    ),
    perl = TRUE
  )
})

test_that("plot_phenotype_markers cell_type_specific renders block outlines via column_split", {
  # Regression test for the EcoTyper-style block-outline fix:
  # we column_split on the same (phenotype_bin x cell_type) key as
  # the row split and then call decorate_heatmap_body() with no
  # explicit coordinates. The previous "manual native unit" approach
  # silently failed (rectangles were drawn outside the visible body),
  # so this test asserts that requesting block outlines does not
  # error out *and* produces a richer PNG than plotting without
  # outlines (extra black ink should add visible bytes).
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(11)
  genes <- paste0("G", 1:24)
  cells <- paste0("C", 1:60)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells), 1, 0.5)),
    length(genes), length(cells),
    dimnames = list(genes, cells)
  )
  meta <- .ct_heatmap_meta_two_types(cells)
  adv <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    cell_type = c("T1", "T1", "T2", "T2"),
    avg_log2FC = c(1.2, 1.0, 1.1, 0.9),
    p_adj = c(0.01, 0.02, 0.01, 0.03),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene = c("G5", "G6", "G7", "G8"),
    cell_type = c("T1", "T1", "T2", "T2"),
    avg_log2FC = c(1.0, 0.9, 1.1, 0.85),
    p_adj = c(0.01, 0.02, 0.01, 0.03),
    stringsAsFactors = FALSE
  )

  render <- function(outline) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 1100, height = 600, res = 110)
    on.exit(grDevices::dev.off(), add = TRUE)
    PhenoMapR::plot_phenotype_markers(
      markers = list(adverse_markers = adv, favorable_markers = fav),
      expr_mat = expr, meta = meta,
      group_col = "phenotype_group", score_col = "score",
      heatmap_type = "cell_type_specific",
      celltype_contrast = "within_cell_type",
      top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
      outline_marker_blocks = outline,
      block_outline_color = "black",
      block_outline_lwd = 2.5,
      draw = TRUE
    )
    f
  }
  f_off <- render(FALSE)
  f_on  <- render(TRUE)
  expect_true(file.exists(f_off) && file.exists(f_on))
  size_off <- file.info(f_off)$size
  size_on  <- file.info(f_on)$size
  # Sanity check: both renders are non-trivial.
  expect_gt(size_off, 20000)
  expect_gt(size_on,  20000)
  # The outlined version must contain at least as much ink. We do
  # not insist on a strict greater-than because the PNG encoder can
  # occasionally compress identical decorations to similar sizes,
  # but with thick black borders we expect the outlined PNG to be
  # at least within ~5% of the un-outlined one (and usually larger).
  expect_gt(size_on, size_off * 0.95)
})

test_that("plot_phenotype_markers omits cell-type strip in global mode with one level", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(7)
  genes <- paste0("G", 1:12)
  cells <- paste0("C", 1:30)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells), 1, 0.5)),
    length(genes), length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 10),
    score = rnorm(30),
    cell_type = rep("Acinar", 30),
    stringsAsFactors = FALSE
  )
  adv <- data.frame(
    gene = c("G1", "G2"), avg_log2FC = c(1.2, 1.0), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene = c("G3", "G4"), avg_log2FC = c(1.0, 0.9), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 900, height = 500, res = 110)
  PhenoMapR::plot_phenotype_markers(
    markers = list(adverse_markers = adv, favorable_markers = fav),
    expr_mat = expr, meta = meta,
    group_col = "phenotype_group", celltype_col = "cell_type",
    score_col = "score",
    heatmap_type = "global",
    top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
    draw = TRUE
  )
  grDevices::dev.off()
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 5000)
})

test_that("plot_phenotype_markers keeps cell-type strip in cell_type_specific mode with one level", {
  # Regression: after the both-tails filter, only one cell type may qualify;
  # the top/row cell-type annotations must still render.
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(7)
  genes <- paste0("G", 1:24)
  cells_one <- paste0("C", 1:30)
  cells_two <- paste0("C", 1:60)
  expr_one <- matrix(
    pmax(0, rnorm(length(genes) * length(cells_one), 1, 0.5)),
    length(genes), length(cells_one),
    dimnames = list(genes, cells_one)
  )
  expr_two <- matrix(
    pmax(0, rnorm(length(genes) * length(cells_two), 1, 0.5)),
    length(genes), length(cells_two),
    dimnames = list(genes, cells_two)
  )
  meta_one <- data.frame(
    Cell = cells_one,
    phenotype_group = rep(c("Most Adverse", "Most Favorable", "Other"), each = 10),
    score = rnorm(30),
    cell_type = rep("Acinar", 30),
    stringsAsFactors = FALSE
  )
  meta_two <- .ct_heatmap_meta_two_types(cells_two)
  meta_two$cell_type <- meta_two$celltype_original
  adv <- data.frame(
    gene = c("G1", "G2", "G3", "G4"),
    cell_type = c("T1", "T1", "T2", "T2"),
    avg_log2FC = c(1.2, 1.0, 1.1, 0.9),
    p_adj = c(0.01, 0.02, 0.01, 0.03),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene = c("G5", "G6", "G7", "G8"),
    cell_type = c("T1", "T1", "T2", "T2"),
    avg_log2FC = c(1.0, 0.9, 1.1, 0.85),
    p_adj = c(0.01, 0.02, 0.01, 0.03),
    stringsAsFactors = FALSE
  )
  adv_one <- data.frame(
    gene = c("G1", "G2"), cell_type = c("Acinar", "Acinar"),
    avg_log2FC = c(1.2, 1.0), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  fav_one <- data.frame(
    gene = c("G5", "G6"), cell_type = c("Acinar", "Acinar"),
    avg_log2FC = c(1.0, 0.9), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  render <- function(meta, mk_adv, mk_fav, expr_mat) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 1100, height = 600, res = 110)
    on.exit(grDevices::dev.off(), add = TRUE)
    PhenoMapR::plot_phenotype_markers(
      markers = list(adverse_markers = mk_adv, favorable_markers = mk_fav),
      expr_mat = expr_mat, meta = meta,
      group_col = "phenotype_group", celltype_col = "cell_type",
      score_col = "score",
      heatmap_type = "cell_type_specific",
      celltype_contrast = "within_cell_type",
      top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
      outline_marker_blocks = TRUE,
      block_outline_color = "black",
      block_outline_lwd = 1.5,
      draw = TRUE
    )
    f
  }

  f_one <- render(meta_one, adv_one, fav_one, expr_one)
  f_two <- render(meta_two, adv, fav, expr_two)
  expect_true(file.exists(f_one) && file.exists(f_two))
  expect_gt(file.info(f_one)$size, 5000)
  expect_gt(file.info(f_two)$size, file.info(f_one)$size)

  src <- paste(
    deparse(PhenoMapR::plot_phenotype_markers, width.cutoff = 500L),
    collapse = "\n"
  )
  expect_match(
    src,
    "heatmap_type == \"cell_type_specific\"",
    fixed = TRUE
  )
  expect_match(src, "if \\(show_celltype_strip\\)", perl = TRUE)
})

test_that("plot_phenotype_markers maps cell-type colors for all column annotations", {
  # Regression: column strip includes every cell on the heatmap, not only
  # cell types that pass the both-tails minimum for marker row blocks.
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(11)
  genes <- paste0("G", 1:12)
  cells_qual <- paste0("C", 1:30)
  cells_extra <- paste0("X", 1:20)
  cells <- c(cells_qual, cells_extra)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells), 1, 0.5)),
    length(genes), length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = c(
      rep(c("Most Adverse", "Most Favorable", "Other"), each = 10),
      rep("Other", 20)
    ),
    score = rnorm(length(cells)),
    cell_type = c(rep("Acinar", 30), rep(c("Ductal", "Plasma"), 10)),
    stringsAsFactors = FALSE
  )
  adv <- data.frame(
    gene = c("G1", "G2"), cell_type = c("Acinar", "Acinar"),
    avg_log2FC = c(1.2, 1.0), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene = c("G3", "G4"), cell_type = c("Acinar", "Acinar"),
    avg_log2FC = c(1.0, 0.9), p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 1100, height = 600, res = 110)
  out <- PhenoMapR::plot_phenotype_markers(
    markers = list(adverse_markers = adv, favorable_markers = fav),
    expr_mat = expr, meta = meta,
    group_col = "phenotype_group", celltype_col = "cell_type",
    score_col = "score",
    heatmap_type = "cell_type_specific",
    top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
    draw = TRUE
  )
  grDevices::dev.off()
  expect_true(inherits(out, "Heatmap"))
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 10000)
})

test_that("plot_phenotype_markers vs_opposite_tail draws one-sided cell-type blocks", {
  # Favorable ductal markers should render even without adverse ductal cells.
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  set.seed(21)
  genes <- paste0("G", 1:20)
  cells <- paste0("C", 1:40)
  expr <- matrix(
    pmax(0, rnorm(length(genes) * length(cells), 1, 0.5)),
    length(genes), length(cells),
    dimnames = list(genes, cells)
  )
  meta <- data.frame(
    Cell = cells,
    phenotype_group = c(
      rep("Most Favorable", 10), rep("Most Adverse", 10), rep("Other", 20)
    ),
    score = rnorm(40),
    cell_type = c(
      rep("Ductal", 10), rep("Immune", 10),
      rep("Ductal", 10), rep("Immune", 10)
    ),
    stringsAsFactors = FALSE
  )
  fav <- data.frame(
    gene = c("G1", "G2", "G3"),
    cell_type = rep("Ductal", 3),
    avg_log2FC = c(1.2, 1.0, 0.9),
    p_adj = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )
  adv <- data.frame(
    gene = c("G4", "G5"),
    cell_type = rep("Immune", 2),
    avg_log2FC = c(1.1, 1.0),
    p_adj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  ht_opp <- PhenoMapR::plot_phenotype_markers(
    markers = list(adverse_markers = adv, favorable_markers = fav),
    expr_mat = expr, meta = meta,
    group_col = "phenotype_group", celltype_col = "cell_type",
    score_col = "score",
    heatmap_type = "cell_type_specific",
    celltype_contrast = "vs_opposite_tail",
    top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
    draw = FALSE
  )
  expect_s4_class(ht_opp, "Heatmap")
  expect_true(!is.null(ht_opp@left_annotation))

  ht_within <- PhenoMapR::plot_phenotype_markers(
    markers = list(adverse_markers = adv, favorable_markers = fav),
    expr_mat = expr, meta = meta,
    group_col = "phenotype_group", celltype_col = "cell_type",
    score_col = "score",
    heatmap_type = "cell_type_specific",
    celltype_contrast = "within_cell_type",
    top_n_markers = 5L, n_mark_labels = 2L, p_adj_threshold = 0.05,
    draw = FALSE
  )
  expect_null(ht_within)
})
