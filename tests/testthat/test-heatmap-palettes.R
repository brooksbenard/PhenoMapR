test_that("list_marker_heatmap_color_palettes returns brewer and viridis names", {
  pal <- list_marker_heatmap_color_palettes()
  expect_type(pal, "list")
  expect_named(pal, c("brewer", "viridis"))
  expect_type(pal$viridis, "character")
  expect_true(length(pal$viridis) >= 8L)
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    expect_true(length(pal$brewer$diverging) > 0L)
    expect_true(length(pal$brewer$qualitative) > 0L)
  }
})

test_that("default marker heatmap palettes match package conventions", {
  pg <- PhenoMapR:::.default_phenotype_palette()
  expect_equal(unname(pg[c("Most Favorable", "Other", "Most Adverse")]),
               c("#2166AC", "#F7F7F7", "#B2182B"))
  sc <- PhenoMapR:::.default_score_colors()
  expect_equal(sc, c("#2166AC", "#FFFFFF", "#B2182B"))
  ex <- PhenoMapR:::.default_expression_colors()
  expect_length(ex, 11L)
})

test_that("resolve phenotype and score palettes from brewer shorthand", {
  skip_if_not_installed("RColorBrewer")
  pg <- PhenoMapR:::.resolve_phenotype_palette("brewer:RdBu")
  expect_named(pg, c("Most Favorable", "Other", "Most Adverse"))
  expect_length(pg, 3L)
  sc <- PhenoMapR:::.resolve_score_colors("brewer:RdBu")
  expect_length(sc, 3L)
})

test_that("manual color_schemes work in plot_phenotype_markers", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  genes <- paste0("G", 1:10)
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
    celltype_original = rep(c("T1", "T2"), 6),
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

  ht <- plot_phenotype_markers(
    markers = markers,
    expr_mat = expr,
    meta = meta,
    group_col = "phenotype_group",
    score_col = "score",
    heatmap_type = "global",
    top_n_markers = 2L,
    n_mark_labels = 1L,
    draw = FALSE,
    color_schemes = list(
      phenotype = list(
        source = "manual",
        colors = c("#111111", "#222222", "#333333")
      ),
      score = list(
        source = "manual",
        colors = c("#444444", "#FFFFFF", "#555555")
      ),
      expression = list(
        source = "manual",
        colors = c("#000000", "#FFFFFF", "#FF0000")
      )
    )
  )
  expect_s4_class(ht, "Heatmap")
})
