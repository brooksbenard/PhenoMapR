# Shiny-side helper tests. The helpers live in inst/shiny/helpers.R and are
# sourced by app.R; here we just source them into a local environment and
# call them as plain functions.

source_shiny_helpers <- function() {
  # Prefer the in-source helpers.R when running tests from a dev tree
  # (so just-edited helpers are exercised without reinstalling), and
  # fall back to the installed copy when the source tree is not
  # available (e.g. R CMD check from a tarball).
  src_helpers <- tryCatch(
    testthat::test_path("..", "..", "inst", "shiny", "helpers.R"),
    error = function(e) NULL
  )
  helpers <- if (!is.null(src_helpers) && file.exists(src_helpers)) {
    src_helpers
  } else {
    system.file("shiny", "helpers.R", package = "PhenoMapR")
  }
  testthat::skip_if(!nzchar(helpers) || !file.exists(helpers),
                    "shiny/helpers.R not found")
  env <- new.env(parent = globalenv())
  sys.source(helpers, envir = env)
  env
}

test_that("detect_metadata_embeddings finds 1-indexed Seurat-style coordinate pairs", {
  e <- source_shiny_helpers()
  md <- data.frame(
    cell_id = paste0("C", 1:3),
    UMAP_1 = c(1.1, 2.2, 3.3),
    UMAP_2 = c(0.1, -0.2, 0.4),
    cell_type = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  pairs <- e$detect_metadata_embeddings(md)
  expect_equal(names(pairs), "UMAP")
  expect_equal(pairs$UMAP$dim1, "UMAP_1")
  expect_equal(pairs$UMAP$dim2, "UMAP_2")
})

test_that("detect_metadata_embeddings handles 0-indexed Scanpy-style coordinate pairs", {
  e <- source_shiny_helpers()
  md <- data.frame(
    cell_id = paste0("C", 1:3),
    X_umap_0 = c(1.1, 2.2, 3.3),
    X_umap_1 = c(0.1, -0.2, 0.4),
    stringsAsFactors = FALSE
  )
  pairs <- e$detect_metadata_embeddings(md)
  expect_equal(names(pairs), "X_umap")
  expect_equal(pairs$X_umap$dim1, "X_umap_0")
  expect_equal(pairs$X_umap$dim2, "X_umap_1")
})

test_that("detect_metadata_embeddings finds multiple distinct embeddings", {
  e <- source_shiny_helpers()
  md <- data.frame(
    cell_id = paste0("C", 1:3),
    UMAP_1 = 1:3, UMAP_2 = 4:6,
    tSNE_1 = 7:9, tSNE_2 = 10:12,
    cluster = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  pairs <- e$detect_metadata_embeddings(md)
  expect_setequal(names(pairs), c("UMAP", "tSNE"))
})

test_that("detect_metadata_embeddings returns empty list when no pair is present", {
  e <- source_shiny_helpers()
  md <- data.frame(
    cell_id = paste0("C", 1:3),
    cell_type = c("a", "b", "c"),
    score = c(1.2, 3.4, 5.6),
    stringsAsFactors = FALSE
  )
  expect_equal(e$detect_metadata_embeddings(md), list())
})

test_that("detect_metadata_embeddings ignores non-numeric coordinate-like columns", {
  e <- source_shiny_helpers()
  md <- data.frame(
    UMAP_1 = c("a", "b", "c"),  # not numeric
    UMAP_2 = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  expect_equal(e$detect_metadata_embeddings(md), list())
})

test_that("extract_embedding_from_metadata builds a clean embedding data.frame", {
  e <- source_shiny_helpers()
  md <- data.frame(
    cell_id = paste0("C", 1:3),
    UMAP_1 = c(1.1, 2.2, 3.3),
    UMAP_2 = c(0.1, -0.2, 0.4),
    stringsAsFactors = FALSE
  )
  emb <- e$extract_embedding_from_metadata(md, "UMAP_1", "UMAP_2",
                                           cell_id_col = "cell_id")
  expect_equal(nrow(emb), 3L)
  expect_equal(emb$cell_id, paste0("C", 1:3))
  expect_equal(emb$dim1, c(1.1, 2.2, 3.3))
  expect_equal(emb$dim2, c(0.1, -0.2, 0.4))
  expect_equal(unique(emb$dim1_name), "UMAP_1")
  expect_equal(unique(emb$dim2_name), "UMAP_2")
})

test_that("extract_embedding_from_metadata falls back to .cell_id / first column when cell_id_col is missing", {
  e <- source_shiny_helpers()
  md <- data.frame(
    .cell_id = paste0("C", 1:3),
    UMAP_1 = 1:3,
    UMAP_2 = 4:6,
    stringsAsFactors = FALSE
  )
  emb <- e$extract_embedding_from_metadata(md, "UMAP_1", "UMAP_2",
                                           cell_id_col = "(none)")
  expect_equal(emb$cell_id, paste0("C", 1:3))
})

# ---- celltype_pairwise_pvalues -------------------------------------------

test_that("celltype_pairwise_pvalues returns significant pairs only with bracket coords", {
  e <- source_shiny_helpers()
  set.seed(123)
  df <- data.frame(
    Score = c(rnorm(40, mean = -1, sd = 0.4),  # CT_A — clearly negative
              rnorm(40, mean =  1, sd = 0.4),  # CT_B — clearly positive
              rnorm(40, mean =  0, sd = 0.4)), # CT_C — middle
    Cell_type = rep(c("CT_A", "CT_B", "CT_C"), each = 40L),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_pairwise_pvalues(df, "Score", "Cell_type")

  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) >= 1L)
  expect_setequal(c("cell_a", "cell_b", "p_val", "label",
                    "xmin", "xmax", "y_pos"),
                  intersect(c("cell_a", "cell_b", "p_val", "label",
                              "xmin", "xmax", "y_pos"), colnames(out)))
  # CT_A vs CT_B should be the most significant pair (different means).
  expect_true(any(
    (out$cell_a == "CT_A" & out$cell_b == "CT_B") |
      (out$cell_a == "CT_B" & out$cell_b == "CT_A")
  ))
  # Bracket y-positions strictly increase so they don't overlap.
  expect_true(all(diff(out$y_pos) >= 0))
  expect_equal(attr(out, "cell_levels"), c("CT_A", "CT_B", "CT_C"))
})

test_that("celltype_pairwise_pvalues returns NULL when no pair is significant", {
  e <- source_shiny_helpers()
  set.seed(7)
  df <- data.frame(
    Score = rnorm(60),
    Cell_type = rep(c("CT_A", "CT_B"), each = 30L),
    stringsAsFactors = FALSE
  )
  expect_null(e$celltype_pairwise_pvalues(df, "Score", "Cell_type",
                                          p_threshold = 1e-12))
})

test_that("celltype_pairwise_pvalues handles single cell type / sparse data", {
  e <- source_shiny_helpers()
  df <- data.frame(
    Score = rnorm(10), Cell_type = rep("CT_only", 10L),
    stringsAsFactors = FALSE
  )
  expect_null(e$celltype_pairwise_pvalues(df, "Score", "Cell_type"))
})

test_that("celltype_pairwise_pvalues honours a caller-supplied cell_levels order", {
  e <- source_shiny_helpers()
  set.seed(321)
  df <- data.frame(
    Score = c(rnorm(40, mean = -1, sd = 0.4),
              rnorm(40, mean =  1, sd = 0.4),
              rnorm(40, mean =  0, sd = 0.4)),
    Cell_type = rep(c("CT_A", "CT_B", "CT_C"), each = 40L),
    stringsAsFactors = FALSE
  )
  custom_order <- c("CT_B", "CT_C", "CT_A")
  out <- e$celltype_pairwise_pvalues(df, "Score", "Cell_type",
                                     cell_levels = custom_order)
  expect_equal(attr(out, "cell_levels"), custom_order)
  # xmin/xmax must reference the custom positions (CT_B = 1, CT_A = 3).
  ab <- out[(out$cell_a == "CT_A" & out$cell_b == "CT_B") |
              (out$cell_a == "CT_B" & out$cell_b == "CT_A"), , drop = FALSE]
  expect_equal(sort(c(ab$xmin[1L], ab$xmax[1L])), c(1L, 3L))
})

# ---- celltype_source_pvalues ---------------------------------------------

test_that("celltype_source_pvalues drops non-significant brackets by default", {
  e <- source_shiny_helpers()
  set.seed(11)
  # CT_A: large source effect (Tumor >> Normal) -> significant
  # CT_B: no source effect                       -> ns, should be dropped
  df <- data.frame(
    Score     = c(rnorm(40, -1, 0.3), rnorm(40, 1, 0.3),
                  rnorm(40,  0, 0.5), rnorm(40, 0, 0.5)),
    Cell_type = rep(c("CT_A", "CT_A", "CT_B", "CT_B"), each = 40L),
    Source    = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = 40L),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_source_pvalues(df, "Score", "Cell_type", "Source")
  expect_true(nrow(out) >= 1L)
  expect_true(all(out$p_val < 0.05))
  expect_false("CT_B" %in% out$Cell_type)

  # significant_only = FALSE keeps the ns row
  out2 <- e$celltype_source_pvalues(df, "Score", "Cell_type", "Source",
                                    significant_only = FALSE)
  expect_true("CT_B" %in% out2$Cell_type)
})

test_that("celltype_anova_pvalue returns one global F-test bracket spanning all cell types", {
  e <- source_shiny_helpers()
  set.seed(101)
  df <- data.frame(
    Score = c(rnorm(40, -1, 0.4),
              rnorm(40,  1, 0.4),
              rnorm(40,  0, 0.4)),
    Cell_type = rep(c("CT_A", "CT_B", "CT_C"), each = 40L),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_anova_pvalue(df, "Score", "Cell_type")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_equal(out$xmin, 1L)
  expect_equal(out$xmax, 3L)
  expect_lt(out$p_val, 1e-10)
  # New label format: bare "ANOVA (p = ...)" with no significance
  # stars and no significance/connector bracket -- the renderer
  # draws this as a centred text annotation.
  expect_match(out$label, "^ANOVA \\(p = ")
  expect_false(grepl("\\*", out$label))
  expect_equal(attr(out, "render_kind"), "annotation")
  expect_true("x_mid" %in% colnames(out))
  expect_equal(attr(out, "cell_levels"), c("CT_A", "CT_B", "CT_C"))
})

test_that("celltype_anova_pvalue returns NULL when there is only one cell type", {
  e <- source_shiny_helpers()
  df <- data.frame(Score = rnorm(20), Cell_type = rep("CT_A", 20L))
  expect_null(e$celltype_anova_pvalue(df, "Score", "Cell_type"))
})

test_that("celltype_anova_pvalue honours caller-supplied cell_levels order", {
  e <- source_shiny_helpers()
  set.seed(202)
  df <- data.frame(
    Score = c(rnorm(40, -1, 0.4),
              rnorm(40,  1, 0.4),
              rnorm(40,  0, 0.4)),
    Cell_type = rep(c("CT_A", "CT_B", "CT_C"), each = 40L),
    stringsAsFactors = FALSE
  )
  custom_order <- c("CT_B", "CT_C", "CT_A")
  out <- e$celltype_anova_pvalue(df, "Score", "Cell_type",
                                 cell_levels = custom_order)
  expect_equal(attr(out, "cell_levels"), custom_order)
  expect_equal(out$xmin, 1L)
  expect_equal(out$xmax, 3L)
})

test_that("celltype_source_anova_pvalues runs per-cell-type F-tests for 3+ sources", {
  e <- source_shiny_helpers()
  set.seed(303)
  # CT_A: 3 separated source means -> significant ANOVA
  # CT_B: 3 identical source means -> NS ANOVA
  df <- data.frame(
    Score = c(rnorm(40, -1, 0.3), rnorm(40,  1, 0.3), rnorm(40, 0, 0.3),
              rnorm(40,  0, 0.5), rnorm(40,  0, 0.5), rnorm(40, 0, 0.5)),
    Cell_type = rep(c("CT_A", "CT_B"), each = 120L),
    Source    = rep(rep(c("S1", "S2", "S3"), each = 40L), times = 2L),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_source_anova_pvalues(df, "Score", "Cell_type", "Source",
                                         significant_only = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_true(all(out$n_sources == 3L))
  expect_lt(out$p_val[out$Cell_type == "CT_A"], 0.001)
  expect_gt(out$p_val[out$Cell_type == "CT_B"], 0.05)
})

test_that("celltype_source_anova_pvalues drops non-significant brackets by default", {
  e <- source_shiny_helpers()
  set.seed(404)
  df <- data.frame(
    Score = c(rnorm(40, -1, 0.3), rnorm(40,  1, 0.3), rnorm(40, 0, 0.3),
              rnorm(40,  0, 0.5), rnorm(40,  0, 0.5), rnorm(40, 0, 0.5)),
    Cell_type = rep(c("CT_A", "CT_B"), each = 120L),
    Source    = rep(rep(c("S1", "S2", "S3"), each = 40L), times = 2L),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_source_anova_pvalues(df, "Score", "Cell_type", "Source")
  expect_equal(nrow(out), 1L)
  expect_equal(out$Cell_type, "CT_A")
})

test_that("celltype_source_anova_pvalues skips cell types with <2 sources", {
  e <- source_shiny_helpers()
  df <- data.frame(
    Score     = rnorm(80),
    Cell_type = rep(c("CT_A", "CT_B"), each = 40L),
    Source    = c(rep("S1", 40L), rep(c("S1", "S2", "S3"), length.out = 40L)),
    stringsAsFactors = FALSE
  )
  out <- e$celltype_source_anova_pvalues(df, "Score", "Cell_type", "Source",
                                         significant_only = FALSE)
  expect_false("CT_A" %in% out$Cell_type)
})

# ---- .plot_radius_unit option-driven default --------------------------------

test_that(".plot_radius_unit picks up the phenomapr.plot_radius_pt option", {
  e <- source_shiny_helpers()
  old <- getOption("phenomapr.plot_radius_pt")
  on.exit(options(phenomapr.plot_radius_pt = old), add = TRUE)

  options(phenomapr.plot_radius_pt = 8)
  expect_equal(as.numeric(e$.plot_radius_unit()), 8)

  options(phenomapr.plot_radius_pt = NULL)
  expect_equal(as.numeric(e$.plot_radius_unit(3)), 3)

  options(phenomapr.plot_radius_pt = "garbage")
  expect_equal(as.numeric(e$.plot_radius_unit(3)), 3)

  options(phenomapr.plot_radius_pt = -1)
  expect_equal(as.numeric(e$.plot_radius_unit(3)), 3)
})

# ---- plot-download modal helpers ----------------------------------------

test_that("phenomapr_modal_dl_btn renders an action-button-shaped button", {
  e <- source_shiny_helpers()
  html <- as.character(htmltools::renderTags(
    e$phenomapr_modal_dl_btn("score_dist_plot")
  )$html)
  expect_match(html, 'id="score_dist_plot_modal_btn"', fixed = TRUE)
  expect_match(html, "action-button", fixed = TRUE)
  expect_match(html, "phenomapr-panel-download-btn", fixed = TRUE)
})

test_that("phenomapr_card_header_modal_dl puts title + button together", {
  e <- source_shiny_helpers()
  html <- as.character(htmltools::renderTags(
    e$phenomapr_card_header_modal_dl(
      shiny::tags$strong("Score distribution"),
      panel_id = "score_dist_plot"
    )
  )$html)
  expect_match(html, "phenomapr-card-header-dl-title", fixed = TRUE)
  expect_match(html, "Score distribution", fixed = TRUE)
  expect_match(html, "score_dist_plot_modal_btn", fixed = TRUE)
})

test_that("phenomapr_plot_download_modal contains all the customisation inputs", {
  e <- source_shiny_helpers()
  html <- as.character(htmltools::renderTags(
    e$phenomapr_plot_download_modal(
      panel_label = "Score distribution",
      defaults    = list(width = 9, height = 6, dpi = 200,
                         format = "pdf", base_size = 16)
    )
  )$html)
  for (needle in c('id="plot_dl_width"', 'id="plot_dl_height"',
                   'id="plot_dl_dpi"', 'id="plot_dl_format"',
                   'id="plot_dl_base_size"', 'id="plot_dl_action"',
                   "Download plot: Score distribution")) {
    expect_match(html, needle, fixed = TRUE)
  }
  # Defaults are pre-populated.
  expect_match(html, 'value="9"',  fixed = TRUE)
  expect_match(html, 'value="6"',  fixed = TRUE)
  expect_match(html, 'value="200"', fixed = TRUE)
})

test_that("phenomapr_plot_download_modal includes the live preview pane", {
  e <- source_shiny_helpers()
  html <- as.character(htmltools::renderTags(
    e$phenomapr_plot_download_modal(panel_label = "Score distribution")
  )$html)
  expect_match(html, 'id="plot_dl_preview"',           fixed = TRUE)
  expect_match(html, "phenomapr-plot-dl-preview",      fixed = TRUE)
  expect_match(html, "phenomapr-plot-dl-preview-title", fixed = TRUE)
  expect_match(html, "phenomapr-plot-dl-preview-frame", fixed = TRUE)
  expect_match(html, "Updates live as you change",     fixed = TRUE)
  # The modal is now size = "l" to give the preview pane breathing room.
  expect_match(html, "modal-lg", fixed = TRUE)
})

test_that("phenomapr_plot_download_modal has the corner-sharpness slider", {
  e <- source_shiny_helpers()
  html <- as.character(htmltools::renderTags(
    e$phenomapr_plot_download_modal(
      panel_label = "Score distribution",
      defaults    = list(radius_pt = 5)
    )
  )$html)
  expect_match(html, 'id="plot_dl_radius_pt"', fixed = TRUE)
  expect_match(html, "corner sharpness",       ignore.case = TRUE)
  expect_match(html, 'data-from="5"',          fixed = TRUE)
})

test_that(".apply_chicklet_radius mutates ggchicklet2 radius in place", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggchicklet2")
  e <- source_shiny_helpers()
  df <- data.frame(x = c("a", "b", "c"), y = c(1, 2, 3))
  make_p <- function() {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = x)) +
      ggchicklet2::geom_chicklet_bar(
        stat = "identity",
        radius = grid::unit(3, "pt")
      )
  }
  p <- make_p()
  expect_equal(format(p$layers[[1]]$geom_params$radius), "3points")
  e$.apply_chicklet_radius(p, 7)
  expect_equal(format(p$layers[[1]]$geom_params$radius), "7points")
})

test_that(".apply_chicklet_radius is a no-op on invalid inputs", {
  e <- source_shiny_helpers()
  # Non-ggplot input returned unchanged.
  expect_identical(e$.apply_chicklet_radius(list(a = 1), 5), list(a = 1))
  expect_identical(e$.apply_chicklet_radius(NULL, 5), NULL)

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggchicklet2")
  df <- data.frame(x = c("a", "b"), y = c(1, 2))
  for (bad in list(-1, NA_real_, "ten", c(2, 3), NaN, Inf)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = x)) +
      ggchicklet2::geom_chicklet_bar(
        stat = "identity", radius = grid::unit(3, "pt")
      )
    e$.apply_chicklet_radius(p, bad)
    expect_equal(format(p$layers[[1]]$geom_params$radius), "3points")
  }
})

test_that(".phenomapr_save_plot honours radius_pt arg for ggchicklet2 plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggchicklet2")
  e <- source_shiny_helpers()
  df <- data.frame(x = c("a", "b", "c"), y = c(1, 2, 3))
  gg_factory <- function() {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = x)) +
      ggchicklet2::geom_chicklet_bar(
        stat = "identity", radius = grid::unit(3, "pt")
      )
  }
  tmp_low  <- tempfile(fileext = ".png")
  tmp_high <- tempfile(fileext = ".png")
  on.exit({
    if (file.exists(tmp_low))  file.remove(tmp_low)
    if (file.exists(tmp_high)) file.remove(tmp_high)
  }, add = TRUE)
  e$.phenomapr_save_plot(tmp_low,  gg_factory(), format = "png",
                          width = 4, height = 3, dpi = 100, radius_pt = 0)
  e$.phenomapr_save_plot(tmp_high, gg_factory(), format = "png",
                          width = 4, height = 3, dpi = 100, radius_pt = 12)
  expect_true(file.exists(tmp_low)  && file.info(tmp_low)$size  > 0L)
  expect_true(file.exists(tmp_high) && file.info(tmp_high)$size > 0L)
  # Different radius should produce different PNG bytes.
  expect_false(identical(
    readBin(tmp_low,  what = "raw", n = file.info(tmp_low)$size),
    readBin(tmp_high, what = "raw", n = file.info(tmp_high)$size)
  ))
})

test_that(".phenomapr_save_plot writes ggplot output for png + pdf + svg", {
  e <- source_shiny_helpers()
  skip_if_not_installed("ggplot2")
  gg <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()
  # SVG output goes through ggplot2::ggsave -> svglite, which is an
  # optional ggplot2 dependency. CI runners (including the
  # test-coverage workflow on Ubuntu) do not always have svglite
  # installed, and ggsave aborts hard when it is not available.
  # Drop the svg branch in that case so we still exercise png/pdf
  # everywhere and the full svg path locally / in environments that
  # do install svglite.
  fmts <- c("png", "pdf")
  if (requireNamespace("svglite", quietly = TRUE)) {
    fmts <- c(fmts, "svg")
  }
  for (fmt in fmts) {
    tmp <- tempfile(fileext = paste0(".", fmt))
    on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
    e$.phenomapr_save_plot(
      file = tmp, plot_obj = gg, format = fmt,
      width = 4, height = 3, dpi = 100, base_size = 11
    )
    expect_true(file.exists(tmp))
    expect_gt(file.info(tmp)$size, 0L)
  }
})

test_that(".phenomapr_save_plot errors on unsupported format", {
  e <- source_shiny_helpers()
  skip_if_not_installed("ggplot2")
  # The dispatch happens in switch() for non-ggplot inputs; for ggplot
  # input ggsave itself rejects unknown devices. Either way we get an
  # error / condition.
  expect_error(
    e$.phenomapr_save_plot(
      file = tempfile(fileext = ".xyz"),
      plot_obj = stats::runif(5),
      format = "xyz", width = 4, height = 3, dpi = 100
    ),
    "Unsupported|unknown"
  )
})

test_that("celltype_source_pvalues honours a caller-supplied cell_levels order", {
  e <- source_shiny_helpers()
  set.seed(13)
  df <- data.frame(
    Score     = c(rnorm(40, -1, 0.3), rnorm(40, 1, 0.3),
                  rnorm(40,  0, 0.3), rnorm(40, 2, 0.3)),
    Cell_type = rep(c("CT_A", "CT_A", "CT_B", "CT_B"), each = 40L),
    Source    = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = 40L),
    stringsAsFactors = FALSE
  )
  custom_order <- c("CT_B", "CT_A")
  out <- e$celltype_source_pvalues(df, "Score", "Cell_type", "Source",
                                   cell_levels = custom_order,
                                   significant_only = FALSE)
  expect_equal(attr(out, "cell_levels"), custom_order)
  # CT_B should now be at x_pos = 1 instead of 2.
  expect_equal(out$x_pos[out$Cell_type == "CT_B"], 1L)
  expect_equal(out$x_pos[out$Cell_type == "CT_A"], 2L)
})

# ---- patient / sample column detection -----------------------------------

test_that("detect_patient_column finds patient/donor/subject style columns", {
  e <- source_shiny_helpers()
  expect_equal(e$detect_patient_column(data.frame(patient = 1:3)),     "patient")
  expect_equal(e$detect_patient_column(data.frame(patient_id = 1:3)),  "patient_id")
  expect_equal(e$detect_patient_column(data.frame(Patient.ID = 1:3)),  "Patient.ID")
  expect_equal(e$detect_patient_column(data.frame(donor_id = 1:3)),    "donor_id")
  expect_equal(e$detect_patient_column(data.frame(subject = 1:3)),     "subject")
  expect_equal(e$detect_patient_column(data.frame(case_id = 1:3)),     "case_id")
  expect_null(e$detect_patient_column(data.frame(cell_type = letters[1:3])))
  expect_null(e$detect_patient_column(NULL))
})

test_that("detect_sample_column finds sample/library/orig.ident style columns", {
  e <- source_shiny_helpers()
  expect_equal(e$detect_sample_column(data.frame(sample = 1:3)),     "sample")
  expect_equal(e$detect_sample_column(data.frame(sample_id = 1:3)),  "sample_id")
  expect_equal(e$detect_sample_column(data.frame(SampleID = 1:3)),   "SampleID")
  expect_equal(e$detect_sample_column(data.frame(library_id = 1:3)), "library_id")
  expect_equal(e$detect_sample_column(data.frame(orig.ident = letters[1:3])),
               "orig.ident")
  expect_null(e$detect_sample_column(data.frame(cell_type = letters[1:3])))
})

test_that("detect_*_column prefers the canonical column when several match", {
  e <- source_shiny_helpers()
  # patient (more canonical) wins over donor when both are present.
  md <- data.frame(donor_id = 1:3, patient_id = 1:3)
  expect_equal(e$detect_patient_column(md), "patient_id")
  # sample wins over orig.ident.
  md2 <- data.frame(orig.ident = letters[1:3], sample = 1:3)
  expect_equal(e$detect_sample_column(md2), "sample")
})

test_that("count_distinct_meta returns NA when column is missing or empty", {
  e <- source_shiny_helpers()
  md <- data.frame(patient = c("P1", "P1", "P2", NA, ""), stringsAsFactors = FALSE)
  expect_equal(e$count_distinct_meta(md, "patient"), 2L)
  expect_true(is.na(e$count_distinct_meta(md, "nope")))
  expect_true(is.na(e$count_distinct_meta(md, NULL)))
  expect_true(is.na(e$count_distinct_meta(NULL, "patient")))
})

# ---- chicklet wrappers ---------------------------------------------------

test_that(".geom_rounded_* wrappers always return a usable ggplot layer", {
  e <- source_shiny_helpers()
  expect_s3_class(e$.geom_rounded_col(),       "Layer")
  expect_s3_class(e$.geom_rounded_stack(),     "Layer")
  expect_s3_class(e$.geom_rounded_histogram(), "Layer")
  expect_s3_class(e$.geom_rounded_boxplot(),   "Layer")
})

# ---- h5ad reader (Python-free, hdf5r-only) -------------------------------
#
# These tests exercise the post-fix `.read_h5ad()` priority order and the
# `.read_h5ad_hdf5r()` reader. We build a tiny synthetic .h5ad on disk via
# hdf5r itself, mirroring the AnnData v0.8 spec layout the reader claims
# to support:
#
#   * /X group, encoding-type="csr_matrix", with /data, /indices, /indptr
#     and a `shape` attribute (n_obs, n_vars).
#   * /var with `_index` -> gene symbols.
#   * /obs with `_index` -> cell IDs, plus a categorical column (encoded
#     as a sub-group with /categories + /codes) and a string column
#     (plain dataset).
#
# This avoids any Python toolchain in the test pipeline (the entire point
# of the fix), so the tests run wherever hdf5r is installed.
make_synthetic_h5ad <- function(path) {
  # CSR layout: 4 cells x 6 genes, 7 non-zero entries.
  n_obs  <- 4L; n_vars <- 6L
  data    <- c(1, 2, 3, 4, 5, 6, 7)
  indices <- c(0L, 2L, 4L, 1L, 3L, 5L, 0L)  # column indices into n_vars
  indptr  <- c(0L, 1L, 3L, 6L, 7L)
  cell_ids <- paste0("Cell_", seq_len(n_obs))
  gene_ids <- c("GENE_A", "GENE_B", "GENE_C",
                "GENE_D", "GENE_E", "GENE_F")
  cell_types  <- c("TypeA", "TypeB", "TypeA", "TypeC")
  cell_source <- c("Normal", "Tumor", "Tumor", "Normal")

  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  # ---- /X (CSR sparse) ---------------------------------------------------
  X <- h5$create_group("X")
  X[["data"]]    <- as.numeric(data)
  X[["indices"]] <- as.integer(indices)
  X[["indptr"]]  <- as.integer(indptr)
  hdf5r::h5attr(X, "encoding-type")    <- "csr_matrix"
  hdf5r::h5attr(X, "encoding-version") <- "0.1.0"
  hdf5r::h5attr(X, "shape")            <- c(n_obs, n_vars)

  # ---- /var (genes) ------------------------------------------------------
  var <- h5$create_group("var")
  var[["_index"]] <- gene_ids
  hdf5r::h5attr(var, "_index")         <- "_index"
  hdf5r::h5attr(var, "encoding-type")  <- "dataframe"
  hdf5r::h5attr(var, "column-order")   <- character(0)

  # ---- /obs (cells) ------------------------------------------------------
  obs <- h5$create_group("obs")
  obs[["_index"]] <- cell_ids
  hdf5r::h5attr(obs, "_index")        <- "_index"
  hdf5r::h5attr(obs, "encoding-type") <- "dataframe"
  hdf5r::h5attr(obs, "column-order")  <- c("cell_type", "source")
  # plain string column
  obs[["source"]] <- cell_source
  # categorical column (AnnData encodes as group with categories + codes)
  ct_grp <- obs$create_group("cell_type")
  hdf5r::h5attr(ct_grp, "encoding-type")    <- "categorical"
  hdf5r::h5attr(ct_grp, "encoding-version") <- "0.2.0"
  hdf5r::h5attr(ct_grp, "ordered")          <- FALSE
  cats   <- c("TypeA", "TypeB", "TypeC")
  codes  <- match(cell_types, cats) - 1L  # 0-based
  ct_grp[["categories"]] <- cats
  ct_grp[["codes"]]      <- as.integer(codes)

  h5$close_all()
  invisible(path)
}

test_that(".read_h5ad_hdf5r() returns genes x cells dgCMatrix with obs metadata", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  make_synthetic_h5ad(path)

  obj <- e$.read_h5ad_hdf5r(path)
  # Sparse, transposed to genes x cells.
  expect_s4_class(obj, "dgCMatrix")
  expect_equal(dim(obj), c(6L, 4L))
  expect_equal(rownames(obj),
               c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E", "GENE_F"))
  expect_equal(colnames(obj),
               c("Cell_1", "Cell_2", "Cell_3", "Cell_4"))

  # The fixture's CSR triple, re-derived from indptr / indices / data:
  #   Row 1 (indptr 0..1): col 0 -> GENE_A = 1
  #   Row 2 (indptr 1..3): col 2 -> GENE_C = 2, col 4 -> GENE_E = 3
  #   Row 3 (indptr 3..6): col 1 -> GENE_B = 4, col 3 -> GENE_D = 5,
  #                         col 5 -> GENE_F = 6
  #   Row 4 (indptr 6..7): col 0 -> GENE_A = 7
  m <- as.matrix(obj)
  expect_equal(m["GENE_A", "Cell_1"], 1)
  expect_equal(m["GENE_C", "Cell_2"], 2)
  expect_equal(m["GENE_E", "Cell_2"], 3)
  expect_equal(m["GENE_B", "Cell_3"], 4)
  expect_equal(m["GENE_D", "Cell_3"], 5)
  expect_equal(m["GENE_F", "Cell_3"], 6)
  expect_equal(m["GENE_A", "Cell_4"], 7)
  # Cells 1 and 4 each carry only GENE_A; rest of the column is zero.
  expect_equal(sum(m[, "Cell_1"] != 0), 1L)
  expect_equal(sum(m[, "Cell_4"] != 0), 1L)

  # phenomapr_obs is populated from /obs (string + categorical columns).
  obs_df <- attr(obj, "phenomapr_obs")
  expect_s3_class(obs_df, "data.frame")
  expect_equal(nrow(obs_df), 4L)
  expect_setequal(colnames(obs_df), c("cell_type", "source"))
  expect_true(is.factor(obs_df$cell_type))
  expect_setequal(levels(obs_df$cell_type), c("TypeA", "TypeB", "TypeC"))
  expect_equal(as.character(obs_df$cell_type),
               c("TypeA", "TypeB", "TypeA", "TypeC"))
  expect_equal(obs_df$source,
               c("Normal", "Tumor", "Tumor", "Normal"))
  # Rownames are obs_names (real cell IDs), NOT 1:N.
  expect_equal(rownames(obs_df), c("Cell_1", "Cell_2", "Cell_3", "Cell_4"))
})

test_that(".read_h5ad() prefers hdf5r over reticulate (Python-free path)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  make_synthetic_h5ad(path)
  # `.read_h5ad()` is the public dispatcher. With hdf5r installed it
  # MUST short-circuit straight to .read_h5ad_hdf5r() and return a
  # dgCMatrix -- no reticulate, no python.builtin.object, no Python at
  # all. This is the regression we're guarding against (previously
  # reticulate had priority, so a user with both installed would land
  # in the AnnData-Python branch and find_phenotype_markers would
  # later error with "expression must be a matrix...").
  obj <- e$.read_h5ad(path)
  expect_s4_class(obj, "dgCMatrix")
  expect_false(inherits(obj, "python.builtin.object"))
})

test_that(".read_h5ad_hdf5r() exposes /obsm entries via the phenomapr_obsm attribute", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  # Synthetic fixture WITH /obsm entries (the standard make_synthetic_h5ad
  # helper above doesn't write /obsm, so we extend it here).
  e <- source_shiny_helpers()
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  make_synthetic_h5ad(path)

  # Patch /obsm into the existing fixture: AnnData stores obsm matrices
  # as (n_obs x k) but HDF5 row-major writes them transposed on disk
  # (k rows, n_obs cols) -- our reader detects-and-flips both
  # orientations, so we use the "(n_obs, k)" orientation that hdf5r
  # writes natively for an R (k, n_obs) matrix.
  h5 <- hdf5r::H5File$new(path, mode = "r+")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  obsm <- h5$create_group("obsm")
  # Build a UMAP matrix with row i = (i, -i). hdf5r will write it as
  # (n_cells x 2) on disk; our reader fingerprints its shape against
  # n_obs == 4 and accepts the (k=2, n_obs=4) orientation.
  umap_mat <- matrix(c(1, -1, 2, -2, 3, -3, 4, -4), nrow = 4L, byrow = TRUE)
  obsm[["X_umap"]] <- umap_mat
  spatial_mat <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80), nrow = 4L, byrow = TRUE)
  obsm[["spatial"]] <- spatial_mat
  h5$close_all()

  obj <- e$.read_h5ad_hdf5r(path)
  obsm_attr <- attr(obj, "phenomapr_obsm")
  expect_type(obsm_attr, "list")
  expect_setequal(names(obsm_attr), c("X_umap", "spatial"))
  expect_equal(nrow(obsm_attr$X_umap), 4L)
  expect_gte(ncol(obsm_attr$X_umap), 2L)
  expect_equal(rownames(obsm_attr$X_umap), colnames(obj))

  # Public dispatcher: the dropdown should list both X_umap and a
  # collapsed "spatial" entry (with "spatial" pinned to the end).
  avail <- e$list_available_embeddings(obj)
  expect_true("X_umap" %in% avail)
  expect_true("spatial" %in% avail)
  expect_equal(tail(avail, 1L), "spatial")

  # extract_embedding() round-trips both keys to a tidy data.frame.
  emb_umap <- e$extract_embedding(obj, "X_umap")
  expect_s3_class(emb_umap, "data.frame")
  expect_equal(nrow(emb_umap), 4L)
  expect_equal(emb_umap$cell_id, colnames(obj))
  emb_spatial <- e$extract_embedding(obj, "spatial")
  expect_true(all(emb_spatial$is_spatial))
})

test_that("clean_matrix_input() preserves phenomapr_obs / phenomapr_obsm attributes", {
  # Regression: when we switched the h5ad reader to hdf5r-only, embeddings
  # ride along on the matrix as the `phenomapr_obsm` attribute. The Shiny
  # cleanup step (`clean_matrix_input()`) used to strip those attributes
  # because rowsum / sweep create fresh objects, which left users with an
  # empty Reduction dropdown after a clean. The fix is in detect_format.R:
  # capture the attrs up front and re-attach to the cleaned matrix.
  set.seed(7)
  m <- matrix(
    pmax(0, rpois(60, 4)),
    nrow = 6L, ncol = 10L,
    dimnames = list(paste0("G", 1:6), paste0("C", 1:10))
  )
  attr(m, "phenomapr_obs") <- data.frame(
    cell_type = rep(c("A", "B"), 5L),
    row.names = colnames(m), stringsAsFactors = FALSE
  )
  attr(m, "phenomapr_obsm") <- list(
    spatial = matrix(seq_len(20L), nrow = 10L, ncol = 2L,
                     dimnames = list(colnames(m), c("x", "y"))),
    X_umap  = matrix(seq_len(20L) - 1, nrow = 10L, ncol = 2L,
                     dimnames = list(colnames(m), c("UMAP_1", "UMAP_2")))
  )
  res <- clean_matrix_input(m, do_hugo = FALSE, do_collapse_dups = FALSE,
                            do_normalize = TRUE, mode = "single_cell",
                            verbose = FALSE)
  expect_true(is.matrix(res$matrix) || inherits(res$matrix, "Matrix"))
  expect_equal(ncol(res$matrix), 10L)
  obs_back <- attr(res$matrix, "phenomapr_obs")
  obsm_back <- attr(res$matrix, "phenomapr_obsm")
  expect_s3_class(obs_back, "data.frame")
  expect_equal(rownames(obs_back), colnames(res$matrix))
  expect_type(obsm_back, "list")
  expect_setequal(names(obsm_back), c("spatial", "X_umap"))
})

test_that("extract_object_metadata() surfaces phenomapr_obs from hdf5r-loaded h5ad", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  make_synthetic_h5ad(path)
  obj <- e$.read_h5ad_hdf5r(path)

  md <- e$extract_object_metadata(obj)
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), 4L)
  # `.cell_id` round-trips against colnames(obj) -- the rowname-fix
  # regression test (we previously got 1:N rownames so .cell_id was
  # numeric strings instead of cell IDs).
  expect_true(".cell_id" %in% colnames(md))
  expect_equal(md$.cell_id, colnames(obj))
  # AnnData /obs columns come through.
  expect_true(all(c("cell_type", "source") %in% colnames(md)))
})

# ---- spatial / segmentation embedding detection -------------------------
#
# Three behaviours we lock down here:
#
#   1. list_available_embeddings() does NOT advertise "segmentation"
#      for plain (non-spatial) Seurat / SCE / matrix inputs -- the
#      dropdown stays clean for the common case.
#   2. AnnData (python.builtin.object): if any obsm key matches the
#      segmentation conventions ("segmentation" / "X_segmentation" /
#      "X_segmentation_centroid"), we surface a single "segmentation"
#      entry (collapsing whichever variant the user happened to pick).
#   3. extract_embedding(adata, "segmentation") reads back the right
#      obsm matrix even when the underlying key is "X_segmentation"
#      rather than the bare "segmentation".

test_that("list_available_embeddings() reports segmentation only when present", {
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  set.seed(1)
  m <- Matrix::Matrix(
    matrix(rpois(50, 1), nrow = 5,
           dimnames = list(paste0("G", 1:5), paste0("C", 1:10))),
    sparse = TRUE
  )
  # Plain matrix has no embeddings at all -- segmentation must not be
  # synthesised out of thin air.
  expect_false("segmentation" %in% e$list_available_embeddings(m))
  # Plain dgCMatrix
  expect_equal(e$list_available_embeddings(m), character(0))
})

test_that("list_available_embeddings() surfaces segmentation for AnnData with obsm['X_segmentation']", {
  skip_if_not_installed("reticulate")
  ad <- tryCatch(reticulate::import("anndata", convert = FALSE),
                 error = function(e) NULL)
  skip_if(is.null(ad),
          "Python anndata module not reachable through reticulate")
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  e <- source_shiny_helpers()
  n_cells <- 8L; n_genes <- 6L
  X <- matrix(seq_len(n_cells * n_genes), n_cells, n_genes,
              dimnames = list(paste0("Cell_", seq_len(n_cells)),
                              paste0("Gene_", seq_len(n_genes))))
  Xpy <- np$asarray(X, dtype = "float32")
  obs_idx <- paste0("Cell_", seq_len(n_cells))
  var_idx <- paste0("Gene_", seq_len(n_genes))
  obs_df <- pd$DataFrame(list(), index = obs_idx)
  var_df <- pd$DataFrame(list(), index = var_idx)
  adata <- ad$AnnData(X = Xpy, obs = obs_df, var = var_df)

  # Squidpy / scanpy-style spatial centroids in obsm["spatial"].
  spatial_xy <- np$asarray(matrix(runif(n_cells * 2L, 0, 100),
                                  n_cells, 2L),
                           dtype = "float32")
  reticulate::py_set_item(adata$obsm, "spatial", spatial_xy)
  # User stored their cell-segmentation centroids under the Scanpy
  # X_-prefixed convention.
  seg_xy <- np$asarray(matrix(runif(n_cells * 2L, 0, 100),
                              n_cells, 2L),
                       dtype = "float32")
  reticulate::py_set_item(adata$obsm, "X_segmentation", seg_xy)

  avail <- e$list_available_embeddings(adata)
  expect_true("spatial" %in% avail)
  # Surfaced as the bare label "segmentation", regardless of which
  # underlying obsm key (segmentation / X_segmentation / ...) was used.
  expect_true("segmentation" %in% avail)
  # The raw "X_segmentation" obsm name is hidden so the dropdown shows
  # one entry per spatial concept, not two near-duplicates.
  expect_false("X_segmentation" %in% avail)
})

test_that("list_available_embeddings() surfaces segmentation for CosMx obs columns", {
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  n <- 20L
  cells <- paste0("C", seq_len(n))
  m <- Matrix::Matrix(
    matrix(rpois(10L * n, 1), 10L, n,
           dimnames = list(paste0("G", seq_len(10L)), cells)),
    sparse = TRUE
  )
  md <- data.frame(
    x_FOV_px = runif(n, 0, 1000),
    y_FOV_px = runif(n, 0, 1000),
    fov = rep(c("A", "B"), each = 10L),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  spatial_xy <- cbind(
    x_slide_mm = runif(n, 4, 6),
    y_slide_mm = runif(n, -9, -8)
  )
  rownames(spatial_xy) <- cells
  attr(m, "phenomapr_obs") <- md
  attr(m, "phenomapr_obsm") <- list(spatial = spatial_xy)

  avail <- e$list_available_embeddings(m)
  expect_true("spatial" %in% avail)
  expect_true("segmentation" %in% avail)

  emb <- e$extract_embedding(m, "segmentation")
  expect_equal(nrow(emb), n)
  expect_true(all(c("x_FOV_px", "y_FOV_px") %in% c(emb$dim1_name[1L], emb$dim2_name[1L])))
  expect_equal(length(unique(emb$sample)), 2L)
})

test_that("spatial_polygons_available detects phenomapr_polygons on matrix", {
  skip_if_not_installed("Matrix")
  e <- source_shiny_helpers()
  cells <- paste0("C", 1:3)
  m <- Matrix::Matrix(
    matrix(1, 4, 3, dimnames = list(paste0("G", 1:4), cells)),
    sparse = TRUE
  )
  poly <- data.frame(
    cell_id = rep(cells, each = 4L),
    x = c(0, 1, 1, 0, 2, 3, 3, 2, 4, 5, 5, 4),
    y = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
    sample = "fov_1",
    stringsAsFactors = FALSE
  )
  attr(m, "phenomapr_polygons") <- poly
  expect_true(e$spatial_polygons_available(m))
  out <- e$extract_spatial_polygons(m)
  expect_equal(sort(unique(out$cell_id)), sort(cells))
  expect_equal(nrow(out), 12L)
  colored <- e$join_spatial_polygon_colors(
    out,
    data.frame(cell_id = cells, score = 1:3, stringsAsFactors = FALSE),
    "score"
  )
  expect_equal(nrow(colored), 12L)
  expect_true(all(colored$score %in% 1:3))
})

test_that("spatial_polygons_available is false for CosMx centroid-only h5ad", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  path <- path.expand("~/Downloads/seurat_test_core_46.h5ad")
  skip_if_not(file.exists(path), "CosMx demo h5ad not available locally")
  e <- source_shiny_helpers()
  obj <- e$.read_h5ad_hdf5r(path)
  expect_false(e$spatial_polygons_available(obj))
})

test_that("apply_global_spatial_coords_for_combined uses spatial frame", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("hdf5r")
  e <- source_shiny_helpers()
  path <- path.expand("~/Downloads/seurat_test_core_46.h5ad")
  skip_if_not(file.exists(path), "CosMx demo h5ad not available locally")
  obj <- e$.read_h5ad_hdf5r(path)
  seg <- e$extract_embedding(obj, "segmentation")
  spatial <- e$extract_embedding(obj, "spatial")
  combined <- e$apply_global_spatial_coords_for_combined(obj, seg)
  expect_false(e$.spatial_coords_are_fov_local(combined))
  expect_equal(combined$dim1, spatial$dim1[match(combined$cell_id, spatial$cell_id)])
  expect_equal(combined$dim2, spatial$dim2[match(combined$cell_id, spatial$cell_id)])
})

test_that("extract_embedding() resolves 'segmentation' to whichever obsm key the user picked", {
  skip_if_not_installed("reticulate")
  ad <- tryCatch(reticulate::import("anndata", convert = FALSE),
                 error = function(e) NULL)
  skip_if(is.null(ad),
          "Python anndata module not reachable through reticulate")
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  e <- source_shiny_helpers()
  n_cells <- 6L; n_genes <- 4L
  X <- matrix(seq_len(n_cells * n_genes), n_cells, n_genes)
  Xpy <- np$asarray(X, dtype = "float32")
  obs_idx <- paste0("Cell_", seq_len(n_cells))
  var_idx <- paste0("Gene_", seq_len(n_genes))
  obs_df <- pd$DataFrame(list(), index = obs_idx)
  var_df <- pd$DataFrame(list(), index = var_idx)
  adata <- ad$AnnData(X = Xpy, obs = obs_df, var = var_df)

  seg_mat <- matrix(c(seq_len(n_cells), seq_len(n_cells) + 100),
                    n_cells, 2L)
  reticulate::py_set_item(adata$obsm, "X_segmentation",
                          np$asarray(seg_mat, dtype = "float32"))

  emb <- e$extract_embedding(adata, "segmentation")
  expect_s3_class(emb, "data.frame")
  expect_equal(nrow(emb), n_cells)
  expect_true(all(c("cell_id", "dim1", "dim2",
                    "is_spatial", "sample") %in% colnames(emb)))
  # Tissue-frame layout flag piggybacks on `is_spatial` so the plot
  # gets coord_fixed() + reversed Y treatment automatically.
  expect_true(all(emb$is_spatial))
  expect_equal(emb$cell_id, obs_idx)
  expect_equal(emb$dim1, as.numeric(seg_mat[, 1L]))
  expect_equal(emb$dim2, as.numeric(seg_mat[, 2L]))
})

test_that(".read_h5ad_hdf5r() error message no longer pushes users to Python", {
  e <- source_shiny_helpers()
  # Stub requireNamespace to claim hdf5r is missing.
  stub_env <- new.env(parent = environment(e$.read_h5ad_hdf5r))
  stub_env$requireNamespace <- function(pkg, ...) {
    if (identical(pkg, "hdf5r")) FALSE else base::requireNamespace(pkg, ...)
  }
  fn <- e$.read_h5ad_hdf5r
  environment(fn) <- stub_env
  expect_error(
    fn(tempfile(fileext = ".h5ad")),
    "hdf5r"
  )
  # We removed the "Alternatively install Python anndata..." sentence so
  # the friendly error doesn't push users into a Python toolchain.
  err <- tryCatch(fn(tempfile(fileext = ".h5ad")),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("anndata", err, fixed = TRUE))
  expect_false(grepl("py_install", err, fixed = TRUE))
})
