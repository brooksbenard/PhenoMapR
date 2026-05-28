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

test_that(".phenomapr_save_plot writes ggplot output for png + pdf + svg", {
  e <- source_shiny_helpers()
  skip_if_not_installed("ggplot2")
  gg <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()
  for (fmt in c("png", "pdf", "svg")) {
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
