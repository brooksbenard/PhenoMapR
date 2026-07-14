#!/usr/bin/env Rscript
# CytoSPACE cell-location overview maps (cell type, PhenoMapR score, prognostic tails).
# Uses row/col from CytoSPACE metadata (Visium spot grid positions for deconvolved cells).
#
# Run from package root:
#   Rscript scripts/render_spatial_cytospace_location_plots.R
#
# Output (scaled + uniform point sizes):
#   inst/figures/spatial_cytospace_locations_scaled.png
#   inst/figures/spatial_cytospace_locations_uniform.png
#   inst/figures/spatial_cytospace_celltype_{scaled,uniform}.png
#   inst/figures/spatial_cytospace_phenomapr_{scaled,uniform}.png
#   inst/figures/spatial_cytospace_prognostic_tails_{scaled,uniform}.png

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[1L]) else normalizePath(".")

if (file.exists(file.path(root, "DESCRIPTION"))) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload or run after devtools::install()")
  }
  pkgload::load_all(root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(PhenoMapR))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Install scales before running this script.", call. = FALSE)
}

source(file.path(root, "scripts", "spatial_plot_helpers.R"), local = TRUE)

vdir <- file.path(root, "vignettes")
rds_cyto <- file.path(vdir, "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

message("Loading CytoSPACE object ...")
seurat <- .spatial_load_cytospace_seurat(rds_cyto, gd_id = gd_id_cyto)
cyto <- .spatial_cytospace_cell_df(seurat)
cell_df <- cyto$cell_df
score_col <- cyto$score_col

jitter_scaled <- .spatial_jitter_params(cell_df, point_range = c(0.5, 1.6), uniform_size = 0.65)

ct_freq <- sort(table(cell_df$CellType, useNA = "no"), decreasing = TRUE)
cell_df$celltype_zorder <- as.numeric(factor(cell_df$CellType, levels = names(ct_freq)))
ct_pal <- PhenoMapR::get_celltype_palette(sort(unique(cell_df$CellType)))

.plot_celltype <- function(point_size_mode, jitter_params) {
  df <- cell_df[order(cell_df$celltype_zorder), , drop = FALSE]
  mapping <- if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    ggplot2::aes(x = .data$row, y = -.data$col, color = .data$CellType, size = .data$points_per_location)
  } else {
    ggplot2::aes(x = .data$row, y = -.data$col, color = .data$CellType)
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$row, y = -.data$col)) +
    .spatial_geom_jitter(
      df,
      mapping = mapping,
      jitter_params = jitter_params,
      point_size_mode = point_size_mode,
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = ct_pal, name = "Cell Type", na.value = "grey90") +
    ggplot2::coord_fixed(ratio = 0.6) +
    .spatial_map_theme(14) +
    ggplot2::labs(title = "CytoSPACE cell types")
  if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    p <- p + ggplot2::scale_size_continuous(
      range = jitter_params$point_range,
      trans = "reverse",
      name = "Cells per spot"
    ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4), ncol = 2))
  }
  p
}

.plot_phenomapr <- function(point_size_mode, jitter_params) {
  df <- cell_df[order(abs(cell_df[[score_col]])), , drop = FALSE]
  mapping <- if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    ggplot2::aes(
      x = .data$row, y = -.data$col,
      color = scale(.data[[score_col]]),
      size = .data$points_per_location
    )
  } else {
    ggplot2::aes(
      x = .data$row, y = -.data$col,
      color = scale(.data[[score_col]])
    )
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$row, y = -.data$col)) +
    .spatial_geom_jitter(
      df,
      mapping = mapping,
      jitter_params = jitter_params,
      point_size_mode = point_size_mode,
      alpha = 0.8
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0,
      trans = scales::pseudo_log_trans(sigma = 0.2),
      name = "PhenoMapR\nz-score"
    ) +
    ggplot2::coord_fixed(ratio = 0.6) +
    .spatial_map_theme(14) +
    ggplot2::labs(title = "PhenoMapR score")
  if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    p <- p + ggplot2::scale_size_continuous(
      range = jitter_params$point_range,
      trans = "reverse",
      name = "Cells per spot"
    )
  }
  p
}

.plot_tails <- function(point_size_mode, jitter_params) {
  df_other <- dplyr::filter(cell_df, .data$prognostic_group == "Other")
  df_extreme <- dplyr::filter(cell_df, .data$prognostic_group != "Other")
  mapping <- if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    ggplot2::aes(
      x = .data$row, y = -.data$col,
      color = .data$prognostic_group,
      size = .data$points_per_location
    )
  } else {
    ggplot2::aes(
      x = .data$row, y = -.data$col,
      color = .data$prognostic_group
    )
  }

  p <- ggplot2::ggplot() +
    .spatial_geom_jitter(
      df_other,
      mapping = mapping,
      jitter_params = jitter_params,
      point_size_mode = point_size_mode,
      alpha = 0.8
    ) +
    .spatial_geom_jitter(
      df_extreme,
      mapping = mapping,
      jitter_params = jitter_params,
      point_size_mode = point_size_mode,
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c(Adverse = "#B2182B", Other = "#f7f7f7", Favorable = "#2166AC"),
      name = "Prognostic group",
      drop = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 0.6) +
    .spatial_map_theme(14) +
    ggplot2::labs(title = "Prognostic tails (5th / 95th percentile)")
  if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    p <- p + ggplot2::scale_size_continuous(
      range = jitter_params$point_range,
      trans = "reverse",
      name = "Cells per spot"
    ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
  }
  p
}

out_dir <- file.path(root, "inst", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (mode in c("scaled", "uniform")) {
  suffix <- if (mode == "scaled") "scaled" else "uniform"
  message("Rendering CytoSPACE location maps (", mode, " points) ...")

  p_ct <- .plot_celltype(mode, jitter_scaled)
  p_pm <- .plot_phenomapr(mode, jitter_scaled)
  p_tl <- .plot_tails(mode, jitter_scaled)

  ggsave(
    file.path(out_dir, paste0("spatial_cytospace_celltype_", suffix, ".png")),
    p_ct, width = 8, height = 7, dpi = 220, bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0("spatial_cytospace_phenomapr_", suffix, ".png")),
    p_pm, width = 8, height = 7, dpi = 220, bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0("spatial_cytospace_prognostic_tails_", suffix, ".png")),
    p_tl, width = 8, height = 7, dpi = 220, bg = "white"
  )

  combo <- (p_ct | p_pm | p_tl) +
    patchwork::plot_annotation(
      title = paste0(
        "CytoSPACE cell locations (",
        if (mode == "scaled") "size scaled by cells per spot" else "uniform point size",
        ")"
      ),
      subtitle = "Visium row/col from CytoSPACE deconvolution | PhenoMapR scored on mapped single cells",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#444444")
      )
    ) +
    patchwork::plot_layout(guides = "collect") &
    theme(aspect.ratio = 1)

  ggsave(
    file.path(out_dir, paste0("spatial_cytospace_locations_", suffix, ".png")),
    combo, width = 22, height = 7.5, dpi = 220, bg = "white"
  )
}

message("Done. CytoSPACE location maps written to ", out_dir)
