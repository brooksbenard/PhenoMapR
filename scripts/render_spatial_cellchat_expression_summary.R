#!/usr/bin/env Rscript
# Summary plots for top spatially co-localized L-R expression axes across prognostic cells.
# CellChat-inspired topic composition + spatial axis maps (from cached pair metrics).
#
# Run from package root:
#   Rscript scripts/render_spatial_cellchat_expression_summary.R
#
# Requires:
#   inst/data/spatial_coloc_cellchat_pairs.rds
#   inst/data/spatial_coloc_cellchat_lr_pairs.rds
# Optional (for spatial axis maps):
#   inst/data/spatial_pair_cell_metrics.rds  (from render_spatial_pair_spatial_maps.R)

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
  library(tidyr)
  library(tibble)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install patchwork for multi-panel summary figures.")
}

source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
lr_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_lr_pairs.rds")
metrics_rds <- file.path(root, "inst", "data", "spatial_pair_cell_metrics.rds")
cellchat_rds <- file.path(root, "inst", "data", "spatial_cellchat_object.rds")

for (f in c(pairs_rds, lr_rds)) {
  if (!file.exists(f)) {
    stop("Missing ", f, ". Run scripts/render_spatial_colocalization_cellchat.R first.")
  }
}

pair_tbl <- readRDS(pairs_rds)
lr_tbl <- readRDS(lr_rds)
metrics <- if (file.exists(metrics_rds)) readRDS(metrics_rds) else NULL

out_dir <- file.path(root, "inst", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal_adv <- "#B2182B"
pal_fav <- "#2166AC"
pal_phys <- "#2166AC"
pal_sec <- "#B2182B"
pal_mixed <- "#9970AB"

dual_pairs <- pair_tbl %>%
  dplyr::filter(.data$dual_spatial_lr) %>%
  dplyr::mutate(
    pair_label = paste0(.data$source, " \u2192 ", .data$target),
    cross_type = .data$source_ct != .data$target_ct
  )

lr_dual <- lr_tbl %>%
  dplyr::inner_join(
    dual_pairs %>% dplyr::select(source, target, source_pg, target_pg, cross_type),
    by = c("source", "target")
  )

ann_vec <- if ("annotation" %in% names(lr_dual)) {
  lr_dual$annotation
} else if ("annotation.x" %in% names(lr_dual)) {
  lr_dual$annotation.x
} else {
  lr_dual$annotation.y
}

mode_vec <- if ("interaction_mode" %in% names(lr_dual) && any(!is.na(lr_dual$interaction_mode))) {
  dplyr::recode(
    lr_dual$interaction_mode,
    physical = "Physical",
    secreted = "Secreted",
    other = "Other",
    .default = lr_dual$interaction_mode
  )
} else {
  dplyr::case_when(
    ann_vec %in% c("Cell-Cell Contact", "ECM-Receptor") ~ "Physical",
    ann_vec == "Secreted Signaling" ~ "Secreted",
    TRUE ~ "Other"
  )
}

lr_dual <- lr_dual %>%
  dplyr::mutate(
    annotation = ann_vec,
    interaction_mode = mode_vec,
    pair_label = paste0(.data$source, " \u2192 ", .data$target)
  )

# --- 1. Topic-style L-R composition (netVisual_TopicComposition analogue) ---
top_lr <- lr_dual %>%
  dplyr::filter(.data$cross_type) %>%
  dplyr::group_by(.data$interaction_name, .data$pathway_name, .data$interaction_mode) %>%
  dplyr::summarise(
    prob_sum = sum(.data$prob, na.rm = TRUE),
    n_pairs = dplyr::n_distinct(.data$pair_label),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$prob_sum)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::mutate(
    lr_label = paste0(.data$interaction_name, "\n(", .data$pathway_name, ")")
  )

top_lr$lr_label <- factor(
  top_lr$lr_label,
  levels = rev(unique(top_lr$lr_label[order(top_lr$prob_sum)]))
)

p_topic <- ggplot(top_lr, aes(x = .data$prob_sum, y = .data$lr_label, fill = .data$interaction_mode)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(Physical = pal_phys, Secreted = pal_sec, Other = pal_mixed, Unknown = "#aaaaaa"),
    name = "Interaction mode"
  ) +
  labs(
    title = "Top spatially co-localized L-R axes",
    subtitle = "Cross-type dual-positive pairs: summed spatial CellChat probability",
    x = "Summed communication probability",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(out_dir, "spatial_coloc_lr_topic_composition.png"),
  p_topic,
  width = 9,
  height = 7,
  dpi = 150,
  bg = "white"
)

# --- 2. Pathway composition by prognostic sender (topic composition across groups) ---
path_comp <- lr_dual %>%
  dplyr::filter(.data$cross_type, .data$source_pg %in% c("Adverse", "Favorable")) %>%
  dplyr::group_by(.data$source_pg, .data$pathway_name, .data$interaction_mode) %>%
  dplyr::summarise(prob_sum = sum(.data$prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(.data$source_pg) %>%
  dplyr::mutate(frac = .data$prob_sum / sum(.data$prob_sum)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(.data$frac > 0.01)

top_pathways <- path_comp %>%
  dplyr::group_by(.data$pathway_name) %>%
  dplyr::summarise(total = sum(.data$prob_sum), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$total)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(.data$pathway_name)

path_comp <- path_comp %>%
  dplyr::filter(.data$pathway_name %in% top_pathways) %>%
  dplyr::mutate(
    pathway_name = factor(.data$pathway_name, levels = rev(top_pathways)),
    source_pg = factor(.data$source_pg, levels = c("Adverse", "Favorable"))
  )

p_path_topic <- ggplot(
  path_comp,
  aes(x = .data$source_pg, y = .data$frac, fill = interaction_mode)
) +
  geom_col(position = "stack", width = 0.6) +
  facet_wrap(~ .data$pathway_name, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c(Physical = pal_phys, Secreted = pal_sec, Other = pal_mixed, Unknown = "#aaaaaa"),
    name = "Mode"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "L-R pathway composition by prognostic sender group",
    subtitle = "Fraction of dual-positive cross-type signaling from adverse vs favorable senders",
    x = "Sender prognostic group",
    y = "Fraction of pathway probability"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 9)
  )

ggsave(
  file.path(out_dir, "spatial_coloc_pathway_topic_composition.png"),
  p_path_topic,
  width = 10,
  height = 8,
  dpi = 150,
  bg = "white"
)

# --- 3. Expression-axis heatmap: top L-R pairs x prognostic cell groups ---
axis_top <- lr_dual %>%
  dplyr::filter(.data$cross_type) %>%
  dplyr::group_by(.data$interaction_name, .data$source) %>%
  dplyr::summarise(axis_score = sum(.data$prob, na.rm = TRUE), .groups = "drop")

top_interactions <- axis_top %>%
  dplyr::group_by(.data$interaction_name) %>%
  dplyr::summarise(total = sum(.data$axis_score), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$total)) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::pull(.data$interaction_name)

axis_mat <- axis_top %>%
  dplyr::filter(.data$interaction_name %in% top_interactions)

if (nrow(axis_mat)) {
  axis_wide <- axis_mat %>%
    dplyr::select(.data$interaction_name, .data$source, .data$axis_score) %>%
    tidyr::pivot_wider(names_from = .data$source, values_from = .data$axis_score, values_fill = 0)
  mat_axis <- as.matrix(axis_wide[, setdiff(names(axis_wide), "interaction_name"), drop = FALSE])
  rownames(mat_axis) <- axis_wide$interaction_name
  mat_axis_log <- log1p(mat_axis)
  axis_long <- as.data.frame(mat_axis_log) %>%
    tibble::rownames_to_column("interaction_name") %>%
    tidyr::pivot_longer(-.data$interaction_name, names_to = "ct_pg", values_to = "log_prob") %>%
    dplyr::mutate(
      interaction_name = factor(
        .data$interaction_name,
        levels = rev(rownames(mat_axis_log)[order(rowSums(mat_axis_log))])
      )
    )

  p_axis_hm <- ggplot(
    axis_long,
    aes(x = .data$ct_pg, y = .data$interaction_name, fill = .data$log_prob)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = .spatial_cellchat_prob_palette_colors(),
      name = "log1p(prob)"
    ) +
    labs(
      title = "Top co-localized expression axes across prognostic groups",
      subtitle = "Rows: L-R pairs among dual-positive cross-type interactions; columns: sender labels",
      x = "Prognostic sender group",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(out_dir, "spatial_coloc_expression_axis_heatmap.png"),
    p_axis_hm,
    width = 11,
    height = 7,
    dpi = 150,
    bg = "white"
  )
}

# --- 4. Spatial axis maps for top dual-positive pairs (spatialTopicPlot / Lee-style) ---
if (!is.null(metrics) && is.list(metrics) && "cells" %in% names(metrics)) {
  cell_df <- metrics$cells
  pair_meta <- metrics$pairs

  theme_sp <- function() {
    theme_void(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9, lineheight = 0.95),
        plot.subtitle = element_text(hjust = 0.5, color = "#555555", size = 8, lineheight = 0.95),
        legend.position = "right",
        legend.key.height = grid::unit(0.35, "cm"),
        plot.margin = margin(2, 2, 2, 2)
      )
  }

  .scale01 <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok)) return(rep(0, length(x)))
    rng <- range(x[ok], na.rm = TRUE)
    if (diff(rng) == 0) return(ifelse(ok, 0.5, 0))
    out <- rep(NA_real_, length(x))
    out[ok] <- (x[ok] - rng[1]) / diff(rng)
    out
  }

  plot_pair_axis <- function(pair_row) {
    slug <- paste0(pair_row$source, "__", pair_row$target)
    slug <- gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", slug))
    sub <- cell_df %>%
      dplyr::filter(.data$pair_slug == slug | (
        .data$pair_source == pair_row$source & .data$pair_target == pair_row$target
      )) %>%
      dplyr::filter(.data$ct_pg == pair_row$source)

    if (!nrow(sub)) return(NULL)

    sub <- sub %>%
      dplyr::mutate(
        axis_score = .scale01(.data$cellchat_sender_score),
        local_axis = .scale01(.data$local_cooccur_score)
      )

    p_lr <- ggplot(sub, aes(x = .data$row, y = -.data$col, color = .data$axis_score)) +
      geom_point(size = 0.35, alpha = 0.85) +
      scale_color_gradientn(
        colours = c("#f7f7f7", "#fdd0bc", "#b2182b"),
        name = "L-R axis\n(sender)"
      ) +
      coord_fixed() +
      labs(
        title = paste0(pair_row$source, " \u2192 ", pair_row$target),
        subtitle = "Sender L-R expression axis on tissue"
      ) +
      theme_sp()

    p_loc <- ggplot(sub, aes(x = .data$row, y = -.data$col, color = .data$local_axis)) +
      geom_point(size = 0.35, alpha = 0.85) +
      .spatial_coloc_score_scale_gg(name = "Local\nco-local.") +
      coord_fixed() +
      labs(
        title = paste0(pair_row$source, " \u2192 ", pair_row$target),
        subtitle = "Local co-localization axis on tissue"
      ) +
      theme_sp()

    patchwork::wrap_plots(p_lr, p_loc, ncol = 2)
  }

  top_sp_pairs <- dual_pairs %>%
    dplyr::filter(.data$cross_type) %>%
    dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
    dplyr::slice_head(n = 4)

  panels <- lapply(seq_len(nrow(top_sp_pairs)), function(i) {
    plot_pair_axis(top_sp_pairs[i, , drop = FALSE])
  })
  panels <- panels[!vapply(panels, is.null, logical(1))]

  if (length(panels)) {
    p_spatial <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_annotation(
        title = "Top spatially co-localized expression axes",
        subtitle = "Dual-positive cross-type pairs: sender L-R axis (left) and local co-localization (right) per row",
        theme = theme(
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, color = "#444444")
        )
      )
    ggsave(
      file.path(out_dir, "spatial_coloc_top_lr_spatial_axes.png"),
      p_spatial,
      width = 9,
      height = 3.2 * length(panels),
      dpi = 150,
      bg = "white",
      limitsize = FALSE
    )
  }
}

# --- 5. Native CellChat spatial plots when cached object is available ---
if (file.exists(cellchat_rds) && requireNamespace("CellChat", quietly = TRUE)) {
  cellchat <- readRDS(cellchat_rds)
  top_path <- lr_dual %>%
    dplyr::filter(.data$cross_type) %>%
    dplyr::group_by(.data$pathway_name) %>%
    dplyr::summarise(s = sum(.data$prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$s)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(.data$pathway_name)

  if (length(top_path) == 1L && nzchar(top_path)) {
    out_net <- file.path(out_dir, "spatial_coloc_cellchat_spatial_network.png")
    message("Writing ", basename(out_net), " via CellChat::netVisual_aggregate ...")
    tryCatch({
      p_net <- CellChat::netVisual_aggregate(
        cellchat,
        signaling = top_path,
        layout = "spatial",
        edge.width.max = 2.5,
        vertex.size.max = 1.2,
        alpha.image = 0.12,
        vertex.label.cex = 0.55,
        point.size = 0.4
      )
      ggplot2::ggsave(
        out_net,
        plot = p_net,
        width = 8,
        height = 7,
        dpi = 150,
        bg = "white"
      )
    }, error = function(e) {
      message("Skipping spatial network plot: ", conditionMessage(e))
    })
  }

  top_lr_name <- top_lr$interaction_name[1]
  if (!is.na(top_lr_name) && nzchar(top_lr_name)) {
    out_feat <- file.path(out_dir, "spatial_coloc_cellchat_spatial_feature.png")
    message("Writing ", basename(out_feat), " via CellChat::spatialFeaturePlot ...")
    tryCatch({
      p_feat <- CellChat::spatialFeaturePlot(
        cellchat,
        pairLR.use = top_lr_name,
        point.size = 0.45,
        color.heatmap = "Reds"
      )
      ggplot2::ggsave(
        out_feat,
        plot = p_feat,
        width = 9,
        height = 4.5,
        dpi = 150,
        bg = "white"
      )
    }, error = function(e) {
      message("Skipping spatial feature plot: ", conditionMessage(e))
    })
  }
}

message("Wrote expression-axis summary figures to ", out_dir)
