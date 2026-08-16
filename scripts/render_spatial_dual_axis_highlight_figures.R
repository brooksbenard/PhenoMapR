#!/usr/bin/env Rscript
# Figure panels for the two highlight axes:
#   Axis A — Ductal_Adverse → Ductal_Adverse (EGFR clinical highlight)
#   Axis B — Fibroblast_Adverse → Ductal_Adverse (integrin αv/β1–ECM validation)
#
# Run from package root:
#   Rscript scripts/render_spatial_dual_axis_highlight_figures.R
#
# Requires:
#   inst/data/spatial_coloc_cellchat_pairs.rds
#   inst/data/spatial_coloc_cellchat_lr_pairs.rds
#   vignettes/HT270P1-..._processed.rds  (for tissue maps)
#
# Writes under inst/figures/:
#   spatial_axis_ductal_egfr_panel.png
#   spatial_axis_fibro_ductal_ecm_panel.png
#   spatial_axis_dual_highlight_overview.png
#   spatial_axis_ductal_egfr_spatial.png
#   spatial_axis_fibro_ductal_ecm_spatial.png
#   spatial_axis_panels/  (individual subpanels)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[1L]) else normalizePath(".")

if (file.exists(file.path(root, "DESCRIPTION"))) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install patchwork before running this script.")
}

source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_plot_helpers.R"), local = TRUE)

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
lr_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_lr_pairs.rds")
if (!file.exists(pairs_rds) || !file.exists(lr_rds)) {
  stop("Missing pair RDS. Run scripts/render_spatial_colocalization_cellchat.R first.")
}

pair_tbl <- readRDS(pairs_rds)
lr_tbl <- readRDS(lr_rds)

out_dir <- file.path(root, "inst", "figures")
panel_dir <- file.path(out_dir, "spatial_axis_panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

pal_adv <- "#B2182B"
pal_fav <- "#2166AC"
pal_egfr <- "#ae017e"
pal_ecm <- "#005C55"
pal_other <- "#bdbdbd"
pal_phys <- "#2166AC"
pal_sec <- "#B2182B"

theme_axis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(color = "#555555", size = base_size - 1),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold")
    )
}

.get_pair <- function(source, target) {
  hit <- pair_tbl %>% dplyr::filter(.data$source == !!source, .data$target == !!target)
  if (!nrow(hit)) stop("Missing pair: ", source, " → ", target)
  hit[1, , drop = FALSE]
}

.get_lr <- function(source, target) {
  lr_tbl %>% dplyr::filter(.data$source == !!source, .data$target == !!target)
}

.metric_card_df <- function(row, role_label) {
  tibble::tibble(
    metric = c(
      "Local coloc",
      "Nhood z",
      "CellChat Σprob",
      "Physical frac",
      "Integrated score"
    ),
    value = c(
      row$local_cooccur,
      row$nhood_z,
      row$cellchat_prob,
      row$cellchat_frac_physical,
      row$integrated_score
    ),
    label = c(
      sprintf("%.2f", row$local_cooccur),
      sprintf("%.2f%s", row$nhood_z, if (isTRUE(row$nhood_q < 0.05)) "*" else ""),
      sprintf("%.2f", row$cellchat_prob),
      sprintf("%.0f%%", 100 * row$cellchat_frac_physical),
      sprintf("%.1f", row$integrated_score)
    ),
    role = role_label
  )
}

.pathway_bars <- function(lr, mode = NULL, top_n = 8L, title, fill) {
  df <- lr
  if (!is.null(mode)) df <- df %>% dplyr::filter(.data$interaction_mode == !!mode)
  pw <- df %>%
    dplyr::group_by(.data$pathway_name) %>%
    dplyr::summarise(prob = sum(.data$prob), n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$prob)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(pathway_name = factor(.data$pathway_name, levels = rev(.data$pathway_name)))
  ggplot(pw, aes(x = .data$prob, y = .data$pathway_name)) +
    geom_col(fill = fill, width = 0.72) +
    geom_text(aes(label = .data$n), hjust = -0.15, size = 3, color = "#444444") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = "Σ CellChat probability", y = NULL) +
    theme_axis()
}

.lr_bars <- function(lr, pathway_keep = NULL, top_n = 10L, title, fill) {
  df <- lr
  if (!is.null(pathway_keep)) {
    df <- df %>% dplyr::filter(.data$pathway_name %in% pathway_keep)
  }
  top <- df %>%
    dplyr::arrange(dplyr::desc(.data$prob)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      pair_lab = paste0(.data$ligand, " → ", .data$receptor),
      pair_lab = factor(.data$pair_lab, levels = rev(unique(.data$pair_lab)))
    )
  ggplot(top, aes(x = .data$prob, y = .data$pair_lab)) +
    geom_col(fill = fill, width = 0.72) +
    labs(title = title, x = "CellChat probability", y = NULL) +
    theme_axis()
}

.mode_donut_df <- function(lr) {
  lr %>%
    dplyr::group_by(.data$interaction_mode) %>%
    dplyr::summarise(prob = sum(.data$prob), .groups = "drop") %>%
    dplyr::mutate(
      frac = .data$prob / sum(.data$prob),
      label = sprintf("%s\n%.0f%%", .data$interaction_mode, 100 * .data$frac)
    )
}

dd <- .get_pair("Ductal_Adverse", "Ductal_Adverse")
fd <- .get_pair("Fibroblast_Adverse", "Ductal_Adverse")
lr_dd <- .get_lr("Ductal_Adverse", "Ductal_Adverse")
lr_fd <- .get_lr("Fibroblast_Adverse", "Ductal_Adverse")
ecm_pathways <- c("COLLAGEN", "LAMININ", "FN1", "THBS", "PERIOSTIN")

message("Building Axis A (Ductal→Ductal EGFR) panels ...")

p_dd_metrics <- ggplot(
  .metric_card_df(dd, "Axis A"),
  aes(x = .data$metric, y = .data$value)
) +
  geom_col(fill = pal_egfr, width = 0.7) +
  geom_text(aes(label = .data$label), vjust = -0.35, size = 3.2) +
  labs(
    title = "Axis A evidence: Ductal_Adverse → Ductal_Adverse",
    subtitle = "Clinical / signaling highlight · * nhood FDR < 0.05",
    x = NULL,
    y = "Value"
  ) +
  theme_axis() +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

p_dd_pw <- .pathway_bars(
  lr_dd,
  top_n = 8L,
  title = "Top pathways (within adverse ductal)",
  fill = "#7a0177"
)

p_dd_egf <- .lr_bars(
  lr_dd,
  pathway_keep = "EGF",
  top_n = 8L,
  title = "EGFR axis L–R (EGF pathway)",
  fill = pal_egfr
)

p_dd_mode <- ggplot(.mode_donut_df(lr_dd), aes(x = "", y = .data$frac, fill = .data$interaction_mode)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c(physical = pal_phys, secreted = pal_sec, other = pal_other),
    name = NULL
  ) +
  labs(title = "Interaction mode share", x = NULL, y = NULL) +
  theme_axis() +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

axis_a <- (p_dd_metrics | p_dd_mode) / (p_dd_pw | p_dd_egf) +
  patchwork::plot_annotation(
    title = "Axis A — Adverse ductal → adverse ductal (EGFR highlight)",
    subtitle = paste0(
      "Dual+ by neighborhood + local co-localization + CellChat; ",
      "EGF Σprob = ", sprintf("%.3f", sum(lr_dd$prob[lr_dd$pathway_name == "EGF"])),
      " (TGFA / EREG / AREG → EGFR)"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "#444444")
    )
  )

ggsave(
  file.path(out_dir, "spatial_axis_ductal_egfr_panel.png"),
  axis_a, width = 12, height = 9, dpi = 200, bg = "white"
)
ggsave(file.path(panel_dir, "axis_a_metrics.png"), p_dd_metrics, width = 6, height = 4, dpi = 180, bg = "white")
ggsave(file.path(panel_dir, "axis_a_pathways.png"), p_dd_pw, width = 6, height = 4.5, dpi = 180, bg = "white")
ggsave(file.path(panel_dir, "axis_a_egfr_lr.png"), p_dd_egf, width = 6, height = 4.5, dpi = 180, bg = "white")

message("Building Axis B (Fibro→Ductal integrin–ECM) panels ...")

p_fd_metrics <- ggplot(
  .metric_card_df(fd, "Axis B"),
  aes(x = .data$metric, y = .data$value)
) +
  geom_col(fill = pal_ecm, width = 0.7) +
  geom_text(aes(label = .data$label), vjust = -0.35, size = 3.2) +
  labs(
    title = "Axis B evidence: Fibroblast_Adverse → Ductal_Adverse",
    subtitle = "PhenoMapR + spatial validation · nhood FDR not significant",
    x = NULL,
    y = "Value"
  ) +
  theme_axis() +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

p_fd_ecm <- .pathway_bars(
  lr_fd %>% dplyr::filter(.data$pathway_name %in% ecm_pathways),
  top_n = 6L,
  title = "Integrin–ECM pathways (validation signature)",
  fill = pal_ecm
)

p_fd_lr <- .lr_bars(
  lr_fd %>% dplyr::filter(.data$pathway_name %in% c("COLLAGEN", "FN1", "LAMININ", "PERIOSTIN")),
  top_n = 10L,
  title = "Top ECM → integrin L–R pairs",
  fill = "#5a9e8f"
)

p_fd_egf <- .lr_bars(
  lr_fd,
  pathway_keep = "EGF",
  top_n = 6L,
  title = "Secondary EGFR overlay (CAF → ductal)",
  fill = pal_egfr
)

axis_b <- (p_fd_metrics | p_fd_ecm) / (p_fd_lr | p_fd_egf) +
  patchwork::plot_annotation(
    title = "Axis B — Adverse fibroblast → adverse ductal (integrin αv/β1–ECM validation)",
    subtitle = paste0(
      "Strongest cross-type dual+ (local coloc + CellChat); ",
      "ECM pathways dominate; EGF Σprob = ",
      sprintf("%.3f", sum(lr_fd$prob[lr_fd$pathway_name == "EGF"])),
      " links to Axis A"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "#444444")
    )
  )

ggsave(
  file.path(out_dir, "spatial_axis_fibro_ductal_ecm_panel.png"),
  axis_b, width = 12, height = 9, dpi = 200, bg = "white"
)
ggsave(file.path(panel_dir, "axis_b_metrics.png"), p_fd_metrics, width = 6, height = 4, dpi = 180, bg = "white")
ggsave(file.path(panel_dir, "axis_b_ecm_pathways.png"), p_fd_ecm, width = 6, height = 4.5, dpi = 180, bg = "white")
ggsave(file.path(panel_dir, "axis_b_ecm_lr.png"), p_fd_lr, width = 6.5, height = 5, dpi = 180, bg = "white")
ggsave(file.path(panel_dir, "axis_b_egfr_overlay.png"), p_fd_egf, width = 6, height = 4, dpi = 180, bg = "white")

message("Building dual-axis overview ...")

cmp <- dplyr::bind_rows(
  dd %>% dplyr::mutate(axis = "A: Ductal→Ductal\n(EGFR highlight)"),
  fd %>% dplyr::mutate(axis = "B: Fibro→Ductal\n(Integrin–ECM validation)")
) %>%
  tidyr::pivot_longer(
    cols = c(
      .data$local_cooccur, .data$nhood_z, .data$cellchat_prob,
      .data$cellchat_frac_physical, .data$integrated_score
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      .data$metric,
      levels = c(
        "local_cooccur", "nhood_z", "cellchat_prob",
        "cellchat_frac_physical", "integrated_score"
      ),
      labels = c(
        "Local coloc", "Nhood z", "CellChat Σprob",
        "Physical fraction", "Integrated score"
      )
    )
  )

p_cmp <- ggplot(cmp, aes(x = .data$axis, y = .data$value, fill = .data$axis)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ .data$metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(
    "A: Ductal→Ductal\n(EGFR highlight)" = pal_egfr,
    "B: Fibro→Ductal\n(Integrin–ECM validation)" = pal_ecm
  )) +
  labs(
    title = "Dual-axis comparison",
    subtitle = "Shared sample: HTAN PDAC Visium + CytoSPACE + PhenoMapR + spatial CellChat",
    x = NULL,
    y = "Value",
    fill = NULL
  ) +
  theme_axis() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

focus_pw <- dplyr::bind_rows(
  lr_dd %>%
    dplyr::filter(.data$pathway_name %in% c("EGF", "LAMININ", "COLLAGEN", "SEMA3", "CEACAM", "EPHA")) %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(prob = sum(prob), .groups = "drop") %>%
    dplyr::mutate(axis = "A: Ductal→Ductal"),
  lr_fd %>%
    dplyr::filter(.data$pathway_name %in% c("EGF", "LAMININ", "COLLAGEN", "FN1", "PERIOSTIN", "SEMA3", "TGFb")) %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(prob = sum(prob), .groups = "drop") %>%
    dplyr::mutate(axis = "B: Fibro→Ductal")
) %>%
  dplyr::mutate(
    highlight = dplyr::case_when(
      .data$axis == "A: Ductal→Ductal" & .data$pathway_name == "EGF" ~ "Primary (EGFR)",
      .data$axis == "B: Fibro→Ductal" & .data$pathway_name %in% c("COLLAGEN", "LAMININ", "FN1", "PERIOSTIN") ~
        "Primary (ECM)",
      .data$pathway_name == "EGF" ~ "EGFR overlay",
      TRUE ~ "Other"
    )
  )

p_focus <- ggplot(
  focus_pw,
  aes(x = reorder(.data$pathway_name, .data$prob), y = .data$prob, fill = .data$highlight)
) +
  geom_col(width = 0.72) +
  facet_wrap(~ .data$axis, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Primary (EGFR)" = pal_egfr,
    "Primary (ECM)" = pal_ecm,
    "EGFR overlay" = "#fa9fb5",
    "Other" = pal_other
  )) +
  labs(
    title = "Pathway focus by axis role",
    x = NULL,
    y = "Σ CellChat probability",
    fill = NULL
  ) +
  theme_axis() +
  theme(legend.position = "bottom")

ann <- tibble::tibble(
  axis = c("Axis A", "Axis B"),
  role = c("Clinical / signaling", "Method validation"),
  molecular = c("EGFR autocrine (TGFA/EREG/AREG)", "ITGAV/ITGB1 ← COL/FN1/LAM/POSTN"),
  spatial = c(
    sprintf("nhood z=%.2f*; local=%.2f", dd$nhood_z, dd$local_cooccur),
    sprintf("nhood z=%.2f (NS); local=%.2f", fd$nhood_z, fd$local_cooccur)
  )
)

p_ann <- ggplot(ann, aes(x = 0, y = reorder(.data$axis, dplyr::desc(.data$axis)))) +
  geom_tile(aes(fill = .data$axis), width = 2.2, height = 0.85, alpha = 0.15) +
  geom_text(
    aes(label = paste0(.data$axis, "  ·  ", .data$role)),
    x = -0.95, hjust = 0, fontface = "bold", size = 3.4
  ) +
  geom_text(
    aes(label = .data$molecular),
    x = -0.95, vjust = 2.2, hjust = 0, size = 3, color = "#333333"
  ) +
  geom_text(
    aes(label = .data$spatial),
    x = -0.95, vjust = 4.0, hjust = 0, size = 2.8, color = "#666666"
  ) +
  scale_fill_manual(values = c("Axis A" = pal_egfr, "Axis B" = pal_ecm), guide = "none") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0.4, 2.6)) +
  labs(title = "Framing", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12))

overview <- (p_ann / p_cmp / p_focus) +
  patchwork::plot_layout(heights = c(0.55, 1.1, 1.35)) +
  patchwork::plot_annotation(
    title = "Two-axis framing for PhenoMapR spatial PDAC analysis",
    subtitle = paste0(
      "Axis A nominates EGFR inside adverse ductal niches; ",
      "Axis B validates CAF–ductal contact biology among adverse tails"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "#444444")
    )
  )

ggsave(
  file.path(out_dir, "spatial_axis_dual_highlight_overview.png"),
  overview, width = 10, height = 12, dpi = 200, bg = "white"
)

# Tissue spatial maps
vdir <- file.path(root, "vignettes")
rds_cyto <- file.path(vdir, "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

.make_axis_spatial <- function(
    source_lab, target_lab, title,
    highlight_genes_lig, highlight_genes_rec, out_png
) {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    message("Skipping spatial maps (Seurat/Matrix not available).")
    return(invisible(NULL))
  }
  if (!file.exists(rds_cyto)) {
    message("Skipping spatial maps (missing CytoSPACE RDS).")
    return(invisible(NULL))
  }

  message("Loading CytoSPACE object for ", source_lab, " → ", target_lab, " ...")
  seurat <- .spatial_load_cytospace_seurat(rds_cyto, gd_id = gd_id_cyto)
  cyto <- .spatial_cytospace_cell_df(seurat)
  cell_df <- cyto$cell_df
  score_col <- cyto$score_col
  jitter <- .spatial_jitter_params(cell_df, point_range = c(0.25, 1.1), uniform_size = 0.45)

  expr_counts <- tryCatch(
    Seurat::GetAssayData(seurat, assay = "Spatial", layer = "counts"),
    error = function(e) Seurat::GetAssayData(seurat, assay = "Spatial", slot = "counts")
  )
  # Prefer Cell barcodes that match assay colnames; fall back to rownames.
  id_candidates <- list(cell_df$Cell, rownames(cell_df), cell_df$Cell_id)
  cells <- character(0)
  id_col <- NULL
  for (cand in id_candidates) {
    if (is.null(cand)) next
    hit <- intersect(as.character(cand), colnames(expr_counts))
    if (length(hit) > length(cells)) {
      cells <- hit
      id_col <- cand
    }
  }
  if (!length(cells)) {
    stop("No overlapping cell IDs between metadata and Spatial assay")
  }
  keep <- match(cells, as.character(id_col))
  cell_df <- cell_df[keep, , drop = FALSE]
  expr_counts <- expr_counts[, cells, drop = FALSE]
  cs <- Matrix::colSums(expr_counts)
  expr <- Matrix::t(Matrix::t(expr_counts) / pmax(cs, 1) * 10000)
  expr <- log1p(expr)

  .mean_genes <- function(genes) {
    genes <- intersect(genes, rownames(expr))
    if (!length(genes)) return(rep(0, nrow(cell_df)))
    mat <- as.matrix(expr[genes, , drop = FALSE])
    if (nrow(mat) == 1L) as.numeric(mat[1, ]) else as.numeric(Matrix::colMeans(mat))
  }

  cell_df$lig_score <- .mean_genes(highlight_genes_lig)
  cell_df$rec_score <- .mean_genes(highlight_genes_rec)
  cell_df$role <- dplyr::case_when(
    identical(source_lab, target_lab) & cell_df$ct_pg == source_lab ~ "Adverse ductal",
    cell_df$ct_pg == source_lab ~ "Sender",
    cell_df$ct_pg == target_lab ~ "Receiver",
    grepl("_Adverse$", cell_df$ct_pg) ~ "Other Adverse",
    grepl("_Favorable$", cell_df$ct_pg) ~ "Favorable",
    TRUE ~ "Other"
  )

  bg <- cell_df %>% dplyr::filter(.data$role %in% c("Other", "Favorable", "Other Adverse"))
  fg <- cell_df %>% dplyr::filter(!.data$role %in% c("Other", "Favorable", "Other Adverse"))

  p_roles <- ggplot() +
    .spatial_geom_jitter(
      bg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#dddddd", alpha = 0.25
    ) +
    .spatial_geom_jitter(
      fg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$role), alpha = 0.9
    ) +
    scale_color_manual(values = c(
      "Adverse ductal" = pal_adv,
      "Sender" = "#E69F00",
      "Receiver" = pal_adv,
      "Other Adverse" = "#f4a582",
      "Favorable" = pal_fav,
      "Other" = "#dddddd"
    )) +
    coord_equal() +
    labs(title = "Sender / receiver on tissue", color = NULL) +
    .spatial_map_theme(10) +
    theme(legend.position = "bottom")

  focus_cells <- if (identical(source_lab, target_lab)) {
    cell_df %>% dplyr::filter(.data$ct_pg == source_lab)
  } else {
    cell_df %>% dplyr::filter(.data$ct_pg %in% c(source_lab, target_lab))
  }

  p_pheno <- ggplot() +
    .spatial_geom_jitter(
      cell_df %>% dplyr::filter(!.data$ct_pg %in% unique(c(source_lab, target_lab))),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.2
    ) +
    .spatial_geom_jitter(
      focus_cells, order_col = score_col, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data[[score_col]]), alpha = 0.9
    ) +
    scale_color_gradient2(
      low = pal_fav, mid = "#f7f7f7", high = pal_adv,
      midpoint = stats::median(focus_cells[[score_col]], na.rm = TRUE),
      name = "PhenoMapR"
    ) +
    coord_equal() +
    labs(title = "PhenoMapR score (axis cells)") +
    .spatial_map_theme(10) +
    theme(legend.position = "bottom")

  sender_cells <- cell_df %>% dplyr::filter(.data$ct_pg == source_lab)
  receiver_cells <- cell_df %>% dplyr::filter(.data$ct_pg == target_lab)

  p_lig <- ggplot() +
    .spatial_geom_jitter(
      cell_df %>% dplyr::filter(.data$ct_pg != source_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.15
    ) +
    .spatial_geom_jitter(
      sender_cells, order_col = "lig_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$lig_score), alpha = 0.9
    ) +
    scale_color_gradientn(colours = .spatial_cellchat_prob_palette_colors(), name = "Ligand\nexpr") +
    coord_equal() +
    labs(title = paste0("Sender ligands: ", paste(highlight_genes_lig, collapse = "/"))) +
    .spatial_map_theme(10) +
    theme(legend.position = "bottom")

  p_rec <- ggplot() +
    .spatial_geom_jitter(
      cell_df %>% dplyr::filter(.data$ct_pg != target_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.15
    ) +
    .spatial_geom_jitter(
      receiver_cells, order_col = "rec_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$rec_score), alpha = 0.9
    ) +
    scale_color_gradientn(colours = .spatial_cellchat_prob_palette_colors(), name = "Receptor\nexpr") +
    coord_equal() +
    labs(title = paste0("Receiver receptors: ", paste(highlight_genes_rec, collapse = "/"))) +
    .spatial_map_theme(10) +
    theme(legend.position = "bottom")

  spat <- (p_roles | p_pheno) / (p_lig | p_rec) +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(out_png, spat, width = 11, height = 10, dpi = 180, bg = "white")
  invisible(TRUE)
}

.make_axis_spatial(
  source_lab = "Ductal_Adverse",
  target_lab = "Ductal_Adverse",
  title = "Axis A tissue maps — Ductal_Adverse EGFR autocrine niche",
  highlight_genes_lig = c("TGFA", "EREG", "AREG"),
  highlight_genes_rec = c("EGFR", "ERBB2"),
  out_png = file.path(out_dir, "spatial_axis_ductal_egfr_spatial.png")
)

.make_axis_spatial(
  source_lab = "Fibroblast_Adverse",
  target_lab = "Ductal_Adverse",
  title = "Axis B tissue maps — Fibroblast_Adverse → Ductal_Adverse integrin–ECM niche",
  highlight_genes_lig = c("COL1A1", "COL1A2", "FN1", "POSTN"),
  highlight_genes_rec = c("ITGAV", "ITGB1", "ITGA2"),
  out_png = file.path(out_dir, "spatial_axis_fibro_ductal_ecm_spatial.png")
)

message("Done. Dual-axis highlight figures written to ", out_dir)
message("  - spatial_axis_ductal_egfr_panel.png")
message("  - spatial_axis_fibro_ductal_ecm_panel.png")
message("  - spatial_axis_dual_highlight_overview.png")
message("  - spatial_axis_ductal_egfr_spatial.png")
message("  - spatial_axis_fibro_ductal_ecm_spatial.png")
message("  - spatial_axis_panels/ (individual subpanels)")
