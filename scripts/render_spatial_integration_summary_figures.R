#!/usr/bin/env Rscript
# Manuscript candidate summary figures for PhenoMapR + spatial co-occurrence + CellChat
# integration. Run from package root:
#   Rscript scripts/render_spatial_integration_summary_figures.R
#
# Requires precomputed pair table from:
#   Rscript scripts/render_spatial_colocalization_cellchat.R

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

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
assoc_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_assoc.rds")
if (!file.exists(pairs_rds)) {
  stop("Missing ", pairs_rds, ". Run scripts/render_spatial_colocalization_cellchat.R first.")
}

pair <- readRDS(pairs_rds)
assoc <- if (file.exists(assoc_rds)) readRDS(assoc_rds) else NULL
pg <- pair %>% dplyr::filter(.data$prognostic_pair)

pal_adv <- "#B2182B"
pal_fav <- "#2166AC"
pal_other <- "#bdbdbd"
pal_phys <- "#2166AC"
pal_sec <- "#B2182B"
pal_mixed <- "#666666"

pg <- pg %>%
  dplyr::mutate(
    spatial_score = scale(.data$nhood_z)[, 1] + scale(log1p(pmax(.data$local_cooccur, 0)))[, 1],
    spatial_score = ifelse(is.finite(.data$spatial_score), .data$spatial_score, 0),
    pair_label = paste0(.data$source, " \u2192 ", .data$target),
    cross_type = .data$source_ct != .data$target_ct,
    evidence_tier = dplyr::case_when(
      .data$dual_spatial_lr ~ "Spatial + CellChat",
      .data$spatial_coloc & !.data$cellchat_sig ~ "Spatial only",
      !.data$spatial_coloc & .data$cellchat_sig ~ "CellChat only",
      TRUE ~ "Neither"
    ),
    dominance = factor(
      .data$interaction_dominance,
      levels = c("physical-dominant", "secreted-dominant", "mixed", "none")
    )
  )

.assoc_rho <- function(tbl, label) {
  if (is.null(tbl) || !nrow(tbl)) return(NA_real_)
  hit <- tbl$comparison == label
  if (!any(hit)) return(NA_real_)
  as.numeric(tbl$spearman_rho[hit][1])
}

.fmt_rho <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%.2f", x)
}

n_prognostic <- nrow(pg)
n_dual <- sum(pg$dual_spatial_lr, na.rm = TRUE)
n_cellchat <- sum(pg$cellchat_prob > 0, na.rm = TRUE)
rho_phys_local <- .assoc_rho(assoc, "Physical prob ~ local co-localization")
rho_sec_neighbor <- .assoc_rho(assoc, "Secreted prob ~ pheno vs neighbor rho")
top_cross_dual <- pg %>%
  dplyr::filter(.data$dual_spatial_lr, .data$cross_type) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
  dplyr::slice_head(n = 1)
top_cross_label <- if (nrow(top_cross_dual)) {
  paste0(top_cross_dual$source, " \u2192 ", top_cross_dual$target)
} else {
  "none"
}
has_adv_fav_dual <- any(
  pg$dual_spatial_lr &
    pg$source_pg == "Adverse" & pg$target_pg == "Favorable",
  na.rm = TRUE
) || any(
  pg$dual_spatial_lr &
    pg$source_pg == "Favorable" & pg$target_pg == "Adverse",
  na.rm = TRUE
)

out_dir <- file.path(root, "inst", "figures", "manuscript_candidates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_ms <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(color = "#444444", size = base_size - 1),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold")
    )
}

# --- Version 1: evidence-tier counts + spatial vs CellChat quadrant ---
tier_counts <- pg %>%
  dplyr::count(.data$evidence_tier) %>%
  dplyr::mutate(
    evidence_tier = factor(
      .data$evidence_tier,
      levels = c("Spatial + CellChat", "Spatial only", "CellChat only", "Neither")
    )
  )

p_tier <- ggplot(tier_counts, aes(x = .data$evidence_tier, y = .data$n, fill = .data$evidence_tier)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = .data$n), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = c(
    "Spatial + CellChat" = "#762A83",
    "Spatial only" = "#5AAE61",
    "CellChat only" = "#9970AB",
    "Neither" = "#dddddd"
  )) +
  labs(
    title = "Evidence overlap among prognostic pairs",
    subtitle = "Neighborhood enrichment or local co-localization vs CellChat L-R",
    x = NULL,
    y = "Number of sender \u2192 receiver pairs"
  ) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pg_plot <- pg %>%
  dplyr::filter(.data$cellchat_prob > 0 | .data$spatial_coloc)

p_quad <- ggplot(
  pg_plot,
  aes(x = .data$spatial_score, y = log1p(.data$cellchat_prob), color = .data$dominance)
) +
  geom_point(aes(size = .data$dual_spatial_lr), alpha = 0.7) +
  scale_color_manual(
    values = c(
      "physical-dominant" = pal_phys,
      "secreted-dominant" = pal_sec,
      "mixed" = pal_mixed,
      "none" = "#aaaaaa"
    ),
    name = "Dominant\nCellChat mode"
  ) +
  scale_size_manual(values = c(`TRUE` = 3.2, `FALSE` = 1.8), name = "Spatial +\nCellChat") +
  labs(
    title = "Integrated spatial score vs CellChat support",
    subtitle = "Among prognostic tail pairs with spatial and/or L-R evidence",
    x = "Combined spatial score (z-scaled nhood z + local co-localization)",
    y = "log1p(CellChat aggregated probability)"
  ) +
  theme_ms()

v1 <- patchwork::wrap_plots(p_tier, p_quad, ncol = 1, heights = c(0.9, 1.3))
ggsave(
  file.path(out_dir, "integration_summary_v1_evidence_quadrant.png"),
  v1, width = 9, height = 9, dpi = 300, bg = "white"
)

# --- Version 2: physical vs secreted association panel (manuscript-style) ---
if (!is.null(assoc)) {
  assoc_plot <- assoc %>%
    dplyr::filter(is.finite(.data$spearman_rho)) %>%
    dplyr::mutate(
      mode = dplyr::if_else(grepl("^Physical", .data$comparison), "Physical", "Secreted"),
      metric = dplyr::case_when(
        grepl("local co-localization", .data$comparison) ~ "Local co-localization",
        grepl("nhood z", .data$comparison) ~ "Neighborhood z",
        grepl("pheno vs neighbor", .data$comparison) ~ "PhenoMapR vs neighbor rho",
        grepl("pheno vs co-localization", .data$comparison) ~ "PhenoMapR vs co-localization rho",
        grepl("Physical fraction", .data$comparison) ~ "Physical fraction",
        TRUE ~ "Other"
      ),
      metric = factor(
        .data$metric,
        levels = c("Neighborhood z", "Local co-localization", "PhenoMapR vs neighbor rho", "Physical fraction")
      )
    ) %>%
    dplyr::filter(.data$metric != "Other")

  p_assoc <- ggplot(assoc_plot, aes(x = .data$metric, y = .data$spearman_rho, fill = .data$mode)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#666666") +
    scale_fill_manual(values = c(Physical = pal_phys, Secreted = pal_sec)) +
    labs(
      title = "Interaction mode tracks spatial proximity differently",
      subtitle = paste0(
        "Spearman rho among prognostic pairs with CellChat support (n = ",
        n_cellchat, ")"
      ),
      x = NULL,
      y = "Spearman rho",
      fill = NULL
    ) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  dom_box <- pg %>%
    dplyr::filter(.data$interaction_dominance %in% c("physical-dominant", "secreted-dominant", "mixed")) %>%
    tidyr::pivot_longer(
      cols = c(.data$nhood_z, .data$local_cooccur),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(.data$metric, nhood_z = "Neighborhood z", local_cooccur = "Local co-localization")
    )

  p_dom <- ggplot(dom_box, aes(x = .data$dominance, y = .data$value, fill = .data$dominance)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85, show.legend = FALSE) +
    facet_wrap(~ .data$metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(
      "physical-dominant" = pal_phys,
      "secreted-dominant" = pal_sec,
      "mixed" = pal_mixed
    )) +
    labs(
      title = "Physical-dominant pairs are more spatially enriched",
      x = NULL,
      y = "Value"
    ) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  v2 <- patchwork::wrap_plots(p_assoc, p_dom, ncol = 1, heights = c(1, 1.1))
  ggsave(
    file.path(out_dir, "integration_summary_v2_mode_spatial_assoc.png"),
    v2, width = 9, height = 8.5, dpi = 300, bg = "white"
  )
}

# --- Version 3: top cross-type dual-positive pairs dot matrix ---
top_pairs <- pg %>%
  dplyr::filter(.data$dual_spatial_lr, .data$cross_type) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(pair_label = factor(.data$pair_label, levels = rev(unique(.data$pair_label))))

pair_long <- top_pairs %>%
  tidyr::pivot_longer(
    cols = c(.data$nhood_z, .data$local_cooccur, .data$cellchat_prob_physical, .data$cellchat_frac_physical),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      .data$metric,
      levels = c("nhood_z", "local_cooccur", "cellchat_prob_physical", "cellchat_frac_physical"),
      labels = c("Nhood z", "Local co-localization", "Physical CellChat prob", "Physical fraction")
    ),
    value_scaled = as.numeric(scale(.data$value))
  )

p_matrix <- ggplot(pair_long, aes(x = .data$metric, y = .data$pair_label, size = abs(.data$value_scaled), color = .data$value_scaled)) +
  geom_point() +
  scale_color_gradient2(low = pal_fav, mid = "#f7f7f7", high = pal_adv, midpoint = 0, name = "Scaled\nvalue") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(
    title = "Top cross-type pairs with spatial co-localization and CellChat support",
    subtitle = paste0("Strongest integrated example: ", top_cross_label),
    x = NULL,
    y = NULL
  ) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(out_dir, "integration_summary_v3_top_cross_type_pairs.png"),
  p_matrix, width = 9, height = 6, dpi = 300, bg = "white"
)

# --- Version 4: compact multi-panel "story" figure (manuscript panel h candidate) ---
findings <- tibble::tibble(
  finding = c(
    "1. Prognostic extremes form like-with-like spatial niches",
    paste0(
      "2. ", n_dual, " / ", n_prognostic,
      " prognostic pairs show spatial + spatial CellChat support"
    ),
    if (has_adv_fav_dual) {
      "3. Adverse\u2013favorable dual-positive pairs present on this sample"
    } else {
      "3. No adverse\u2013favorable dual-positive pairs on this sample"
    },
    paste0(
      "4. Physical L-R tracks local co-localization (\u03c1 \u2248 ",
      .fmt_rho(rho_phys_local), ")"
    ),
    paste0(
      "5. Secreted L-R tracks neighbor PhenoMapR concordance (\u03c1 \u2248 ",
      .fmt_rho(rho_sec_neighbor), ")"
    ),
    paste0("6. ", top_cross_label, ": top cross-type integrated candidate")
  ),
  y = seq(6, 1)
)

p_findings <- ggplot(findings, aes(x = 0, y = .data$y, label = .data$finding)) +
  geom_text(hjust = 0, size = 3.6, lineheight = 0.95) +
  xlim(-0.05, 1) +
  labs(title = "Major integration findings", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0))

mini_flow <- tibble::tibble(
  step = factor(
    c("PhenoMapR\ntail labels", "Spatial co-\noccurrence", "CellChat\nL-R inference", "Integrated\nevidence"),
    levels = c("PhenoMapR\ntail labels", "Spatial co-\noccurrence", "CellChat\nL-R inference", "Integrated\nevidence")
  ),
  x = 1:4,
  y = 1
)
p_flow <- ggplot(mini_flow, aes(x = .data$x, y = .data$y)) +
  geom_tile(aes(fill = .data$step), width = 0.9, height = 0.55, color = "white", linewidth = 1) +
  geom_text(aes(label = .data$step), size = 3.2, lineheight = 0.9) +
  geom_segment(
    data = data.frame(x = c(1.45, 2.45, 3.45), xend = c(1.55, 2.55, 3.55), y = 1, yend = 1),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  scale_fill_manual(
    values = c(
      "PhenoMapR\ntail labels" = "#deebf7",
      "Spatial co-\noccurrence" = "#c7e9c0",
      "CellChat\nL-R inference" = "#fdd0bc",
      "Integrated\nevidence" = "#e7d4e8"
    )
  ) +
  xlim(0.5, 4.5) +
  ylim(0.7, 1.3) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5))

dual_highlight <- pg %>%
  dplyr::filter(.data$dual_spatial_lr) %>%
  dplyr::mutate(
    pg_combo = paste0(.data$source_pg, " \u2192 ", .data$target_pg),
    highlight = .data$cross_type & .data$source_pg == "Adverse" & .data$target_pg == "Adverse"
  )

p_pg <- ggplot(dual_highlight, aes(x = .data$pg_combo)) +
  geom_bar(fill = "#cccccc") +
  geom_bar(
    data = dplyr::filter(dual_highlight, .data$highlight),
    fill = pal_adv
  ) +
  labs(
    title = "Dual-positive pairs by prognostic combination",
    subtitle = "Red = cross-type adverse\u2013adverse",
    x = NULL,
    y = "Count"
  ) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

v4 <- patchwork::wrap_plots(
  p_flow,
  p_findings,
  p_tier + theme(plot.title = element_text(size = 11)),
  p_pg,
  ncol = 2,
  nrow = 2,
  heights = c(0.55, 1)
) +
  patchwork::plot_annotation(
    title = "Integration of PhenoMapR, spatial co-occurrence, and CellChat (candidate panel)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  file.path(out_dir, "integration_summary_v4_multipanel_story.png"),
  v4, width = 12, height = 9, dpi = 300, bg = "white"
)

# --- Version 5: bubble link diagram (spatial metric vs mode, sized by dual status) ---
pg_cc <- pg %>% dplyr::filter(.data$cellchat_prob > 0)
p_bubble <- ggplot(
  pg_cc,
  aes(
    x = .data$local_cooccur,
    y = .data$cellchat_frac_physical,
    size = .data$cellchat_prob,
    color = .data$source_pg
  )
) +
  geom_point(alpha = 0.55) +
  scale_color_manual(
    values = c(Adverse = pal_adv, Favorable = pal_fav, Other = pal_other),
    name = "Sender\nprognostic group"
  ) +
  scale_size_continuous(range = c(1, 7), name = "CellChat\nprobability") +
  labs(
    title = "Local co-localization vs physical interaction fraction",
    subtitle = "Each point = prognostic sender \u2192 receiver pair with CellChat support",
    x = "Local co-localization score",
    y = "Fraction of CellChat probability from physical interactions"
  ) +
  theme_ms()

ggsave(
  file.path(out_dir, "integration_summary_v5_local_vs_physical_frac.png"),
  p_bubble, width = 8.5, height = 6.5, dpi = 300, bg = "white"
)

# --- Version 6: three-axis summary for top adverse cross-type pairs ---
top_adv <- pg %>%
  dplyr::filter(.data$source_pg == "Adverse", .data$target_pg == "Adverse", .data$cross_type) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::mutate(pair_short = paste0(.data$source_ct, "\u2192", .data$target_ct))

radar_long <- top_adv %>%
  tidyr::pivot_longer(
    cols = c(.data$nhood_z, .data$local_cooccur, .data$cellchat_prob, .data$pheno_neighbor_rho),
    names_to = "axis",
    values_to = "raw"
  ) %>%
  dplyr::group_by(.data$axis) %>%
  dplyr::mutate(norm = (.data$raw - min(.data$raw, na.rm = TRUE)) / pmax(diff(range(.data$raw, na.rm = TRUE)), 1e-9)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    axis = factor(
      .data$axis,
      levels = c("nhood_z", "local_cooccur", "cellchat_prob", "pheno_neighbor_rho"),
      labels = c("Nhood z", "Local co-localization", "CellChat prob", "Pheno vs neighbor \u03c1")
    )
  )

p_radar <- ggplot(radar_long, aes(x = .data$axis, y = .data$norm, group = .data$pair_short, color = .data$pair_short)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    title = "Adverse cross-type pairs: multi-evidence profiles",
    subtitle = "Axes min\u2013max normalized across top 8 adverse sender \u2192 adverse receiver pairs",
    x = NULL,
    y = "Normalized score",
    color = "Pair"
  ) +
  theme_ms() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

ggsave(
  file.path(out_dir, "integration_summary_v6_adverse_cross_type_profiles.png"),
  p_radar, width = 10, height = 6, dpi = 300, bg = "white"
)

message("Wrote manuscript candidate figures to ", out_dir)
message(
  "Files:\n",
  paste(list.files(out_dir, pattern = "integration_summary_v.*\\.png$"), collapse = "\n")
)
