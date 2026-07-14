#!/usr/bin/env Rscript
# Spatial maps of PhenoMapR score, local co-localization, and CellChat potential
# for top integrated cell-type pairs (PhenoMapR + spatial co-occurrence + CellChat).
#
# Run from package root:
#   Rscript scripts/render_spatial_pair_spatial_maps.R
#
# Requires:
#   Rscript scripts/render_spatial_colocalization_cellchat.R
#
# Output (scaled + uniform point sizes):
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>.png
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>_uniform.png
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>_compact.png
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>_compact_uniform.png
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>_tails.png
#   inst/figures/spatial_pair_maps/spatial_pair_<source>__<target>_tails_uniform.png
#   inst/data/spatial_pair_cell_metrics.rds

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
  library(Matrix)
  library(Seurat)
  library(tibble)
  library(tidyr)
})

for (pkg in c("patchwork", "RANN")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install ", pkg, " before running this script.", call. = FALSE)
  }
}

source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_plot_helpers.R"), local = TRUE)

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
lr_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_lr_pairs.rds")
if (!file.exists(pairs_rds) || !file.exists(lr_rds)) {
  stop(
    "Missing pair RDS. Run scripts/render_spatial_colocalization_cellchat.R first.",
    call. = FALSE
  )
}

vdir <- file.path(root, "vignettes")
rds_cyto <- file.path(vdir, "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

.pair_slug <- function(source, target) {
  clean <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
  paste0(clean(source), "__", clean(target))
}

.parse_gene_set <- function(gene_str) {
  unique(unlist(strsplit(as.character(gene_str), "_", fixed = TRUE)))
}

.gene_expr_score <- function(expr_mat, gene_str, cells) {
  genes <- .parse_gene_set(gene_str)
  genes <- intersect(genes, rownames(expr_mat))
  if (!length(genes)) return(rep(0, length(cells)))
  mat <- as.matrix(expr_mat[genes, cells, drop = FALSE])
  if (nrow(mat) == 1L) {
    return(as.numeric(mat[1L, ]))
  }
  colMeans(mat)
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

message("Loading CytoSPACE object ...")
seurat <- .spatial_load_cytospace_seurat(rds_cyto, gd_id = gd_id_cyto)
cyto <- .spatial_cytospace_cell_df(seurat)
cell_df <- cyto$cell_df
score_col <- cyto$score_col
spatial_jitter <- .spatial_jitter_params(cell_df, point_range = c(0.35, 1.4), uniform_size = 0.55)

scoc_df <- data.frame(
  x = cell_df$row,
  y = cell_df$col,
  ct_pg = cell_df$ct_pg,
  stringsAsFactors = FALSE
)
rownames(scoc_df) <- cell_df$Cell

xy_mat <- as.matrix(scoc_df[, c("x", "y"), drop = FALSE])
if (any(duplicated(xy_mat) | duplicated(xy_mat, fromLast = TRUE))) {
  rng <- max(
    diff(range(xy_mat[, 1], na.rm = TRUE)),
    diff(range(xy_mat[, 2], na.rm = TRUE)),
    na.rm = TRUE
  )
  eps <- if (is.finite(rng) && rng > 0) rng * 1e-5 else 1e-6
  set.seed(4L)
  scoc_df$x <- scoc_df$x + stats::runif(nrow(scoc_df), 0, eps)
  scoc_df$y <- scoc_df$y + stats::runif(nrow(scoc_df), 0, eps)
}

k_sc <- min(20L, max(5L, nrow(scoc_df) - 1L))
local_maxnsteps <- 3L
rad_sc <- max(
  3,
  0.04 * max(
    diff(range(scoc_df$x, na.rm = TRUE)),
    diff(range(scoc_df$y, na.rm = TRUE)),
    na.rm = TRUE
  )
)

.build_spatial_neighbors <- function() {
  coords_local <<- as.matrix(scoc_df[, c("x", "y"), drop = FALSE])
  neighbors_local <<- Seurat::FindNeighbors(coords_local, k.param = k_sc, verbose = FALSE)
  adj_local <<- neighbors_local$nn
  rownames(adj_local) <<- colnames(adj_local) <<- rownames(scoc_df)
  degrees_local <<- Matrix::colSums(adj_local) + 1
  res_local <<- RANN::nn2(
    data = coords_local,
    query = coords_local,
    searchtype = "radius",
    radius = rad_sc,
    k = k_sc
  )
  label_vec <<- as.character(scoc_df$ct_pg)
  neighbor_labels <<- vector("list", nrow(scoc_df))
  for (i in seq_len(nrow(scoc_df))) {
    ni <- res_local$nn.idx[i, ]
    ni <- ni[ni > 0 & ni != i]
    neighbor_labels[[i]] <<- if (length(ni)) label_vec[ni] else character(0)
  }
}

.build_spatial_neighbors()

.spatial_cooccur_local_scores <- function(cluster_x, cluster_y) {
  binary <- vapply(
    neighbor_labels,
    function(nl) {
      if (!length(nl)) return(0)
      as.numeric(any(nl == cluster_x) & any(nl == cluster_y))
    },
    numeric(1)
  )
  if (local_maxnsteps <= 0L) return(binary)
  s_norm <- matrix(binary) / degrees_local
  as.numeric((adj_local %*% s_norm) + s_norm)
}

message("Normalizing expression for CellChat-style L-R scoring ...")
expr_counts <- Seurat::GetAssayData(seurat, assay = "Spatial", layer = "counts")
rownames(cell_df) <- cell_df$Cell
cells_all <- intersect(rownames(cell_df), colnames(expr_counts))
cell_df <- cell_df[cells_all, , drop = FALSE]
scoc_df <- scoc_df[cells_all, , drop = FALSE]
rownames(scoc_df) <- rownames(cell_df)
.build_spatial_neighbors()
expr_counts <- expr_counts[, cells_all, drop = FALSE]
cs <- Matrix::colSums(expr_counts)
expr_norm <- Matrix::t(Matrix::t(expr_counts) / pmax(cs, 1) * 10000)
expr_norm <- log1p(expr_norm)

.compute_sender_lr_score <- function(source_label, target_label, lr_rows, top_n = 8L) {
  lr_top <- lr_rows %>%
    dplyr::arrange(dplyr::desc(.data$prob)) %>%
    dplyr::slice_head(n = min(top_n, nrow(lr_rows)))
  if (!nrow(lr_top)) return(rep(0, nrow(cell_df)))

  sender_idx <- which(cell_df$ct_pg == source_label)
  if (!length(sender_idx)) return(rep(0, nrow(cell_df)))

  sender_cells <- cell_df$Cell[sender_idx]
  scores <- rep(0, length(sender_cells))
  names(scores) <- sender_cells

  for (k in seq_len(nrow(lr_top))) {
    lig <- lr_top$ligand[k]
    rec <- lr_top$receptor[k]
    w <- lr_top$prob[k]
    lig_expr <- .gene_expr_score(expr_norm, lig, sender_cells)
    rec_neigh <- vapply(
      sender_cells,
      function(cid) {
        i <- match(cid, rownames(scoc_df))
        ni <- res_local$nn.idx[i, ]
        ni <- ni[ni > 0 & ni != i]
        if (!length(ni)) return(0)
        tgt <- rownames(scoc_df)[ni][label_vec[ni] == target_label]
        if (!length(tgt)) return(0)
        mean(.gene_expr_score(expr_norm, rec, tgt))
      },
      numeric(1)
    )
    scores <- scores + w * lig_expr * rec_neigh
  }

  out <- rep(0, nrow(cell_df))
  out[match(names(scores), cell_df$Cell)] <- as.numeric(scores)
  out
}

.compute_target_neighbor_frac <- function(source_label, target_label) {
  sender_idx <- which(cell_df$ct_pg == source_label)
  out <- rep(0, nrow(cell_df))
  if (!length(sender_idx)) return(out)

  for (idx in sender_idx) {
    cid <- cell_df$Cell[idx]
    i <- match(cid, rownames(scoc_df))
    ni <- res_local$nn.idx[i, ]
    ni <- ni[ni > 0 & ni != i]
    if (!length(ni)) next
    out[idx] <- mean(label_vec[ni] == target_label)
  }
  out
}

pair_tbl <- readRDS(pairs_rds)
lr_tbl <- readRDS(lr_rds)

pairs_to_plot <- pair_tbl %>%
  dplyr::filter(.data$dual_spatial_lr, .data$source_ct != .data$target_ct) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
  dplyr::slice_head(n = 5L)

if (!nrow(pairs_to_plot)) {
  stop("No cross-type dual-positive pairs found in ", pairs_rds)
}

out_dir <- file.path(root, "inst", "figures", "spatial_pair_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_spatial <- function(base_size = 11) {
  .spatial_map_theme(base_size)
}

.order_by_abs <- .spatial_order_by_abs

.geom_spatial_jitter <- function(
  data,
  order_col = NULL,
  mapping = NULL,
  alpha = 0.85,
  point_size_mode = "scaled",
  ...
) {
  if (is.null(mapping)) {
    mapping <- aes(x = .data$row, y = -.data$col)
  }
  .spatial_geom_jitter(
    data = data,
    mapping = mapping,
    order_col = order_col,
    jitter_params = spatial_jitter,
    point_size_mode = point_size_mode,
    alpha = alpha,
    ...
  )
}

.spatial_base <- function(
    data,
    aes_color,
    aes_size = "points_per_location",
    order_col = aes_color,
    point_size_mode = "scaled") {
  mapping <- if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    aes(
      x = .data$row,
      y = -.data$col,
      color = .data[[aes_color]],
      size = .data[[aes_size]]
    )
  } else {
    aes(
      x = .data$row,
      y = -.data$col,
      color = .data[[aes_color]]
    )
  }
  p <- ggplot(data, aes(x = .data$row, y = -.data$col)) +
    .geom_spatial_jitter(
      data = data,
      order_col = order_col,
      mapping = mapping,
      point_size_mode = point_size_mode
    ) +
    coord_fixed(ratio = 0.6) +
    theme_spatial()
  .spatial_add_size_scale(p, point_size_mode, spatial_jitter, guide = "none")
}

.tail_palette <- function(source_ct, target_ct) {
  labs <- c(
    paste0(source_ct, "_Adverse"),
    paste0(source_ct, "_Favorable"),
    paste0(target_ct, "_Adverse"),
    paste0(target_ct, "_Favorable")
  )
  cols <- c(
    "#B2182B", "#2166AC",
    "#E6550D", "#74ADD1"
  )
  if (source_ct == target_ct) {
    labs <- labs[1:2]
    cols <- cols[1:2]
  }
  stats::setNames(cols, labs)
}

.plot_tail_overview <- function(df, source_ct, target_ct, point_size_mode = "scaled") {
  ct_keep <- unique(c(source_ct, target_ct))
  tail_labs <- c(
    paste0(source_ct, "_Adverse"),
    paste0(source_ct, "_Favorable"),
    paste0(target_ct, "_Adverse"),
    paste0(target_ct, "_Favorable")
  )
  if (source_ct == target_ct) {
    tail_labs <- tail_labs[1:2]
  }

  df <- df %>%
    dplyr::mutate(
      tail_label = dplyr::if_else(
        .data$CellType %in% ct_keep & .data$prognostic_group %in% c("Adverse", "Favorable"),
        paste0(.data$CellType, "_", .data$prognostic_group),
        NA_character_
      ),
      in_analysis_tail = !is.na(.data$tail_label),
      tail_highlight = dplyr::if_else(
        .data$tail_label %in% tail_labs,
        .data$tail_label,
        "Other tail cells"
      )
    )

  df_bg <- df %>%
    dplyr::filter(.data$in_analysis_tail, !(.data$tail_label %in% tail_labs)) %>%
    .order_by_abs(., "pheno_z")

  df_hi <- df %>%
    dplyr::filter(.data$tail_label %in% tail_labs) %>%
    .order_by_abs(., "pheno_z")

  pal <- .tail_palette(source_ct, target_ct)
  pal["Other tail cells"] <- "#dddddd"

  lev <- c("Other tail cells", tail_labs)
  df_hi$tail_highlight <- factor(df_hi$tail_highlight, levels = lev)
  df_bg$tail_highlight <- factor(df_bg$tail_highlight, levels = lev)

  size_mapping <- if (.spatial_point_size_mode(point_size_mode) == "scaled") {
    aes(size = .data$points_per_location)
  } else {
    NULL
  }

  p <- ggplot(df, aes(x = .data$row, y = -.data$col)) +
    .geom_spatial_jitter(
      data = df_bg,
      order_col = "pheno_z",
      mapping = size_mapping,
      color = "#dddddd",
      alpha = 0.35,
      point_size_mode = point_size_mode
    ) +
    .geom_spatial_jitter(
      data = df_hi,
      order_col = "pheno_z",
      mapping = if (.spatial_point_size_mode(point_size_mode) == "scaled") {
        aes(color = .data$tail_highlight, size = .data$points_per_location)
      } else {
        aes(color = .data$tail_highlight)
      },
      alpha = 0.92,
      point_size_mode = point_size_mode
    ) +
    scale_color_manual(
      values = pal,
      drop = FALSE,
      name = "Prognostic tail"
    ) +
    coord_fixed(ratio = 0.6) +
    labs(
      title = "A  Prognostic tail cells (both types)",
      subtitle = paste(
        paste(ct_keep, collapse = " & "),
        "| 5th percentile Adverse / Favorable"
      )
    ) +
    theme_spatial()
  .spatial_add_size_scale(p, point_size_mode, spatial_jitter, guide = "none")
}

pal_adv <- "#B2182B"
pal_fav <- "#2166AC"
pal_bg <- "#e8e8e8"

.plot_pair_panels <- function(pair_row, cell_metrics, point_size_mode = "scaled") {
  scaled_pts <- .spatial_point_size_mode(point_size_mode) == "scaled"
  src <- pair_row$source
  tgt <- pair_row$target
  src_ct <- pair_row$source_ct
  tgt_ct <- pair_row$target_ct
  pair_title <- paste0(src, " \u2192 ", tgt)

  df <- cell_metrics %>%
    dplyr::mutate(
      pair_role = dplyr::case_when(
        .data$ct_pg == src ~ "Sender",
        .data$ct_pg == tgt ~ "Receiver",
        TRUE ~ "Other"
      ),
      pair_role = factor(.data$pair_role, levels = c("Sender", "Receiver", "Other")),
      pair_role_order = dplyr::case_when(
        .data$pair_role == "Sender" ~ abs(.data$cellchat_sender_score),
        .data$pair_role == "Receiver" ~ abs(.data$pheno_z),
        TRUE ~ 0
      ),
      integrated_hotspot = .scale01(.data$local_cooccur_score) +
        .scale01(.data$cellchat_sender_score)
    )

  p_tails <- .plot_tail_overview(df, src_ct, tgt_ct, point_size_mode = point_size_mode)

  df_other <- df %>%
    dplyr::filter(.data$pair_role == "Other") %>%
    .order_by_abs(., "pheno_z")
  df_pair <- df %>%
    dplyr::filter(.data$pair_role != "Other") %>%
    .order_by_abs(., "pair_role_order")

  p_roles <- ggplot(df, aes(x = .data$row, y = -.data$col)) +
    .geom_spatial_jitter(
      data = df_other,
      order_col = "pheno_z",
      mapping = if (scaled_pts) aes(size = .data$points_per_location) else NULL,
      color = pal_bg,
      alpha = 0.35,
      point_size_mode = point_size_mode
    ) +
    .geom_spatial_jitter(
      data = df_pair,
      order_col = "pair_role_order",
      mapping = if (scaled_pts) {
        aes(color = .data$pair_role, size = .data$points_per_location)
      } else {
        aes(color = .data$pair_role)
      },
      alpha = 0.9,
      point_size_mode = point_size_mode
    ) +
    scale_color_manual(
      values = c(Sender = pal_adv, Receiver = pal_fav, Other = pal_bg),
      drop = FALSE,
      name = "Pair role"
    ) +
    coord_fixed(ratio = 0.6) +
    labs(title = "B  Sender / receiver locations") +
    theme_spatial()
  p_roles <- .spatial_add_size_scale(p_roles, point_size_mode, spatial_jitter, guide = "none")

  p_pheno <- .spatial_base(df, "pheno_z", order_col = "pheno_z", point_size_mode = point_size_mode) +
    scale_color_gradient2(
      low = pal_fav, mid = "#F7F7F7", high = pal_adv,
      midpoint = 0,
      name = "PhenoMapR\nz-score"
    ) +
    labs(title = "C  PhenoMapR score")

  p_coloc <- .spatial_base(
    df, "local_cooccur_score", order_col = "local_cooccur_score", point_size_mode = point_size_mode
  ) +
    .spatial_coloc_score_scale_gg(name = "Local\nco-localization") +
    labs(
      title = "D  Local co-localization niche",
      subtitle = sprintf(
        "Pair mean = %.2f | nhood z = %.2f",
        pair_row$local_cooccur,
        pair_row$nhood_z
      )
    )

  df_bg <- df %>%
    dplyr::filter(.data$ct_pg != src) %>%
    .order_by_abs(., "pheno_z")
  df_sender <- df %>%
    dplyr::filter(.data$ct_pg == src) %>%
    .order_by_abs(., "cellchat_sender_score")

  p_cellchat <- ggplot(df, aes(x = .data$row, y = -.data$col)) +
    .geom_spatial_jitter(
      data = df_bg,
      order_col = "pheno_z",
      mapping = if (scaled_pts) aes(size = .data$points_per_location) else NULL,
      color = pal_bg,
      alpha = 0.25,
      point_size_mode = point_size_mode
    ) +
    .geom_spatial_jitter(
      data = df_sender,
      order_col = "cellchat_sender_score",
      mapping = if (scaled_pts) {
        aes(color = .data$cellchat_sender_score, size = .data$points_per_location)
      } else {
        aes(color = .data$cellchat_sender_score)
      },
      alpha = 0.95,
      point_size_mode = point_size_mode
    ) +
    scale_color_viridis_c(option = "magma", name = "Sender L-R\npotential") +
    coord_fixed(ratio = 0.6) +
    labs(
      title = "E  CellChat sender potential",
      subtitle = sprintf(
        "Top L-R pairs | agg prob = %.1f",
        pair_row$cellchat_prob
      )
    ) +
    theme_spatial()
  p_cellchat <- .spatial_add_size_scale(p_cellchat, point_size_mode, spatial_jitter, guide = "none")

  p_hotspot <- .spatial_base(
    df, "integrated_hotspot", order_col = "integrated_hotspot", point_size_mode = point_size_mode
  ) +
    scale_color_gradient(
      low = "#ffffcc",
      high = "#800026",
      name = "Integrated\nhotspot"
    ) +
    labs(
      title = "F  Integrated spatial + signaling hotspot",
      subtitle = "Scaled co-localization + sender L-R score"
    )

  p_nbr <- ggplot(df, aes(x = .data$row, y = -.data$col)) +
    .geom_spatial_jitter(
      data = df_bg,
      order_col = "pheno_z",
      mapping = if (scaled_pts) aes(size = .data$points_per_location) else NULL,
      color = pal_bg,
      alpha = 0.25,
      point_size_mode = point_size_mode
    ) +
    .geom_spatial_jitter(
      data = df_sender,
      order_col = "target_neighbor_frac",
      mapping = if (scaled_pts) {
        aes(color = .data$target_neighbor_frac, size = .data$points_per_location)
      } else {
        aes(color = .data$target_neighbor_frac)
      },
      alpha = 0.95,
      point_size_mode = point_size_mode
    ) +
    scale_color_viridis_c(option = "D", name = "Target\nneighbor frac", limits = c(0, 1)) +
    coord_fixed(ratio = 0.6) +
    labs(
      title = "G  Target density near senders",
      subtitle = paste0("Fraction of ", tgt, " in local neighborhood")
    ) +
    theme_spatial()
  p_nbr <- .spatial_add_size_scale(p_nbr, point_size_mode, spatial_jitter, guide = "none")

  layout <- patchwork::wrap_plots(
    p_tails, p_roles, p_pheno,
    p_coloc, p_cellchat, p_hotspot,
    p_nbr,
    ncol = 3, nrow = 3
  ) +
    patchwork::plot_annotation(
      title = paste0("Spatial evidence: ", pair_title),
      subtitle = paste0(
        "Integrated score = ", round(pair_row$integrated_score, 2),
        " | dual spatial + CellChat support | ",
        if (point_size_mode == "uniform") "uniform point size" else "size scaled by cells per spot"
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#444444")
      )
    )

  layout
}

all_metrics <- list()

for (i in seq_len(nrow(pairs_to_plot))) {
  pr <- pairs_to_plot[i, ]
  slug <- .pair_slug(pr$source, pr$target)
  message("Plotting pair ", i, "/", nrow(pairs_to_plot), ": ", pr$source, " -> ", pr$target)

  lr_rows <- lr_tbl %>%
    dplyr::filter(.data$source == pr$source, .data$target == pr$target)

  local_scores <- .spatial_cooccur_local_scores(pr$source, pr$target)
  lr_sender <- .compute_sender_lr_score(pr$source, pr$target, lr_rows)
  nbr_frac <- .compute_target_neighbor_frac(pr$source, pr$target)

  metrics <- cell_df %>%
    dplyr::mutate(
      local_cooccur_score = local_scores[match(.data$Cell, rownames(scoc_df))],
      cellchat_sender_score = lr_sender,
      target_neighbor_frac = nbr_frac,
      pair_source = pr$source,
      pair_target = pr$target
    )

  all_metrics[[slug]] <- metrics

  for (point_size_mode in c("scaled", "uniform")) {
    suffix <- .spatial_output_suffix(point_size_mode)
    scaled_pts <- .spatial_point_size_mode(point_size_mode) == "scaled"
    uni <- spatial_jitter$uniform_size

    fig <- .plot_pair_panels(pr, metrics, point_size_mode = point_size_mode)
    ggsave(
      file.path(out_dir, paste0("spatial_pair_", slug, suffix, ".png")),
      fig,
      width = 15,
      height = 14,
      dpi = 220,
      bg = "white"
    )

    p_tails_only <- .plot_tail_overview(
      metrics, pr$source_ct, pr$target_ct, point_size_mode = point_size_mode
    ) +
      patchwork::plot_annotation(
        title = paste0("Prognostic tail cells: ", pr$source_ct, " & ", pr$target_ct),
        subtitle = paste0("Pair under study: ", pr$source, " \u2192 ", pr$target),
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#444444")
        )
      )
    ggsave(
      file.path(out_dir, paste0("spatial_pair_", slug, "_tails", suffix, ".png")),
      p_tails_only,
      width = 8,
      height = 7,
      dpi = 220,
      bg = "white"
    )

    df <- metrics %>%
      dplyr::mutate(
        pair_role_order = dplyr::case_when(
          .data$ct_pg == pr$source ~ abs(.data$cellchat_sender_score),
          .data$ct_pg == pr$target ~ abs(.data$pheno_z),
          TRUE ~ 0
        )
      )
    df_bg <- df %>%
      dplyr::filter(.data$ct_pg != pr$source, .data$ct_pg != pr$target) %>%
      .order_by_abs(., "pheno_z")
    df_tgt <- df %>%
      dplyr::filter(.data$ct_pg == pr$target) %>%
      .order_by_abs(., "pheno_z")
    df_src <- df %>%
      dplyr::filter(.data$ct_pg == pr$source) %>%
      .order_by_abs(., "cellchat_sender_score")

    p_tails <- .plot_tail_overview(
      df, pr$source_ct, pr$target_ct, point_size_mode = point_size_mode
    ) +
      labs(title = "Tail cells")

    sz_bg <- if (point_size_mode == "uniform") uni else 0.4
    sz_tgt <- if (point_size_mode == "uniform") uni else 0.7
    sz_src <- if (point_size_mode == "uniform") uni else 0.9

    p1 <- ggplot(df, aes(x = .data$row, y = -.data$col)) +
      .geom_spatial_jitter(
        data = df_bg,
        order_col = "pheno_z",
        color = pal_bg,
        alpha = 0.3,
        size = sz_bg,
        point_size_mode = point_size_mode
      ) +
      .geom_spatial_jitter(
        data = df_tgt,
        order_col = "pheno_z",
        color = pal_fav,
        alpha = 0.85,
        size = sz_tgt,
        point_size_mode = point_size_mode
      ) +
      .geom_spatial_jitter(
        data = df_src,
        order_col = "cellchat_sender_score",
        mapping = aes(color = .data$cellchat_sender_score),
        alpha = 0.95,
        size = sz_src,
        point_size_mode = point_size_mode
      ) +
      scale_color_viridis_c(option = "magma", name = "Sender\nL-R") +
      coord_fixed(ratio = 0.6) +
      labs(title = "Roles + CellChat") +
      theme_spatial()
    p1 <- .spatial_add_size_scale(p1, point_size_mode, spatial_jitter, guide = "none")

    p2 <- .spatial_base(
      metrics, "local_cooccur_score", order_col = "local_cooccur_score", point_size_mode = point_size_mode
    ) +
      .spatial_coloc_score_scale_gg(name = "Co-local") +
      labs(title = "Co-localization")

    p3 <- .spatial_base(metrics, "pheno_z", order_col = "pheno_z", point_size_mode = point_size_mode) +
      scale_color_gradient2(low = pal_fav, mid = "#F7F7F7", high = pal_adv, midpoint = 0) +
      labs(title = "PhenoMapR")

    fig_compact <- patchwork::wrap_plots(p_tails, p1, p2, p3, ncol = 4) +
      patchwork::plot_annotation(
        title = paste0(pr$source, " \u2192 ", pr$target),
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
      )

    ggsave(
      file.path(out_dir, paste0("spatial_pair_", slug, "_compact", suffix, ".png")),
      fig_compact,
      width = 18,
      height = 4.8,
      dpi = 220,
      bg = "white"
    )
  }
}

metrics_rds <- list(
  pairs = pairs_to_plot,
  cells = dplyr::bind_rows(all_metrics, .id = "pair_slug")
)
saveRDS(metrics_rds, file.path(root, "inst", "data", "spatial_pair_cell_metrics.rds"))

message("Wrote spatial pair maps to ", out_dir)
message("Files:")
print(list.files(out_dir, pattern = "\\.png$", full.names = FALSE))
