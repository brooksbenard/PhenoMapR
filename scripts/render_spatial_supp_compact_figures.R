#!/usr/bin/env Rscript
# Supplemental compact panels for prognostic-cell spatial neighborhood,
# co-localization, and cell–cell communication.
#
#   1. Three-metric heatmap (nhood z, local coloc, CellChat) for all
#      Adverse/Favorable cell-type groups
#   2. Pathway × co-localized-pair heatmap (top L-R programs supporting coloc)
#   3. Compact tissue maps of the EGFR and integrin–ECM highlight axes
#   4. Compact 2×2 tissue maps of the top physical CCC + co-localization axes
#      (Fib Adv → Fib Adv COLLAGEN, Ductal Adv → Ductal Adv LAMININ,
#       Fib Adv → Ductal Adv COLLAGEN)
#
# Run from package root:
#   Rscript scripts/render_spatial_supp_compact_figures.R
#
# Requires:
#   inst/data/spatial_coloc_cellchat_pairs.rds
#   inst/data/spatial_coloc_cellchat_lr_pairs.rds
# Optional (tissue maps):
#   vignettes/HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds

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
  library(patchwork)
})

for (pkg in c("ComplexHeatmap", "circlize", "grid", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install ", pkg, " before running this script.")
  }
}
suppressPackageStartupMessages(library(ComplexHeatmap))

source(file.path(root, "scripts", "spatial_colocal_heatmap_helpers.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_plot_helpers.R"), local = TRUE)

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
lr_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_lr_pairs.rds")
if (!file.exists(pairs_rds) || !file.exists(lr_rds)) {
  stop("Missing pair RDS. Run scripts/render_spatial_colocalization_cellchat.R first.")
}

pair_df <- readRDS(pairs_rds)
lr_tbl <- readRDS(lr_rds)

out_dir <- file.path(root, "inst", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suffix_pg <- c("Adverse", "Favorable", "Other")
.parse_ct_pg <- function(lab) {
  for (s in suffix_pg) {
    suff <- paste0("_", s)
    if (nzchar(lab) && endsWith(lab, suff)) {
      ct <- substr(lab, 1L, nchar(lab) - nchar(suff))
      return(c(ct = ct, pg = s))
    }
  }
  c(ct = lab, pg = NA_character_)
}

.abbr_lab <- function(x) {
  x <- gsub("Endothelial", "Endo", x)
  x <- gsub("Fibroblast", "Fibro", x)
  x <- gsub("Macrophage", "Mac", x)
  x <- gsub("Schwann", "Schw", x)
  x <- gsub("_Adverse", " Adv", x)
  x <- gsub("_Favorable", " Fav", x)
  x
}

.pair_label <- function(source, target) {
  paste0(.abbr_lab(source), " \u2192 ", .abbr_lab(target))
}

all_labs <- sort(unique(c(pair_df$source, pair_df$target)))
meta_labs <- as.data.frame(
  do.call(rbind, lapply(all_labs, function(l) {
    p <- .parse_ct_pg(l)
    data.frame(label = l, CellType = p[["ct"]], Prognostic = p[["pg"]], stringsAsFactors = FALSE)
  })),
  stringsAsFactors = FALSE
)
meta_labs <- meta_labs[!is.na(meta_labs$Prognostic), , drop = FALSE]
labs_ord <- meta_labs$label[order(meta_labs$CellType, meta_labs$Prognostic)]

pal_pg <- c(Adverse = "#B2182B", Favorable = "#2166AC", Other = "#f7f7f7")
pal_ct <- PhenoMapR::get_celltype_palette(sort(unique(meta_labs$CellType)))
# Distinct from Adverse/Favorable RdBu so mode is not read as prognostic group.
pal_phys <- "#1B9E77"
pal_sec <- "#E69F00"
pal_other_mode <- "#7570B3"
pal_adv <- "#B2182B"
pal_fav <- "#2166AC"
dual_outline_color <- "#762A83"
col_ncells <- .spatial_load_label_ncells(root, labs_ord)

mat_nhood_z <- .build_pair_matrix(labs_ord, pair_df, "nhood_z", fill = NA_real_)
mat_nhood_q <- .build_pair_matrix(labs_ord, pair_df, "nhood_q", fill = NA_real_)
mat_local <- .build_pair_matrix(labs_ord, pair_df, "local_cooccur", fill = 0)
mat_local[!is.finite(mat_local)] <- 0
mat_cc_prob <- .build_pair_matrix(labs_ord, pair_df, "cellchat_prob", fill = 0)
mat_cc_log <- log1p(pmax(mat_cc_prob, 0))
mat_dual <- .build_pair_matrix(labs_ord, pair_df, "dual_spatial_lr", fill = FALSE)

mabs_nhood <- suppressWarnings(max(abs(as.numeric(mat_nhood_z)), na.rm = TRUE))
if (!is.finite(mabs_nhood) || mabs_nhood <= 0) mabs_nhood <- 1
col_fun_nhood <- circlize::colorRamp2(
  c(-mabs_nhood, 0, mabs_nhood),
  c("#7F312FFF", "#f7f7f7", "#005C55FF")
)
mabs_local <- suppressWarnings(max(as.numeric(mat_local), na.rm = TRUE))
if (!is.finite(mabs_local) || mabs_local <= 0) mabs_local <- 1
col_fun_local <- .spatial_coloc_score_col_fun(mabs_local)
cc_lim <- max(mat_cc_log[is.finite(mat_cc_log)], 0.01, na.rm = TRUE)
col_fun_cc <- .spatial_cellchat_prob_col_fun(cc_lim)

cell_mm_supp <- .spatial_hm_cell_mm_compact

ha_cluster <- .spatial_make_ha(
  mat_local, meta_labs, col_ncells, pal_ct, pal_pg, compact = TRUE
)
hm_ord <- .cluster_heatmap_order(
  mat_local,
  col_fun_local,
  ha_cluster$row,
  ha_cluster$col,
  distance = "euclidean"
)

ha_first <- .spatial_make_ha(
  mat_local, meta_labs, col_ncells, pal_ct, pal_pg,
  compact = TRUE,
  anno_legend_ncol = 2L,
  legend_title_size = 6L,
  show_col_anno_names = TRUE,
  show_col_legend = TRUE,
  show_ncells_axis = TRUE
)
ha_rest <- .spatial_make_ha(
  mat_local, meta_labs, col_ncells, pal_ct, pal_pg,
  compact = TRUE,
  show_col_anno_names = FALSE,
  show_col_legend = FALSE,
  show_ncells_axis = FALSE
)

.cell_fun_nhood_q <- function(mat_q) {
  function(j, i, x, y, w, h, fill) {
    qv <- mat_q[i, j]
    if (is.finite(qv) && qv < 0.05) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(col = "black", lwd = 1.1, fill = NA))
    }
  }
}

.cell_fun_local_dual <- function(mat_local, dual_flag, score_thr, dual_col) {
  function(j, i, x, y, w, h, fill) {
    val <- mat_local[i, j]
    if (is.finite(val) && val >= score_thr) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(col = "black", lwd = 1.1, fill = NA))
    }
    if (isTRUE(dual_flag[i, j])) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(col = dual_col, lwd = 1.6, fill = NA))
    }
  }
}

.cell_fun_dual <- function(dual_flag, dual_col) {
  function(j, i, x, y, w, h, fill) {
    if (isTRUE(dual_flag[i, j])) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(col = dual_col, lwd = 1.6, fill = NA))
    }
  }
}

.ha_under_colorbar <- function(name) {
  args <- list(
    ComplexHeatmap::anno_empty(border = FALSE, height = grid::unit(20, "mm")),
    show_annotation_name = FALSE
  )
  names(args)[1] <- name
  do.call(ComplexHeatmap::HeatmapAnnotation, args)
}

ht_nhood <- .build_compact_integrated_ht(
  mat_nhood_z,
  "Nhood z",
  col_fun_nhood,
  hm_ord,
  ha_row = ha_first$row,
  ha_col = ha_first$col,
  cell_fun = .cell_fun_nhood_q(mat_nhood_q),
  title_line1 = "Neighborhood enrichment",
  title_line2 = "z; black q<0.05",
  show_legend = FALSE,
  cell_mm = cell_mm_supp,
  bottom_annotation = .ha_under_colorbar("colorbar_nhood")
)
ht_coloc <- .build_compact_integrated_ht(
  mat_local,
  "Local coloc",
  col_fun_local,
  hm_ord,
  ha_row = NULL,
  ha_col = ha_rest$col,
  cell_fun = .cell_fun_local_dual(mat_local, mat_dual, 0.75, dual_outline_color),
  title_line1 = "Local co-localization",
  title_line2 = "black score\u22650.75, purple dual+",
  show_legend = FALSE,
  cell_mm = cell_mm_supp,
  bottom_annotation = .ha_under_colorbar("colorbar_coloc")
)
ht_ccc <- .build_compact_integrated_ht(
  mat_cc_log,
  "log1p(CC)",
  col_fun_cc,
  hm_ord,
  ha_row = NULL,
  ha_col = ha_rest$col,
  cell_fun = .cell_fun_dual(mat_dual, dual_outline_color),
  title_line1 = "Spatial CellChat",
  title_line2 = "log1p aggregated L-R; purple dual+",
  show_legend = FALSE,
  cell_mm = cell_mm_supp,
  bottom_annotation = .ha_under_colorbar("colorbar_ccc")
)

outline_leg <- .spatial_outline_legend(
  dual_col = dual_outline_color,
  score_label = "q<0.05 or score \u2265 0.75",
  dual_label = "Spatial + CellChat",
  legend_title_size = 6,
  point_size_mm = c(2.4, 2.8)
)

.draw_under_heatmap_colorbar <- function(anno_name, col_fun, title, at) {
  at <- at[is.finite(at)]
  if (!length(at)) at <- c(0, 1)
  lo <- min(at)
  hi <- max(at)
  if (hi <= lo) hi <- lo + 1
  labs <- format(at, trim = TRUE, scientific = FALSE, digits = 2)
  ComplexHeatmap::decorate_annotation(anno_name, {
    grid::grid.text(
      title,
      y = grid::unit(1, "npc") - grid::unit(0.6, "mm"),
      just = c("center", "top"),
      gp = grid::gpar(fontsize = 6, fontface = "bold")
    )
    pal <- matrix(col_fun(seq(lo, hi, length.out = 80)), nrow = 1)
    grid::grid.raster(
      pal,
      interpolate = TRUE,
      width = grid::unit(0.86, "npc"),
      height = grid::unit(2.2, "mm"),
      y = grid::unit(0.55, "npc")
    )
    tx <- 0.07 + 0.86 * (at - lo) / (hi - lo)
    grid::grid.text(
      labs,
      x = grid::unit(tx, "npc"),
      y = grid::unit(0.55, "npc") - grid::unit(3.6, "mm"),
      just = c("center", "top"),
      gp = grid::gpar(fontsize = 6, col = "black")
    )
  })
}

out_three <- file.path(out_dir, "spatial_supp_nhood_coloc_ccc_allgroups.png")
message("Writing ", basename(out_three), " ...")
n_lab <- ncol(mat_local)
hm_in <- n_lab * cell_mm_supp / 25.4
left_in <- 28 / 25.4
grDevices::png(
  out_three,
  width = 3 * hm_in + left_in + 1.25,
  height = hm_in + 2.15,
  units = "in",
  res = 160
)
ComplexHeatmap::draw(
  ht_nhood + ht_coloc + ht_ccc,
  column_title = .spatial_two_line_title(
    "Neighborhood, co-localization, and cell–cell communication",
    "All Adverse / Favorable cell-type groups; shared clustering from local co-localization"
  ),
  column_title_gp = .spatial_title_gp(11),
  ht_gap = grid::unit(2.5, "mm"),
  merge_legend = FALSE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  annotation_legend_list = list(outline_leg),
  padding = grid::unit(c(8, 16, 10, 6), "mm")
)
nhood_at <- unique(round(c(-mabs_nhood, 0, mabs_nhood), 0))
local_at <- unique(round(c(0, mabs_local / 2, mabs_local), 2))
cc_at <- unique(round(c(0, cc_lim / 2, cc_lim), 1))
.draw_under_heatmap_colorbar("colorbar_nhood", col_fun_nhood, "Neighborhood z", nhood_at)
.draw_under_heatmap_colorbar("colorbar_coloc", col_fun_local, "Co-localization", local_at)
.draw_under_heatmap_colorbar("colorbar_ccc", col_fun_cc, "log1p(CellChat)", cc_at)
grDevices::dev.off()

# --- Pathway heatmap: top L-R programs supporting co-localized pairs ----------
dual_pairs <- pair_df %>%
  dplyr::filter(.data$dual_spatial_lr %in% TRUE) %>%
  dplyr::mutate(
    pair_id = paste0(.data$source, "||", .data$target),
    pair_lab = .pair_label(.data$source, .data$target),
    same_type = .data$source_ct == .data$target_ct,
    pair_class = ifelse(.data$same_type, "Same type", "Cross-type")
  ) %>%
  dplyr::arrange(dplyr::desc(.data$local_cooccur), dplyr::desc(.data$cellchat_prob))

lr_dual <- lr_tbl %>%
  dplyr::inner_join(
    dual_pairs %>% dplyr::select(source, target, pair_id, pair_lab, pair_class, source_pg, target_pg),
    by = c("source", "target")
  )

mode_vec <- if ("interaction_mode" %in% names(lr_dual) && any(!is.na(lr_dual$interaction_mode))) {
  dplyr::recode(
    as.character(lr_dual$interaction_mode),
    physical = "Physical",
    secreted = "Secreted",
    other = "Other",
    .default = as.character(lr_dual$interaction_mode)
  )
} else {
  ann <- if ("annotation.x" %in% names(lr_dual)) lr_dual$annotation.x else lr_dual$annotation
  dplyr::case_when(
    ann %in% c("Cell-Cell Contact", "ECM-Receptor") ~ "Physical",
    ann == "Secreted Signaling" ~ "Secreted",
    TRUE ~ "Other"
  )
}
lr_dual$interaction_mode <- mode_vec

pw_pair <- lr_dual %>%
  dplyr::group_by(.data$pathway_name, .data$pair_lab, .data$pair_class, .data$source_pg, .data$target_pg) %>%
  dplyr::summarise(prob = sum(.data$prob, na.rm = TRUE), .groups = "drop")

pw_mode_tot <- lr_dual %>%
  dplyr::group_by(.data$pathway_name, .data$interaction_mode) %>%
  dplyr::summarise(total = sum(.data$prob, na.rm = TRUE), .groups = "drop")

.top_mode <- function(mode, n) {
  hits <- pw_mode_tot %>%
    dplyr::filter(.data$interaction_mode == !!mode) %>%
    dplyr::arrange(dplyr::desc(.data$total))
  head(hits$pathway_name, n)
}
# Physical L-R mass dominates CellChat; keep secreted programs visible so
# axes such as EGF are not dropped behind laminin/collagen.
top_pw <- unique(c(
  .top_mode("Physical", 12L),
  .top_mode("Secreted", 8L),
  "EGF"
))

pw_mode <- lr_dual %>%
  dplyr::filter(.data$pathway_name %in% top_pw) %>%
  dplyr::group_by(.data$pathway_name, .data$interaction_mode) %>%
  dplyr::summarise(prob = sum(.data$prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(.data$pathway_name) %>%
  dplyr::slice_max(.data$prob, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

pair_levels <- unique(dual_pairs$pair_lab)
pw_sub <- pw_pair %>%
  dplyr::filter(.data$pathway_name %in% top_pw, .data$pair_lab %in% pair_levels)

mat_pw <- matrix(
  0,
  nrow = length(top_pw),
  ncol = length(pair_levels),
  dimnames = list(top_pw, pair_levels)
)
for (k in seq_len(nrow(pw_sub))) {
  mat_pw[pw_sub$pathway_name[k], pw_sub$pair_lab[k]] <- pw_sub$prob[k]
}
mat_pw_log <- log1p(pmax(mat_pw, 0))

pair_meta <- dual_pairs %>%
  dplyr::distinct(.data$pair_lab, .keep_all = TRUE)
pair_meta <- pair_meta[match(colnames(mat_pw_log), pair_meta$pair_lab), , drop = FALSE]

pw_mode <- pw_mode[match(rownames(mat_pw_log), pw_mode$pathway_name), , drop = FALSE]
pw_mode$interaction_mode[is.na(pw_mode$interaction_mode)] <- "Other"

pal_mode <- c(Physical = pal_phys, Secreted = pal_sec, Other = pal_other_mode)
pw_lim <- max(mat_pw_log, 0.01, na.rm = TRUE)
col_fun_pw <- .spatial_cellchat_prob_col_fun(pw_lim)

ha_pw_col <- ComplexHeatmap::HeatmapAnnotation(
  `Local coloc` = pair_meta$local_cooccur,
  `Sender type` = pair_meta$source_ct,
  `Receiver type` = pair_meta$target_ct,
  Sender = pair_meta$source_pg,
  Receiver = pair_meta$target_pg,
  col = list(
    `Local coloc` = .spatial_coloc_score_col_fun(
      max(pair_meta$local_cooccur, 0.01, na.rm = TRUE)
    ),
    `Sender type` = pal_ct,
    `Receiver type` = pal_ct,
    Sender = pal_pg[c("Adverse", "Favorable")],
    Receiver = pal_pg[c("Adverse", "Favorable")]
  ),
  show_legend = c(
    `Local coloc` = TRUE,
    `Sender type` = TRUE,
    `Receiver type` = FALSE,
    Sender = TRUE,
    Receiver = FALSE
  ),
  annotation_name_gp = grid::gpar(fontsize = 7),
  annotation_legend_param = list(
    `Local coloc` = list(
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 7)
    ),
    `Sender type` = list(
      title = "Cell type",
      ncol = 2L,
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 7)
    ),
    Sender = list(
      title = "Prognostic group",
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 7)
    )
  ),
  simple_anno_size = grid::unit(2.4, "mm"),
  gap = grid::unit(0.6, "mm")
)
ha_pw_row <- ComplexHeatmap::rowAnnotation(
  Mode = pw_mode$interaction_mode,
  col = list(Mode = pal_mode),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Mode = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7))
  ),
  simple_anno_size = grid::unit(2.6, "mm")
)

col_split <- factor(pair_meta$pair_class, levels = c("Cross-type", "Same type"))
row_split <- factor(pw_mode$interaction_mode, levels = c("Physical", "Secreted", "Other"))

ht_pw <- ComplexHeatmap::Heatmap(
  mat_pw_log,
  name = "log1p(prob)",
  col = col_fun_pw,
  na_col = "#f7f7f7",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = .spatial_cor_dist,
  clustering_distance_columns = .spatial_cor_dist,
  row_split = row_split,
  column_split = col_split,
  row_gap = grid::unit(1.4, "mm"),
  column_gap = grid::unit(1.8, "mm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 8),
  top_annotation = ha_pw_col,
  left_annotation = ha_pw_row,
  border = TRUE,
  heatmap_legend_param = list(
    title = "log1p(CellChat)",
    title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 7)
  ),
  column_title = NULL,
  row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  width = grid::unit(ncol(mat_pw_log) * 4.2, "mm"),
  height = grid::unit(nrow(mat_pw_log) * 4.6, "mm")
)

out_pw <- file.path(out_dir, "spatial_supp_coloc_pathway_heatmap.png")
message("Writing ", basename(out_pw), " ...")
grDevices::png(out_pw, width = 11.2, height = 6.6, units = "in", res = 160)
ComplexHeatmap::draw(
  ht_pw,
  column_title = .spatial_two_line_title(
    "Top pathways supporting co-localization of prognostic cell types",
    "Dual-positive pairs; top physical and secreted programs (EGF retained); fill = log1p summed L-R probability"
  ),
  column_title_gp = .spatial_title_gp(11),
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = grid::unit(c(8, 4, 4, 4), "mm")
)
grDevices::dev.off()

# --- Compact tissue maps for highlight axes -----------------------------------
# Match CytoSPACE overview maps (e.g. spatial_cytospace_celltype_scaled.png).
.cyto_coord <- function() ggplot2::coord_fixed(ratio = 0.6)

.compact_map_theme <- function(base_size = 6) {
  .spatial_map_theme(base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.height = grid::unit(0.18, "cm"),
      legend.key.width = grid::unit(0.5, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-8, 0, 0, 0),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.margin = ggplot2::margin(0, 1, 0, 1)
    )
}

.int_tick_labels <- function(x) {
  sprintf("%d", as.integer(round(x)))
}

.whole_breaks <- function(vals, min_hi = 1L) {
  hi <- suppressWarnings(max(vals, na.rm = TRUE))
  if (!is.finite(hi) || hi <= 0) hi <- min_hi
  hi <- as.integer(max(min_hi, ceiling(hi)))
  list(limits = c(0, hi), breaks = seq(0, hi, by = 1))
}

.compact_colorbar <- function(name) {
  ggplot2::guides(
    colour = ggplot2::guide_colorbar(
      title = name,
      barwidth = grid::unit(0.95, "cm"),
      barheight = grid::unit(0.14, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  )
}

.build_axis_plots <- function(
    cell_df, expr, score_col, jitter,
    source_lab, target_lab,
    highlight_genes_lig, highlight_genes_rec,
    lig_title, rec_title
) {
  .mean_genes <- function(genes) {
    genes <- intersect(genes, rownames(expr))
    if (!length(genes)) return(rep(0, nrow(cell_df)))
    mat <- as.matrix(expr[genes, , drop = FALSE])
    if (nrow(mat) == 1L) as.numeric(mat[1, ]) else as.numeric(Matrix::colMeans(mat))
  }

  df <- cell_df
  df$lig_score <- .mean_genes(highlight_genes_lig)
  df$rec_score <- .mean_genes(highlight_genes_rec)
  df$role <- dplyr::case_when(
    identical(source_lab, target_lab) & df$ct_pg == source_lab ~ "Adverse ductal",
    df$ct_pg == source_lab ~ "Sender",
    df$ct_pg == target_lab ~ "Receiver",
    grepl("_Adverse$", df$ct_pg) ~ "Other Adverse",
    grepl("_Favorable$", df$ct_pg) ~ "Favorable",
    TRUE ~ "Other"
  )
  bg <- df %>% dplyr::filter(.data$role %in% c("Other", "Favorable", "Other Adverse"))
  fg <- df %>% dplyr::filter(!.data$role %in% c("Other", "Favorable", "Other Adverse"))
  focus_cells <- if (identical(source_lab, target_lab)) {
    df %>% dplyr::filter(.data$ct_pg == source_lab)
  } else {
    df %>% dplyr::filter(.data$ct_pg %in% c(source_lab, target_lab))
  }
  sender_cells <- df %>% dplyr::filter(.data$ct_pg == source_lab)
  receiver_cells <- df %>% dplyr::filter(.data$ct_pg == target_lab)

  p_roles <- ggplot() +
    .spatial_geom_jitter(
      bg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#dddddd", alpha = 0.22
    ) +
    .spatial_geom_jitter(
      fg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$role), alpha = 0.9
    ) +
    scale_color_manual(
      values = c(
        "Adverse ductal" = pal_adv,
        "Sender" = "#E69F00",
        "Receiver" = pal_adv,
        "Other Adverse" = "#f4a582",
        "Favorable" = pal_fav,
        "Other" = "#dddddd"
      ),
      breaks = unique(fg$role)
    ) +
    .cyto_coord() +
    labs(title = "Sender / receiver", color = NULL) +
    .compact_map_theme(8) +
    theme(legend.key.width = grid::unit(0.35, "cm"))

  p_pheno <- ggplot() +
    .spatial_geom_jitter(
      df %>% dplyr::filter(!.data$ct_pg %in% unique(c(source_lab, target_lab))),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.18
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
    .cyto_coord() +
    labs(title = "PhenoMapR (axis cells)") +
    .compact_map_theme(8) +
    .compact_colorbar("PhenoMapR")

  p_lig <- ggplot() +
    .spatial_geom_jitter(
      df %>% dplyr::filter(.data$ct_pg != source_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.12
    ) +
    .spatial_geom_jitter(
      sender_cells, order_col = "lig_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$lig_score), alpha = 0.9
    ) +
    scale_color_gradientn(colours = .spatial_cellchat_prob_palette_colors(), name = "Ligand") +
    .cyto_coord() +
    labs(title = lig_title) +
    .compact_map_theme(8) +
    .compact_colorbar("Ligand")

  p_rec <- ggplot() +
    .spatial_geom_jitter(
      df %>% dplyr::filter(.data$ct_pg != target_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.12
    ) +
    .spatial_geom_jitter(
      receiver_cells, order_col = "rec_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$rec_score), alpha = 0.9
    ) +
    scale_color_gradientn(colours = .spatial_cellchat_prob_palette_colors(), name = "Receptor") +
    .cyto_coord() +
    labs(title = rec_title) +
    .compact_map_theme(8) +
    .compact_colorbar("Receptor")

  list(roles = p_roles, pheno = p_pheno, lig = p_lig, rec = p_rec)
}

.mean_expr_genes <- function(expr, n_cells, genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(0, n_cells))
  mat <- as.matrix(expr[genes, , drop = FALSE])
  if (nrow(mat) == 1L) as.numeric(mat[1, ]) else as.numeric(Matrix::colMeans(mat))
}

.local_cooccur_cell_scores <- function(neighbor_labels, adj, degrees, source_lab, target_lab) {
  binary <- vapply(
    neighbor_labels,
    function(nl) {
      if (!length(nl)) return(0)
      as.numeric(any(nl == source_lab) & any(nl == target_lab))
    },
    numeric(1)
  )
  s_norm <- matrix(binary / pmax(as.numeric(degrees), 1), ncol = 1L)
  as.numeric((adj %*% s_norm) + s_norm)
}

.role_color <- function(lab) {
  ct <- .parse_ct_pg(lab)[["ct"]]
  col <- unname(pal_ct[ct])
  if (!length(col) || is.na(col) || !nzchar(col)) pal_adv else col
}

.build_top3_axis_panels <- function(
    cell_df, expr, jitter, coloc_score,
    source_lab, target_lab,
    lig_genes, rec_genes,
    lig_title, rec_title
) {
  df <- cell_df
  df$lig_score <- .mean_expr_genes(expr, nrow(df), lig_genes)
  df$rec_score <- .mean_expr_genes(expr, nrow(df), rec_genes)
  df$coloc <- pmin(pmax(as.numeric(coloc_score), 0), 1)
  same <- identical(source_lab, target_lab)
  src_lab <- .abbr_lab(source_lab)
  tgt_lab <- .abbr_lab(target_lab)
  sender_name <- if (same) paste0(src_lab, " (autocrine)") else paste0("Sender (", src_lab, ")")
  receiver_name <- if (same) sender_name else paste0("Receiver (", tgt_lab, ")")
  df$role <- dplyr::case_when(
    df$ct_pg == source_lab ~ sender_name,
    !same & df$ct_pg == target_lab ~ receiver_name,
    TRUE ~ "Other"
  )
  bg <- df %>% dplyr::filter(.data$role == "Other")
  fg <- df %>% dplyr::filter(.data$role != "Other")
  coloc_fg <- fg %>%
    dplyr::filter(is.finite(.data$coloc), .data$coloc > 0)
  jitter_coloc <- jitter
  jitter_coloc$uniform_size <- max(jitter$uniform_size, 0.42)
  sender_cells <- df %>% dplyr::filter(.data$ct_pg == source_lab)
  receiver_cells <- df %>% dplyr::filter(.data$ct_pg == target_lab)
  role_vals <- if (same) {
    stats::setNames(pal_adv, sender_name)
  } else {
    vals <- c(
      stats::setNames(.role_color(source_lab), sender_name),
      stats::setNames(.role_color(target_lab), receiver_name)
    )
    vals[!duplicated(names(vals))]
  }
  lig_sc <- .whole_breaks(sender_cells$lig_score)
  rec_sc <- .whole_breaks(receiver_cells$rec_score)
  map_theme <- .compact_map_theme(5.5) +
    ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-8, 0, -1, 0),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 5.5)
    )

  p_coloc <- ggplot() +
    .spatial_geom_jitter(
      bg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.18
    ) +
    .spatial_geom_jitter(
      coloc_fg, order_col = "coloc", jitter_params = jitter_coloc, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$coloc), alpha = 1
    ) +
    .spatial_coloc_score_scale_gg(
      name = "Coloc",
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = .int_tick_labels
    ) +
    .cyto_coord() +
    labs(title = "Co-localization") +
    map_theme +
    .compact_colorbar("Coloc")

  p_roles <- ggplot() +
    .spatial_geom_jitter(
      bg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#dddddd", alpha = 0.22
    ) +
    .spatial_geom_jitter(
      fg, order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$role), alpha = 0.9
    ) +
    scale_color_manual(values = role_vals, breaks = names(role_vals), name = NULL) +
    .cyto_coord() +
    labs(title = "Sender / receiver") +
    map_theme +
    theme(
      legend.key.width = grid::unit(0.22, "cm"),
      legend.key.height = grid::unit(0.22, "cm"),
      legend.text = element_text(size = 5.5)
    )

  p_lig <- ggplot() +
    .spatial_geom_jitter(
      df %>% dplyr::filter(.data$ct_pg != source_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.12
    ) +
    .spatial_geom_jitter(
      sender_cells, order_col = "lig_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$lig_score), alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = .spatial_pubu_palette_colors(),
      limits = lig_sc$limits,
      breaks = lig_sc$breaks,
      labels = .int_tick_labels,
      name = "Ligand"
    ) +
    .cyto_coord() +
    labs(title = lig_title) +
    map_theme +
    .compact_colorbar("Ligand")

  p_rec <- ggplot() +
    .spatial_geom_jitter(
      df %>% dplyr::filter(.data$ct_pg != target_lab),
      order_col = NULL, jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col), color = "#eeeeee", alpha = 0.12
    ) +
    .spatial_geom_jitter(
      receiver_cells, order_col = "rec_score", jitter_params = jitter, point_size_mode = "uniform",
      mapping = aes(x = .data$row, y = -.data$col, color = .data$rec_score), alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = .spatial_cellchat_prob_palette_colors(),
      limits = rec_sc$limits,
      breaks = rec_sc$breaks,
      labels = .int_tick_labels,
      name = "Receptor"
    ) +
    .cyto_coord() +
    labs(title = rec_title) +
    map_theme +
    .compact_colorbar("Receptor")

  list(coloc = p_coloc, roles = p_roles, lig = p_lig, rec = p_rec)
}

vdir <- file.path(root, "vignettes")
rds_cyto <- file.path(vdir, "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

skip_maps <- identical(tolower(Sys.getenv("SPATIAL_SUPP_SKIP_MAPS", unset = "0")), "1")
if (
  !skip_maps &&
    file.exists(rds_cyto) &&
    requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE) &&
    requireNamespace("RANN", quietly = TRUE)
) {
  message("Loading CytoSPACE object for compact spatial maps ...")
  seurat <- .spatial_load_cytospace_seurat(rds_cyto, gd_id = gd_id_cyto)
  cyto <- .spatial_cytospace_cell_df(seurat)
  cell_df <- cyto$cell_df
  score_col <- cyto$score_col
  jitter <- .spatial_jitter_params(cell_df, point_range = c(0.18, 0.7), uniform_size = 0.28)

  expr_counts <- tryCatch(
    Seurat::GetAssayData(seurat, assay = "Spatial", layer = "counts"),
    error = function(e) Seurat::GetAssayData(seurat, assay = "Spatial", slot = "counts")
  )
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
  rm(seurat, expr_counts, cs)
  invisible(gc())

  xy <- as.matrix(cell_df[, c("row", "col"), drop = FALSE])
  rownames(xy) <- as.character(seq_len(nrow(xy)))
  if (any(duplicated(xy) | duplicated(xy, fromLast = TRUE))) {
    rng <- max(diff(range(xy[, 1], na.rm = TRUE)), diff(range(xy[, 2], na.rm = TRUE)), na.rm = TRUE)
    eps <- if (is.finite(rng) && rng > 0) rng * 1e-5 else 1e-6
    set.seed(4L)
    xy <- xy + cbind(
      stats::runif(nrow(xy), 0, eps),
      stats::runif(nrow(xy), 0, eps)
    )
  }
  k_nn <- min(20L, max(5L, nrow(xy) - 1L))
  rad_nn <- max(
    3,
    0.04 * max(diff(range(xy[, 1], na.rm = TRUE)), diff(range(xy[, 2], na.rm = TRUE)), na.rm = TRUE)
  )
  nn_idx <- RANN::nn2(
    data = xy,
    query = xy,
    searchtype = "radius",
    radius = rad_nn,
    k = k_nn
  )$nn.idx
  labels <- as.character(cell_df$ct_pg)
  neighbor_labels <- vector("list", nrow(cell_df))
  for (i in seq_len(nrow(cell_df))) {
    ni <- nn_idx[i, ]
    ni <- ni[ni > 0L & ni != i]
    neighbor_labels[[i]] <- if (length(ni)) labels[ni] else character(0)
  }
  neighbors_local <- Seurat::FindNeighbors(xy, k.param = k_nn, verbose = FALSE)
  adj_local <- neighbors_local$nn
  degrees_local <- Matrix::colSums(adj_local) + 1

  top3_axes <- list(
    list(
      source = "Fibroblast_Adverse",
      target = "Fibroblast_Adverse",
      pathway = "COLLAGEN",
      lig = c("COL1A1", "COL1A2", "COL6A3"),
      rec = c("ITGA1", "ITGA2", "ITGB1"),
      lig_title = "Ligands (COL1A1/A2, COL6A3)",
      rec_title = "Receptors (ITGA1/A2/B1)",
      title = "Fibro Adv \u2192 Fibro Adv \u00b7 COLLAGEN",
      file = "spatial_axis_fibro_collagen_spatial_compact.png"
    ),
    list(
      source = "Ductal_Adverse",
      target = "Ductal_Adverse",
      pathway = "LAMININ",
      lig = c("LAMC2", "LAMB3", "LAMA3"),
      rec = c("ITGA2", "ITGA3", "ITGA6", "ITGB1"),
      lig_title = "Ligands (LAMC2/B3/A3)",
      rec_title = "Receptors (ITGA2/A3/A6/B1)",
      title = "Ductal Adv \u2192 Ductal Adv \u00b7 LAMININ",
      file = "spatial_axis_ductal_laminin_spatial_compact.png"
    ),
    list(
      source = "Fibroblast_Adverse",
      target = "Ductal_Adverse",
      pathway = "COLLAGEN",
      lig = c("COL1A1", "COL1A2", "COL6A3"),
      rec = c("ITGA2", "ITGA3", "ITGB1"),
      lig_title = "Ligands (COL1A1/A2, COL6A3)",
      rec_title = "Receptors (ITGA2/A3/B1)",
      title = "Fibro Adv \u2192 Ductal Adv \u00b7 COLLAGEN",
      file = "spatial_axis_fibro_ductal_collagen_spatial_compact.png"
    )
  )
  for (ax in top3_axes) {
    message("Writing ", ax$file, " ...")
    panels <- .build_top3_axis_panels(
      cell_df, expr, jitter,
      coloc_score = .local_cooccur_cell_scores(
        neighbor_labels, adj_local, degrees_local, ax$source, ax$target
      ),
      source_lab = ax$source,
      target_lab = ax$target,
      lig_genes = ax$lig,
      rec_genes = ax$rec,
      lig_title = ax$lig_title,
      rec_title = ax$rec_title
    )
    spat_ax <- (panels$coloc | panels$roles) / (panels$lig | panels$rec) +
      patchwork::plot_annotation(
        title = ax$title,
        theme = theme(
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1)
        )
      )
    ggsave(
      file.path(out_dir, ax$file),
      spat_ax, width = 4.25, height = 4.2, dpi = 240, bg = "white"
    )
  }

  plots_a <- .build_axis_plots(
    cell_df, expr, score_col, jitter,
    source_lab = "Ductal_Adverse",
    target_lab = "Ductal_Adverse",
    highlight_genes_lig = c("TGFA", "EREG", "AREG"),
    highlight_genes_rec = c("EGFR", "ERBB2"),
    lig_title = "Ligands (TGFA/EREG/AREG)",
    rec_title = "Receptors (EGFR/ERBB2)"
  )
  plots_b <- .build_axis_plots(
    cell_df, expr, score_col, jitter,
    source_lab = "Fibroblast_Adverse",
    target_lab = "Ductal_Adverse",
    highlight_genes_lig = c("COL1A1", "COL1A2", "FN1", "POSTN"),
    highlight_genes_rec = c("ITGAV", "ITGB1", "ITGA2"),
    lig_title = "Ligands (COL1A1/2, FN1, POSTN)",
    rec_title = "Receptors (ITGAV/B1/A2)"
  )

  spat_a <- (plots_a$roles | plots_a$pheno) / (plots_a$lig | plots_a$rec) +
    patchwork::plot_annotation(
      title = "Axis A — Ductal Adv EGFR autocrine niche",
      theme = theme(plot.title = element_text(face = "bold", size = 11))
    )
  spat_b <- (plots_b$roles | plots_b$pheno) / (plots_b$lig | plots_b$rec) +
    patchwork::plot_annotation(
      title = "Axis B — Fibro Adv \u2192 Ductal Adv integrin–ECM",
      theme = theme(plot.title = element_text(face = "bold", size = 11))
    )

  ggsave(
    file.path(out_dir, "spatial_axis_ductal_egfr_spatial_compact.png"),
    spat_a, width = 6.6, height = 6.0, dpi = 200, bg = "white"
  )
  ggsave(
    file.path(out_dir, "spatial_axis_fibro_ductal_ecm_spatial_compact.png"),
    spat_b, width = 6.6, height = 6.0, dpi = 200, bg = "white"
  )

  spat_a_row <- plots_a$roles | plots_a$pheno | plots_a$lig | plots_a$rec
  spat_b_row <- plots_b$roles | plots_b$pheno | plots_b$lig | plots_b$rec
  spat_combined <- (spat_a_row / spat_b_row) +
    patchwork::plot_annotation(
      title = "Compact spatial co-localization of prognostic highlight axes",
      subtitle = "Top: Axis A Ductal Adv EGFR  |  Bottom: Axis B Fibro Adv \u2192 Ductal Adv ECM",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "#444444")
      )
    )
  ggsave(
    file.path(out_dir, "spatial_supp_axis_spatial_compact.png"),
    spat_combined, width = 10.6, height = 6.4, dpi = 200, bg = "white"
  )
} else if (skip_maps) {
  message("Skipping compact spatial maps (SPATIAL_SUPP_SKIP_MAPS=1).")
} else {
  message("Skipping compact spatial maps (missing CytoSPACE RDS, Seurat/Matrix, or RANN).")
}

message("Done. Supplemental compact figures written to ", out_dir)
message("  - spatial_supp_nhood_coloc_ccc_allgroups.png")
message("  - spatial_supp_coloc_pathway_heatmap.png")
message("  - spatial_axis_fibro_collagen_spatial_compact.png")
message("  - spatial_axis_ductal_laminin_spatial_compact.png")
message("  - spatial_axis_fibro_ductal_collagen_spatial_compact.png")
message("  - spatial_axis_ductal_egfr_spatial_compact.png")
message("  - spatial_axis_fibro_ductal_ecm_spatial_compact.png")
message("  - spatial_supp_axis_spatial_compact.png")
