#!/usr/bin/env Rscript
# Regenerate integrated Spearman + spatial CellChat heatmaps from cached pair table.
# Run after render_spatial_colocalization_cellchat.R (or when only figures need refresh):
#   Rscript scripts/render_spatial_integrated_heatmaps.R

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
})

for (pkg in c("ComplexHeatmap", "circlize")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install ", pkg, " before running this script.")
  }
}
suppressPackageStartupMessages(library(ComplexHeatmap))

source(file.path(root, "scripts", "spatial_colocal_heatmap_helpers.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)

pairs_rds <- file.path(root, "inst", "data", "spatial_coloc_cellchat_pairs.rds")
if (!file.exists(pairs_rds)) {
  stop("Missing ", pairs_rds, ". Run scripts/render_spatial_colocalization_cellchat.R first.")
}
pair_df <- readRDS(pairs_rds)

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
col_ncells <- .spatial_load_label_ncells(root, labs_ord)

mat_cor_pheno_neighbor <- .build_pair_matrix(
  labs_ord, pair_df, "pheno_neighbor_rho", fill = NA_real_
)
mat_cor_pheno_cooccur <- .build_pair_matrix(
  labs_ord, pair_df, "pheno_cooccur_rho", fill = NA_real_
)

out_dir <- file.path(root, "inst", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cor_outline_threshold <- 0.5
dual_outline_color <- "#762A83"

mat_dual <- .build_pair_matrix(labs_ord, pair_df, "dual_spatial_lr", fill = FALSE)
mat_cc_prob <- .build_pair_matrix(labs_ord, pair_df, "cellchat_prob", fill = 0)
mat_cc_log <- log1p(pmax(mat_cc_prob, 0))
mat_integrated <- .build_pair_matrix(labs_ord, pair_df, "integrated_score", fill = 0)

col_fun_cor <- .spatial_spearman_col_fun()
cc_lim <- max(mat_cc_log[is.finite(mat_cc_log)], 0.01, na.rm = TRUE)
col_fun_cc <- .spatial_cellchat_prob_col_fun(cc_lim)
int_lim <- max(abs(mat_integrated[is.finite(mat_integrated)]), 0.01, na.rm = TRUE)
col_fun_int <- circlize::colorRamp2(
  c(-int_lim, 0, int_lim),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

ha_cor <- .spatial_make_ha(
  mat_cor_pheno_neighbor, meta_labs, col_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = 9L
)
hm_ord <- .cluster_heatmap_order(
  mat_cor_pheno_neighbor,
  col_fun_cor,
  ha_cor$row,
  ha_cor$col
)

.cell_fun_rho_and_dual <- function(mat_rho, mat_dual_flag, rho_thr, dual_col) {
  function(j, i, x, y, w, h, fill) {
    val <- mat_rho[i, j]
    if (is.finite(val) && abs(val) > rho_thr) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
      )
    }
    if (isTRUE(mat_dual_flag[i, j])) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = dual_col, lwd = 2.2, fill = NA)
      )
    }
  }
}

cell_fun_integrated <- .cell_fun_rho_and_dual(
  mat_cor_pheno_neighbor,
  mat_dual,
  cor_outline_threshold,
  dual_outline_color
)
cell_fun_cooccur <- .cell_fun_rho_and_dual(
  mat_cor_pheno_cooccur,
  mat_dual,
  cor_outline_threshold,
  dual_outline_color
)

wh_full <- .spatial_hm_wh(mat_cor_pheno_neighbor)
png_sz <- .spatial_png_inches(mat_cor_pheno_neighbor, pad_h = 3.8)

ht_neighbor_outline <- ComplexHeatmap::Heatmap(
  mat_cor_pheno_neighbor,
  name = "Spearman",
  col = col_fun_cor,
  width = wh_full$width,
  height = wh_full$height,
  na_col = "#eeeeee",
  row_order = hm_ord$row,
  column_order = hm_ord$col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_cor$row,
  top_annotation = ha_cor$col,
  border = TRUE,
  heatmap_legend_param = list(title = "Spearman rho"),
  cell_fun = cell_fun_integrated
)

.spatial_draw_ht_png(
  ht_neighbor_outline,
  file.path(
    out_dir,
    "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined_integrated.png"
  ),
  "PhenoMapR vs neighbor PhenoMapR (clustered)",
  paste0(
    "|rho| > ", cor_outline_threshold, " black; spatial + CellChat purple"
  ),
  width_in = png_sz["width"],
  height_in = png_sz["height"],
  legend_mode = "right_1col",
  annotation_legend_list = list(
    .spatial_outline_legend(
      dual_col = dual_outline_color,
      score_label = paste0("|rho| > ", cor_outline_threshold),
      dual_label = "Spatial + CellChat"
    )
  )
)

ht_cooccur_outline <- ComplexHeatmap::Heatmap(
  mat_cor_pheno_cooccur,
  name = "Spearman",
  col = col_fun_cor,
  width = wh_full$width,
  height = wh_full$height,
  na_col = "#eeeeee",
  row_order = hm_ord$row,
  column_order = hm_ord$col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_cor$row,
  top_annotation = ha_cor$col,
  border = TRUE,
  heatmap_legend_param = list(title = "Spearman rho"),
  cell_fun = cell_fun_cooccur
)

.spatial_draw_ht_png(
  ht_cooccur_outline,
  file.path(
    out_dir,
    "spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined_integrated.png"
  ),
  "PhenoMapR vs local co-localization (clustered)",
  paste0(
    "|rho| > ", cor_outline_threshold, " black; spatial + CellChat purple"
  ),
  width_in = png_sz["width"],
  height_in = png_sz["height"],
  legend_mode = "right_1col",
  annotation_legend_list = list(
    .spatial_outline_legend(
      dual_col = dual_outline_color,
      score_label = paste0("|rho| > ", cor_outline_threshold),
      dual_label = "Spatial + CellChat"
    )
  )
)

ht_cc <- ComplexHeatmap::Heatmap(
  mat_cc_log,
  name = "log1p(CellChat)",
  col = col_fun_cc,
  width = wh_full$width,
  height = wh_full$height,
  na_col = "#eeeeee",
  row_order = hm_ord$row,
  column_order = hm_ord$col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_cor$row,
  top_annotation = ha_cor$col,
  border = TRUE
)

.spatial_draw_ht_png(
  ht_cc,
  file.path(out_dir, "spatial_coloc_cellchat_prob_heatmap_clustered.png"),
  "Spatial CellChat L-R probability (log1p)",
  "Same row/column order as neighbor PhenoMapR rho",
  width_in = png_sz["width"],
  height_in = png_sz["height"],
  legend_mode = "right_1col"
)

# --- Compact 2x2 integrated panel (fits vignette width) ---
ha_compact <- .spatial_make_ha(
  mat_cor_pheno_neighbor, meta_labs, col_ncells, pal_ct, pal_pg, compact = TRUE
)

ht_p1 <- .build_compact_integrated_ht(
  mat_cor_pheno_neighbor,
  "Spearman",
  col_fun_cor,
  hm_ord,
  ha_row = ha_compact$row,
  ha_col = ha_compact$col,
  cell_fun = cell_fun_integrated,
  title_line1 = "Neighbor PhenoMapR rho",
  title_line2 = paste0("|rho|> ", cor_outline_threshold, " black; dual+ purple"),
  legend_param = list(title = "Spearman rho")
)
ht_p2 <- .build_compact_integrated_ht(
  mat_cc_log,
  "log1p(CC)",
  col_fun_cc,
  hm_ord,
  ha_row = ha_compact$row,
  ha_col = ha_compact$col,
  title_line1 = "Spatial CellChat prob",
  title_line2 = "log1p aggregated L-R",
  show_legend = TRUE
)
ht_p3 <- .build_compact_integrated_ht(
  mat_integrated,
  "Integrated",
  col_fun_int,
  hm_ord,
  ha_row = ha_compact$row,
  ha_col = ha_compact$col,
  title_line1 = "Integrated score",
  title_line2 = "Spatial + PhenoMapR + CellChat",
  show_legend = TRUE
)
ht_p4 <- .build_compact_integrated_ht(
  mat_dual * 1,
  "Dual+",
  c(`0` = "#f0f0f0", `1` = dual_outline_color),
  hm_ord,
  ha_row = ha_compact$row,
  ha_col = ha_compact$col,
  title_line1 = "Dual-positive pairs",
  title_line2 = "Spatial co-localization + CellChat",
  legend_param = list(
    title = "Spatial +\nCellChat",
    at = c(0, 1),
    labels = c("No", "Yes")
  )
)

out_four <- file.path(out_dir, "spatial_coloc_integrated_four_panel.png")
png_four <- .spatial_png_inches(
  mat_cor_pheno_neighbor,
  pad_w = 2.5,
  pad_h = 1.8,
  cell_mm = .spatial_hm_cell_mm_compact
)
message("Writing ", basename(out_four), " ...")
grDevices::png(
  out_four,
  width = min(8.5, png_four["width"] + 0.5),
  height = png_four["height"] * 4.1,
  units = "in",
  res = 150
)
ComplexHeatmap::draw(
  ht_p1 %v% ht_p2 %v% ht_p3 %v% ht_p4,
  column_title = .spatial_two_line_title(
    "Integrated spatial co-localization, PhenoMapR, and spatial CellChat",
    "Prognostic sender \u2192 receiver pairs"
  ),
  column_title_gp = .spatial_title_gp(11),
  ht_gap = grid::unit(2, "mm"),
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  padding = grid::unit(c(12, 22, 2, 12), "mm")
)
grDevices::dev.off()

pair_prognostic <- pair_df %>% dplyr::filter(.data$prognostic_pair)
tier_df <- pair_prognostic %>%
  dplyr::mutate(
    evidence_tier = dplyr::case_when(
      .data$dual_spatial_lr ~ "Spatial + CellChat",
      .data$spatial_coloc & !.data$cellchat_sig ~ "Spatial only",
      !.data$spatial_coloc & .data$cellchat_sig ~ "CellChat only",
      TRUE ~ "Neither"
    )
  ) %>%
  dplyr::count(.data$evidence_tier) %>%
  dplyr::mutate(
    evidence_tier = factor(
      .data$evidence_tier,
      levels = c("Spatial + CellChat", "Spatial only", "CellChat only", "Neither")
    )
  )
p_ev <- ggplot(tier_df, aes(x = .data$evidence_tier, y = .data$n, fill = .data$evidence_tier)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = .data$n), vjust = -0.35, size = 3.8) +
  scale_fill_manual(values = c(
    "Spatial + CellChat" = dual_outline_color,
    "Spatial only" = "#5AAE61",
    "CellChat only" = "#9970AB",
    "Neither" = "#dddddd"
  )) +
  labs(
    title = "Integrated evidence among prognostic sender \u2192 receiver pairs",
    subtitle = paste0(
      sum(pair_prognostic$dual_spatial_lr), " / ",
      nrow(pair_prognostic),
      " pairs with spatial co-localization and spatial CellChat L-R support"
    ),
    x = NULL,
    y = "Pair count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 22, hjust = 1)
  )
out_ev <- file.path(out_dir, "spatial_coloc_integration_evidence.png")
message("Writing ", basename(out_ev), " ...")
ggsave(out_ev, p_ev, width = 8.5, height = 5.5, dpi = 150, bg = "white")

# --- Adverse / Favorable only: neighborhood enrichment + coloc × CellChat ---
# Pair table is already filtered to Adverse/Favorable tails (no *_Other).
af_labs <- labs_ord
af_meta <- meta_labs[meta_labs$label %in% af_labs, , drop = FALSE]
af_ncells <- col_ncells[af_labs]

mat_nhood_z <- .build_pair_matrix(af_labs, pair_df, "nhood_z", fill = NA_real_)
mat_nhood_q <- .build_pair_matrix(af_labs, pair_df, "nhood_q", fill = NA_real_)
mat_local_af <- .build_pair_matrix(af_labs, pair_df, "local_cooccur", fill = 0)
mat_local_af[!is.finite(mat_local_af)] <- 0
mat_dual_af <- .build_pair_matrix(af_labs, pair_df, "dual_spatial_lr", fill = FALSE)
mat_cc_sig_af <- .build_pair_matrix(af_labs, pair_df, "cellchat_sig", fill = FALSE)

mabs_af <- suppressWarnings(max(abs(as.numeric(mat_nhood_z)), na.rm = TRUE))
if (!is.finite(mabs_af) || mabs_af <= 0) mabs_af <- 1
col_fun_nhood_af <- circlize::colorRamp2(
  c(-mabs_af, 0, mabs_af),
  c("#7F312FFF", "#f7f7f7", "#005C55FF")
)

ha_af <- .spatial_make_ha(
  mat_nhood_z, af_meta, af_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = 9L
)
wh_af <- .spatial_hm_wh(mat_nhood_z)
png_af <- .spatial_png_inches(mat_nhood_z, pad_h = 4.4)

.cell_fun_nhood_stars <- function(mat_q) {
  function(j, i, x, y, w, h, fill) {
    qv <- mat_q[i, j]
    if (!is.finite(qv)) return(invisible(NULL))
    sym <- if (qv < 0.001) {
      "***"
    } else if (qv < 0.01) {
      "**"
    } else if (qv < 0.05) {
      "*"
    } else {
      return(invisible(NULL))
    }
    grid::grid.text(sym, x, y, gp = grid::gpar(fontsize = 9))
  }
}

ht_nhood_af <- ComplexHeatmap::Heatmap(
  mat_nhood_z,
  name = "Z",
  col = col_fun_nhood_af,
  width = wh_af$width,
  height = wh_af$height,
  na_col = "#eeeeee",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_af$row,
  top_annotation = ha_af$col,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Neighborhood\nenrichment z",
    9L
  ),
  cell_fun = .cell_fun_nhood_stars(mat_nhood_q)
)

.spatial_draw_ht_png(
  ht_nhood_af,
  file.path(out_dir, "spatial_colocalization_nhood_enrichment_adverse_favorable.png"),
  "Neighborhood enrichment of Adverse / Favorable cell types",
  "Excludes Other (non-tail) groups; FDR stars from full-tissue nhood test",
  width_in = png_af["width"],
  height_in = png_af["height"],
  legend_mode = "right_1col"
)

nhood_af_ord <- .cluster_heatmap_order(
  mat_nhood_z,
  col_fun_nhood_af,
  ha_af$row,
  ha_af$col,
  distance = "euclidean"
)
mat_nhood_z_cl <- mat_nhood_z[nhood_af_ord$row, nhood_af_ord$col, drop = FALSE]
mat_nhood_q_cl <- mat_nhood_q[nhood_af_ord$row, nhood_af_ord$col, drop = FALSE]
ha_af_cl <- .spatial_make_ha(
  mat_nhood_z_cl, af_meta, af_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = 9L
)
ht_nhood_af_cluster <- ComplexHeatmap::Heatmap(
  mat_nhood_z_cl,
  name = "Z",
  col = col_fun_nhood_af,
  width = wh_af$width,
  height = wh_af$height,
  na_col = "#eeeeee",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_af_cl$row,
  top_annotation = ha_af_cl$col,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Neighborhood\nenrichment z",
    9L
  ),
  cell_fun = .cell_fun_nhood_stars(mat_nhood_q_cl)
)

.spatial_draw_ht_png(
  ht_nhood_af_cluster,
  file.path(
    out_dir,
    "spatial_colocalization_nhood_enrichment_adverse_favorable_clustered.png"
  ),
  "Neighborhood enrichment of Adverse / Favorable cell types (clustered)",
  "Excludes Other (non-tail) groups; FDR stars from full-tissue nhood test",
  width_in = png_af["width"],
  height_in = png_af["height"],
  legend_mode = "right_1col"
)

# Local co-localization scores with high-score + dual CellChat outlines
local_outline_threshold <- 0.75
mabs_local_af <- suppressWarnings(max(as.numeric(mat_local_af), na.rm = TRUE))
if (!is.finite(mabs_local_af) || mabs_local_af <= 0) mabs_local_af <- 1
col_fun_local_af <- .spatial_coloc_score_col_fun(mabs_local_af)

ha_local_af <- .spatial_make_ha(
  mat_local_af, af_meta, af_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = 9L
)
local_af_ord <- .cluster_heatmap_order(
  mat_local_af,
  col_fun_local_af,
  ha_local_af$row,
  ha_local_af$col,
  distance = "euclidean"
)
mat_local_af_cl <- mat_local_af[local_af_ord$row, local_af_ord$col, drop = FALSE]
mat_dual_af_cl <- mat_dual_af[local_af_ord$row, local_af_ord$col, drop = FALSE]
ha_local_af_cl <- .spatial_make_ha(
  mat_local_af_cl, af_meta, af_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = 9L
)

.cell_fun_local_and_dual <- function(
    mat_local,
    mat_dual_flag,
    score_thr,
    dual_col
) {
  function(j, i, x, y, w, h, fill) {
    val <- mat_local[i, j]
    if (is.finite(val) && val >= score_thr) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
      )
    }
    if (isTRUE(mat_dual_flag[i, j])) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = dual_col, lwd = 2.2, fill = NA)
      )
    }
  }
}

ht_local_af_integrated <- ComplexHeatmap::Heatmap(
  mat_local_af_cl,
  name = "Mean score",
  col = col_fun_local_af,
  width = wh_af$width,
  height = wh_af$height,
  na_col = "#eeeeee",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_local_af_cl$row,
  top_annotation = ha_local_af_cl$col,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Co-localization\nscore",
    9L
  ),
  cell_fun = .cell_fun_local_and_dual(
    mat_local_af_cl,
    mat_dual_af_cl,
    local_outline_threshold,
    dual_outline_color
  )
)

.spatial_draw_ht_png(
  ht_local_af_integrated,
  file.path(
    out_dir,
    "spatial_colocalization_colocalization_scores_clustered_outlined_integrated.png"
  ),
  "Local co-localization of Adverse / Favorable cell types (clustered)",
  paste0(
    "score \u2265 ", local_outline_threshold,
    " black; spatial co-localization + CellChat purple"
  ),
  width_in = png_af["width"],
  height_in = png_af["height"],
  legend_mode = "right_1col",
  annotation_legend_list = list(
    .spatial_outline_legend(
      dual_col = dual_outline_color,
      score_label = paste0("score \u2265 ", local_outline_threshold),
      dual_label = "Spatial + CellChat"
    )
  )
)

message(
  "Adverse/Favorable nhood + coloc/CellChat overlays: ",
  sum(mat_dual_af, na.rm = TRUE), " dual-positive pairs; ",
  sum(mat_cc_sig_af, na.rm = TRUE), " CellChat-supported pairs."
)

message(
  "Done. Integrated heatmaps written to ", out_dir,
  " (", sum(pair_prognostic$dual_spatial_lr), " dual-positive / ",
  nrow(pair_prognostic), " prognostic pairs)."
)
