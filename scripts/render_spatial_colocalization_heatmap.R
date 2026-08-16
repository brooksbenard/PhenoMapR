#!/usr/bin/env Rscript
# Pre-render spatialCooccur figures for the spatial-transcriptomics vignette.
# Uses the full CytoSPACE cell set (no subsampling):
#   - spatialCooccur::nhood_enrichment() neighborhood heatmap
#   - spatialCooccur::cooccur_local() pairwise local co-localization heatmaps
#     (ordered and clustered)
#
# Run from package root:
#   Rscript scripts/render_spatial_colocalization_heatmap.R
#   Rscript scripts/render_spatial_colocalization_heatmap.R /path/to/PhenoMapR
#
# Output:
#   inst/figures/spatial_colocalization_nhood_enrichment.png
#   inst/figures/spatial_colocalization_nhood_enrichment_clustered.png
#   inst/figures/spatial_colocalization_colocalization_scores.png
#   inst/figures/spatial_colocalization_colocalization_scores_clustered.png
#   inst/figures/spatial_colocalization_colocalization_scores_clustered_outlined.png
#   inst/figures/spatial_colocalization_pheno_vs_cooccur_spearman.png
#   inst/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered.png
#   inst/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined.png
#   inst/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman.png
#   inst/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered.png
#   inst/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined.png

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[1L])
} else {
  normalizePath(".")
}

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
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
})

if (!requireNamespace("spatialCooccur", quietly = TRUE)) {
  stop(
    "Install spatialCooccur: remotes::install_github(\"juninamo/spatialCooccur\")",
    call. = FALSE
  )
}

source(file.path(root, "scripts", "spatial_colocal_heatmap_helpers.R"), local = TRUE)
source(file.path(root, "scripts", "spatial_colocal_palette.R"), local = TRUE)

.spatial_coloc_legend_title_size <- 9L

vdir <- file.path(root, "vignettes")
rds_cyto <- file.path(vdir, "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds")
gd_id_cyto <- "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll"

.download_spatial_rds <- function(dest, env_var = "PHENOMAPR_SPATIAL_RDS_URL") {
  if (file.exists(dest)) return(invisible(dest))
  url <- Sys.getenv(env_var, unset = "")
  if (nzchar(url)) {
    message("Downloading CytoSPACE RDS from ", env_var, " ...")
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    if (file.exists(dest)) return(invisible(dest))
  }
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    stop("Place ", basename(dest), " in vignettes/ or set ", env_var)
  }
  options(googledrive_quiet = TRUE)
  googledrive::drive_deauth()
  message("Downloading CytoSPACE RDS from Google Drive ...")
  googledrive::drive_download(googledrive::as_id(gd_id_cyto), dest, overwrite = TRUE)
  invisible(dest)
}

.download_spatial_rds(rds_cyto)
if (!file.exists(rds_cyto)) {
  stop("CytoSPACE RDS not found: ", rds_cyto)
}

message("Loading CytoSPACE object ...")
seurat <- readRDS(rds_cyto)
score_col <- grep("weighted_sum_score", names(seurat@meta.data), value = TRUE)[1]
if (is.na(score_col)) {
  stop("No weighted_sum_score column in CytoSPACE metadata")
}

cell_locations <- seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::select("Cell", "row", "col", "CellType", dplyr::all_of(score_col)) %>%
  dplyr::group_by(.data$row, .data$col) %>%
  dplyr::mutate(points_per_location = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    percentile = dplyr::percent_rank(.data[[score_col]]),
    prognostic_group = dplyr::case_when(
      .data$percentile < 0.05 ~ "Favorable",
      .data$percentile >= 0.95 ~ "Adverse",
      TRUE ~ "Other"
    )
  )

dfb <- cell_locations[
  stats::complete.cases(cell_locations[, c("Cell", "row", "col", "CellType", "prognostic_group")]),
  ,
  drop = FALSE
]
dfb$Cell <- as.character(dfb$Cell)
dfb$row <- as.numeric(dfb$row)
dfb$col <- as.numeric(dfb$col)
dfb$CellType <- as.character(dfb$CellType)
dfb$prognostic_group <- as.character(dfb$prognostic_group)
dfb <- dfb[dfb$prognostic_group %in% c("Adverse", "Favorable", "Other"), , drop = FALSE]
dfb <- dfb[!is.na(dfb$CellType) & nzchar(trimws(dfb$CellType)), , drop = FALSE]
dfb$ct_pg <- paste0(dfb$CellType, "_", dfb$prognostic_group)

scoc_df <- data.frame(
  x = dfb$row,
  y = dfb$col,
  ct_pg = dfb$ct_pg,
  stringsAsFactors = FALSE
)
rownames(scoc_df) <- dfb$Cell
scoc_df$cell_type <- scoc_df$ct_pg

message("Running nhood_enrichment on ", nrow(scoc_df), " cells (no subsampling) ...")

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
n_perm_sc <- 100L
n_grp <- length(unique(scoc_df$ct_pg))
if (nrow(scoc_df) < (k_sc + 3L) || n_grp < 2L) {
  stop("Not enough cells or label groups for nhood_enrichment")
}

scoc_nhood <- spatialCooccur::nhood_enrichment(
  scoc_df,
  cluster_key = "ct_pg",
  neighbors.k = k_sc,
  connectivity_key = "nn",
  transformation = TRUE,
  n_perms = n_perm_sc,
  seed = 42L,
  n_jobs = 1L
)

if (is.null(scoc_nhood$zscore)) {
  stop("nhood_enrichment returned no zscore matrix")
}

zm <- as.matrix(scoc_nhood$zscore)
lab_clean <- function(nm) gsub("^Cluster", "", nm)
rn <- lab_clean(rownames(zm))
cn <- lab_clean(colnames(zm))

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

all_labs <- sort(unique(c(rn, cn)))
meta_labs <- as.data.frame(
  do.call(rbind, lapply(all_labs, function(l) {
    p <- .parse_ct_pg(l)
    data.frame(label = l, CellType = p[["ct"]], Prognostic = p[["pg"]], stringsAsFactors = FALSE)
  })),
  stringsAsFactors = FALSE
)
meta_labs$Prognostic <- factor(meta_labs$Prognostic, levels = suffix_pg)
meta_labs <- meta_labs[!is.na(meta_labs$Prognostic), , drop = FALSE]
ord <- order(meta_labs$CellType, meta_labs$Prognostic)
labs_ord <- meta_labs$label[ord]
labs_ord <- labs_ord[labs_ord %in% rn & labs_ord %in% cn]
if (length(labs_ord) < 2L) {
  stop("Not enough overlapping labels in zscore matrix for heatmap")
}

mat <- zm[match(labs_ord, rn), match(labs_ord, cn), drop = FALSE]
dimnames(mat) <- list(labs_ord, labs_ord)

row_ct <- as.character(meta_labs$CellType[match(rownames(mat), meta_labs$label)])
row_pg <- as.character(meta_labs$Prognostic[match(rownames(mat), meta_labs$label)])
col_ct <- as.character(meta_labs$CellType[match(colnames(mat), meta_labs$label)])
col_pg <- as.character(meta_labs$Prognostic[match(colnames(mat), meta_labs$label)])

pal_pg <- c(Adverse = "#B2182B", Favorable = "#2166AC", Other = "#f7f7f7")
uct <- sort(unique(c(row_ct, col_ct)))
pal_ct <- PhenoMapR::get_celltype_palette(uct)

col_ncells <- rep(0L, ncol(mat))
names(col_ncells) <- colnames(mat)
tab_ct <- table(scoc_df$ct_pg)
hit <- names(col_ncells) %in% names(tab_ct)
col_ncells[hit] <- as.integer(tab_ct[names(col_ncells)[hit]])
.spatial_save_label_ncells(root, tab_ct)

mat_p <- 2 * stats::pnorm(-abs(mat))
mat_p[!is.finite(mat)] <- NA_real_
mat_p_adj <- mat_p
pv_flat <- as.vector(mat_p)
ok_flat <- is.finite(pv_flat)
mat_p_adj[] <- NA_real_
mat_p_adj[ok_flat] <- stats::p.adjust(pv_flat[ok_flat], method = "BH")

mabs <- suppressWarnings(max(abs(as.numeric(mat)), na.rm = TRUE))
if (!is.finite(mabs) || mabs <= 0) mabs <- 1
col_fun <- circlize::colorRamp2(c(-mabs, 0, mabs), c("#7F312FFF", "#f7f7f7", "#005C55FF"))

ha_nhood <- .spatial_make_ha(
  mat, meta_labs, col_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = .spatial_coloc_legend_title_size
)
row_ha <- ha_nhood$row
col_ha <- ha_nhood$col
hm_wh_nhood <- .spatial_hm_wh(mat)

ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Z",
  col = col_fun,
  width = hm_wh_nhood$width,
  height = hm_wh_nhood$height,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha,
  top_annotation = col_ha,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Neighborhood\nenrichment z",
    .spatial_coloc_legend_title_size
  ),
  border = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    qv <- mat_p_adj[i, j]
    if (!is.finite(qv)) return(invisible(NULL))
    sym <- if (qv < 0.001) "***" else if (qv < 0.01) "**" else if (qv < 0.05) "*" else return(invisible(NULL))
    grid::grid.text(sym, x, y, gp = grid::gpar(fontsize = 9))
  }
)

out_dir <- file.path(root, "inst", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(out_dir, "spatial_colocalization_nhood_enrichment.png")

message("Writing ", out_png, " ...")
png_sz <- .spatial_png_inches(mat)
grDevices::png(out_png, width = png_sz["width"], height = png_sz["height"], units = "in", res = 150)
.spatial_draw_heatmap(
  ht,
  column_title = "Neighborhood enrichment\nof prognostic cell types",
  column_title_gp = .spatial_title_gp(14),
  legend_mode = "right_1col"
)
grDevices::dev.off()

ht_nhood_cluster <- ComplexHeatmap::Heatmap(
  mat,
  name = "Z",
  col = col_fun,
  width = hm_wh_nhood$width,
  height = hm_wh_nhood$height,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha,
  top_annotation = col_ha,
  row_title = "Reference Cell Type",
  column_title = "Neighborhood Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Neighborhood\nenrichment z",
    .spatial_coloc_legend_title_size
  ),
  border = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    qv <- mat_p_adj[i, j]
    if (!is.finite(qv)) return(invisible(NULL))
    sym <- if (qv < 0.001) "***" else if (qv < 0.01) "**" else if (qv < 0.05) "*" else return(invisible(NULL))
    grid::grid.text(sym, x, y, gp = grid::gpar(fontsize = 9))
  }
)

out_png_nhood_cluster <- file.path(out_dir, "spatial_colocalization_nhood_enrichment_clustered.png")
message("Writing ", out_png_nhood_cluster, " ...")
grDevices::png(
  out_png_nhood_cluster,
  width = png_sz["width"],
  height = png_sz["height"],
  units = "in",
  res = 150
)
.spatial_draw_heatmap(
  ht_nhood_cluster,
  column_title = "Neighborhood enrichment\nof prognostic cell types (clustered)",
  column_title_gp = .spatial_title_gp(14),
  legend_mode = "right_1col"
)
grDevices::dev.off()

# Adverse / Favorable only (exclude *_Other) — same z-scores, filtered labels
af_keep <- rownames(mat)[
  as.character(meta_labs$Prognostic[match(rownames(mat), meta_labs$label)]) %in%
    c("Adverse", "Favorable")
]
if (length(af_keep) >= 2L) {
  mat_af <- mat[af_keep, af_keep, drop = FALSE]
  mat_p_adj_af <- mat_p_adj[af_keep, af_keep, drop = FALSE]
  col_ncells_af <- col_ncells[af_keep]
  meta_af <- meta_labs[meta_labs$label %in% af_keep, , drop = FALSE]
  ha_af <- .spatial_make_ha(
    mat_af, meta_af, col_ncells_af, pal_ct, pal_pg,
    anno_legend_ncol = 1L,
    legend_title_size = .spatial_coloc_legend_title_size
  )
  hm_wh_af <- .spatial_hm_wh(mat_af)
  png_sz_af <- .spatial_png_inches(mat_af)
  .cell_fun_af_stars <- function(j, i, x, y, w, h, fill) {
    qv <- mat_p_adj_af[i, j]
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
  ht_af <- ComplexHeatmap::Heatmap(
    mat_af,
    name = "Z",
    col = col_fun,
    width = hm_wh_af$width,
    height = hm_wh_af$height,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ha_af$row,
    top_annotation = ha_af$col,
    row_title = "Reference Cell Type",
    column_title = "Neighborhood Cell Type",
    row_title_gp = grid::gpar(fontsize = 12),
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    column_title_side = "bottom",
    heatmap_legend_param = .spatial_heatmap_legend_param(
      "Neighborhood\nenrichment z",
      .spatial_coloc_legend_title_size
    ),
    border = TRUE,
    cell_fun = .cell_fun_af_stars
  )
  out_png_af <- file.path(
    out_dir, "spatial_colocalization_nhood_enrichment_adverse_favorable.png"
  )
  message("Writing ", out_png_af, " ...")
  grDevices::png(
    out_png_af,
    width = png_sz_af["width"],
    height = png_sz_af["height"],
    units = "in",
    res = 150
  )
  .spatial_draw_heatmap(
    ht_af,
    column_title = "Neighborhood enrichment of Adverse / Favorable cell types",
    column_title_gp = .spatial_title_gp(14),
    legend_mode = "right_1col"
  )
  grDevices::dev.off()

  ht_af_cluster <- ComplexHeatmap::Heatmap(
    mat_af,
    name = "Z",
    col = col_fun,
    width = hm_wh_af$width,
    height = hm_wh_af$height,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ha_af$row,
    top_annotation = ha_af$col,
    row_title = "Reference Cell Type",
    column_title = "Neighborhood Cell Type",
    row_title_gp = grid::gpar(fontsize = 12),
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    column_title_side = "bottom",
    heatmap_legend_param = .spatial_heatmap_legend_param(
      "Neighborhood\nenrichment z",
      .spatial_coloc_legend_title_size
    ),
    border = TRUE,
    cell_fun = .cell_fun_af_stars
  )
  out_png_af_cluster <- file.path(
    out_dir,
    "spatial_colocalization_nhood_enrichment_adverse_favorable_clustered.png"
  )
  message("Writing ", out_png_af_cluster, " ...")
  grDevices::png(
    out_png_af_cluster,
    width = png_sz_af["width"],
    height = png_sz_af["height"],
    units = "in",
    res = 150
  )
  .spatial_draw_heatmap(
    ht_af_cluster,
    column_title = "Neighborhood enrichment of Adverse / Favorable cell types (clustered)",
    column_title_gp = .spatial_title_gp(14),
    legend_mode = "right_1col"
  )
  grDevices::dev.off()
}

# --- cooccur_local heatmap ---
message("Running cooccur_local on ", nrow(scoc_df), " cells ...")
rad_sc <- max(
  3,
  0.04 * max(
    diff(range(scoc_df$x, na.rm = TRUE)),
    diff(range(scoc_df$y, na.rm = TRUE)),
    na.rm = TRUE
  )
)
local_maxnsteps <- 3L
local_min_cells <- 3L

coords_local <- as.matrix(scoc_df[, c("x", "y"), drop = FALSE])
neighbors_local <- Seurat::FindNeighbors(coords_local, k.param = k_sc, verbose = FALSE)
adj_local <- neighbors_local$nn
rownames(adj_local) <- colnames(adj_local) <- rownames(scoc_df)
degrees_local <- Matrix::colSums(adj_local) + 1
res_local <- RANN::nn2(
  data = coords_local,
  query = coords_local,
  searchtype = "radius",
  radius = rad_sc,
  k = k_sc
)
label_vec <- as.character(scoc_df$cell_type)
neighbor_labels <- vector("list", nrow(scoc_df))
for (i in seq_len(nrow(scoc_df))) {
  ni <- res_local$nn.idx[i, ]
  ni <- ni[ni > 0 & ni != i]
  neighbor_labels[[i]] <- if (length(ni)) label_vec[ni] else character(0)
}

.spatial_cooccur_local_scores <- function(cluster_x, cluster_y) {
  binary <- vapply(
    neighbor_labels,
    function(nl) {
      if (!length(nl)) return(0)
      as.numeric(any(nl == cluster_x) & any(nl == cluster_y))
    },
    numeric(1)
  )
  if (local_maxnsteps <= 0L) {
    return(binary)
  }
  s_norm <- matrix(binary) / degrees_local
  as.numeric((adj_local %*% s_norm) + s_norm)
}

all_ref_labels <- labs_ord
mat_local <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)
ref_idx_by_label <- split(seq_len(nrow(scoc_df)), scoc_df$ct_pg)

for (ref_label in all_ref_labels) {
  ref_idx <- ref_idx_by_label[[ref_label]]
  if (is.null(ref_idx) || length(ref_idx) < local_min_cells) next

  for (target_label in all_ref_labels) {
    scores <- .spatial_cooccur_local_scores(ref_label, target_label)
    s_ref <- scores[ref_idx]
    s_nz <- s_ref[is.finite(s_ref) & s_ref > 0]
    mat_local[ref_label, target_label] <- if (length(s_nz)) mean(s_nz) else 0
  }
}

if (all(!is.finite(mat_local))) {
  stop("cooccur_local returned no scores for any reference group")
}

mabs_local <- suppressWarnings(max(as.numeric(mat_local), na.rm = TRUE))
if (!is.finite(mabs_local) || mabs_local <= 0) mabs_local <- 1
col_fun_local <- .spatial_coloc_score_col_fun(mabs_local)

ha_local <- .spatial_make_ha(
  mat_local, meta_labs, col_ncells, pal_ct, pal_pg,
  anno_legend_ncol = 1L,
  legend_title_size = .spatial_coloc_legend_title_size
)
row_ha_local <- ha_local$row
col_ha_local <- ha_local$col
hm_wh_local <- .spatial_hm_wh(mat_local)
png_sz_local <- .spatial_png_inches(mat_local)

ht_local <- ComplexHeatmap::Heatmap(
  mat_local,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  na_col = "#eeeeee",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Co-localization\nscore",
    .spatial_coloc_legend_title_size
  )
)

out_png_local <- file.path(out_dir, "spatial_colocalization_colocalization_scores.png")
message("Writing ", out_png_local, " ...")
grDevices::png(
  out_png_local,
  width = png_sz_local["width"],
  height = png_sz_local["height"],
  units = "in",
  res = 150
)
.spatial_draw_heatmap(
  ht_local,
  column_title = "Pairwise local co-localization scores\n(all prognostic cell-type groups)",
  column_title_gp = .spatial_title_gp(14),
  legend_mode = "right_1col"
)
grDevices::dev.off()

mat_local_cluster <- mat_local
mat_local_cluster[!is.finite(mat_local_cluster)] <- 0

ht_local_cluster <- ComplexHeatmap::Heatmap(
  mat_local_cluster,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Co-localization\nscore",
    .spatial_coloc_legend_title_size
  )
)
grDevices::pdf(NULL)
ht_local_cluster_drawn <- ComplexHeatmap::draw(ht_local_cluster)
local_row_ord <- ComplexHeatmap::row_order(ht_local_cluster_drawn)
local_col_ord <- ComplexHeatmap::column_order(ht_local_cluster_drawn)
grDevices::dev.off()

local_outline_threshold <- 0.75
ht_local_cluster_outline <- ComplexHeatmap::Heatmap(
  mat_local_cluster,
  name = "Mean score",
  col = col_fun_local,
  width = hm_wh_local$width,
  height = hm_wh_local$height,
  row_order = local_row_ord,
  column_order = local_col_ord,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 7),
  left_annotation = row_ha_local,
  top_annotation = col_ha_local,
  row_title = "Reference Cell Type",
  column_title = "Target Cell Type",
  row_title_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(fontsize = 12),
  row_title_side = "left",
  column_title_side = "bottom",
  border = TRUE,
  heatmap_legend_param = .spatial_heatmap_legend_param(
    "Co-localization\nscore",
    .spatial_coloc_legend_title_size
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat_local_cluster[i, j]
    if (is.finite(val) && val >= local_outline_threshold) {
      grid::grid.rect(
        x, y, w, h,
        gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
      )
    }
  }
)

out_png_local_cluster <- file.path(
  out_dir, "spatial_colocalization_colocalization_scores_clustered.png"
)
message("Writing ", out_png_local_cluster, " ...")
grDevices::png(
  out_png_local_cluster,
  width = png_sz_local["width"],
  height = png_sz_local["height"],
  units = "in",
  res = 150
)
.spatial_draw_heatmap(
  ht_local_cluster,
  column_title = "Pairwise local co-localization scores\n(all prognostic cell-type groups; clustered)",
  column_title_gp = .spatial_title_gp(14),
  legend_mode = "right_1col"
)
grDevices::dev.off()

out_png_local_cluster_outline <- file.path(
  out_dir,
  "spatial_colocalization_colocalization_scores_clustered_outlined.png"
)
message("Writing ", out_png_local_cluster_outline, " ...")
grDevices::png(
  out_png_local_cluster_outline,
  width = png_sz_local["width"],
  height = png_sz_local["height"],
  units = "in",
  res = 150
)
.spatial_draw_heatmap(
  ht_local_cluster_outline,
  column_title = paste0(
    "Pairwise local co-localization scores\n(clustered; score \u2265 ",
    local_outline_threshold, " outlined)"
  ),
  column_title_gp = .spatial_title_gp(14),
  legend_mode = "right_1col"
)
grDevices::dev.off()

# --- PhenoMapR score vs local co-localization ---
message("Computing PhenoMapR vs co-localization correlations ...")
cor_min_cells <- 10L
pheno_vec <- stats::setNames(as.numeric(dfb[[score_col]]), dfb$Cell)
pheno_vec <- pheno_vec[rownames(scoc_df)]

neighbor_idx <- vector("list", nrow(scoc_df))
for (i in seq_len(nrow(scoc_df))) {
  ni <- res_local$nn.idx[i, ]
  neighbor_idx[[i]] <- ni[ni > 0 & ni != i]
}

mat_nbr_pheno <- matrix(
  NA_real_,
  nrow = nrow(scoc_df),
  ncol = length(all_ref_labels),
  dimnames = list(rownames(scoc_df), all_ref_labels)
)
for (i in seq_len(nrow(scoc_df))) {
  ni <- neighbor_idx[[i]]
  if (!length(ni)) next
  nl <- label_vec[ni]
  pv <- pheno_vec[ni]
  for (target_label in unique(nl)) {
    hit <- nl == target_label
    if (any(hit)) {
      mat_nbr_pheno[i, target_label] <- mean(pv[hit], na.rm = TRUE)
    }
  }
}

mat_cor_pheno_cooccur <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)
mat_cor_pheno_neighbor <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)

.spatial_spearman_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < cor_min_cells) return(NA_real_)
  x_ok <- x[ok]
  y_ok <- y[ok]
  if (stats::sd(x_ok) == 0 || stats::sd(y_ok) == 0) return(NA_real_)
  stats::cor(x_ok, y_ok, method = "spearman")
}

for (ref_label in all_ref_labels) {
  ref_idx <- ref_idx_by_label[[ref_label]]
  if (is.null(ref_idx) || length(ref_idx) < cor_min_cells) next
  pheno_ref <- pheno_vec[ref_idx]
  for (target_label in all_ref_labels) {
    cooc_ref <- .spatial_cooccur_local_scores(ref_label, target_label)[ref_idx]
    mat_cor_pheno_cooccur[ref_label, target_label] <- .spatial_spearman_cor(pheno_ref, cooc_ref)
    nbr_ref <- mat_nbr_pheno[ref_idx, target_label]
    mat_cor_pheno_neighbor[ref_label, target_label] <- .spatial_spearman_cor(pheno_ref, nbr_ref)
  }
}

col_fun_cor <- .spatial_spearman_col_fun()

.write_cor_heatmap <- function(
    mat_cor,
    out_path,
    title,
    clustered = FALSE,
    outline_abs_threshold = NULL
) {
  ha_cor <- .spatial_make_ha(
    mat_cor, meta_labs, col_ncells, pal_ct, pal_pg,
    anno_legend_ncol = 1L,
    legend_title_size = .spatial_coloc_legend_title_size
  )
  hm_wh <- .spatial_hm_wh(mat_cor)
  png_sz <- .spatial_png_inches(mat_cor)
  row_order <- NULL
  col_order <- NULL
  if (isTRUE(clustered) && !is.null(outline_abs_threshold)) {
    ht_cluster <- ComplexHeatmap::Heatmap(
      mat_cor,
      name = "Spearman",
      col = col_fun_cor,
      width = hm_wh$width,
      height = hm_wh$height,
      na_col = "#eeeeee",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      clustering_distance_rows = .spatial_cor_dist,
      clustering_distance_columns = .spatial_cor_dist,
      show_row_names = FALSE,
      show_column_names = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      left_annotation = ha_cor$row,
      top_annotation = ha_cor$col,
      border = TRUE
    )
    grDevices::pdf(NULL)
    ht_cluster_drawn <- ComplexHeatmap::draw(ht_cluster)
    row_order <- ComplexHeatmap::row_order(ht_cluster_drawn)
    col_order <- ComplexHeatmap::column_order(ht_cluster_drawn)
    grDevices::dev.off()
    clustered <- FALSE
  }
  cell_fun_outline <- NULL
  if (!is.null(outline_abs_threshold)) {
    thr <- outline_abs_threshold
    cell_fun_outline <- function(j, i, x, y, w, h, fill) {
      val <- mat_cor[i, j]
      if (is.finite(val) && abs(val) > thr) {
        grid::grid.rect(
          x, y, w, h,
          gp = grid::gpar(col = "black", lwd = 1.5, fill = NA)
        )
      }
    }
  }
  ht_cor <- ComplexHeatmap::Heatmap(
    mat_cor,
    name = "Spearman",
    col = col_fun_cor,
    width = hm_wh$width,
    height = hm_wh$height,
    na_col = "#eeeeee",
    cluster_rows = clustered,
    cluster_columns = clustered,
    row_order = row_order,
    column_order = col_order,
    clustering_distance_rows = if (clustered) .spatial_cor_dist else "euclidean",
    clustering_distance_columns = if (clustered) .spatial_cor_dist else "euclidean",
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
    left_annotation = ha_cor$row,
    top_annotation = ha_cor$col,
    row_title = "Reference Cell Type",
    column_title = "Target Cell Type",
    row_title_gp = grid::gpar(fontsize = 12),
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    column_title_side = "bottom",
    border = TRUE,
    heatmap_legend_param = .spatial_heatmap_legend_param(
      "Spearman rho",
      .spatial_coloc_legend_title_size
    ),
    cell_fun = cell_fun_outline
  )
  message("Writing ", out_path, " ...")
  grDevices::png(out_path, width = png_sz["width"], height = png_sz["height"], units = "in", res = 150)
  .spatial_draw_heatmap(
    ht_cor,
    column_title = title,
    column_title_gp = .spatial_title_gp(14),
    legend_mode = "right_1col"
  )
  grDevices::dev.off()
}

cor_outline_threshold <- 0.5

.write_cor_heatmap(
  mat_cor_pheno_cooccur,
  file.path(out_dir, "spatial_colocalization_pheno_vs_cooccur_spearman.png"),
  "PhenoMapR vs local co-localization\n(Spearman rho)",
  clustered = FALSE
)
.write_cor_heatmap(
  mat_cor_pheno_cooccur,
  file.path(out_dir, "spatial_colocalization_pheno_vs_cooccur_spearman_clustered.png"),
  "PhenoMapR vs local co-localization\n(Spearman rho, clustered)",
  clustered = TRUE
)
.write_cor_heatmap(
  mat_cor_pheno_cooccur,
  file.path(out_dir, "spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined.png"),
  paste0(
    "PhenoMapR vs local co-localization\n(clustered; |rho| > ",
    cor_outline_threshold, " outlined)"
  ),
  clustered = TRUE,
  outline_abs_threshold = cor_outline_threshold
)
.write_cor_heatmap(
  mat_cor_pheno_neighbor,
  file.path(out_dir, "spatial_colocalization_pheno_vs_neighbor_pheno_spearman.png"),
  "PhenoMapR vs neighbor PhenoMapR\n(Spearman rho)",
  clustered = FALSE
)
.write_cor_heatmap(
  mat_cor_pheno_neighbor,
  file.path(out_dir, "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered.png"),
  "PhenoMapR vs neighbor PhenoMapR\n(Spearman rho, clustered)",
  clustered = TRUE
)
.write_cor_heatmap(
  mat_cor_pheno_neighbor,
  file.path(out_dir, "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined.png"),
  paste0(
    "PhenoMapR vs neighbor PhenoMapR\n(clustered; |rho| > ",
    cor_outline_threshold, " outlined)"
  ),
  clustered = TRUE,
  outline_abs_threshold = cor_outline_threshold
)

message("Done. Re-knit vignettes/spatial-transcriptomics.Rmd to pick up the figures.")
