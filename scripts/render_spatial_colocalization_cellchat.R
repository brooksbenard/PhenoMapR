#!/usr/bin/env Rscript
# Integrate spatial co-localization / PhenoMapR correlations with spatial
# CellChat ligand-receptor inference for CytoSPACE-mapped prognostic cells.
#
# Run from package root:
#   Rscript scripts/render_spatial_colocalization_cellchat.R
#
# Requires: spatialCooccur, CellChat, Seurat, RANN, ComplexHeatmap (optional for heatmap)
# Output:
#   inst/figures/spatial_coloc_cellchat_bubble.png
#   inst/figures/spatial_coloc_cellchat_dual_heatmap.png
#   inst/figures/spatial_coloc_cellchat_chord.png
#   inst/figures/spatial_coloc_cellchat_lr_bubble.png
#   inst/figures/spatial_coloc_cellchat_type_spatial.png
#   inst/figures/spatial_coloc_cellchat_type_pheno.png
#   inst/figures/spatial_coloc_cellchat_type_assoc.png
#   inst/figures/spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined_integrated.png
#   inst/figures/spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined_integrated.png
#   inst/figures/spatial_coloc_integrated_four_panel.png
#   inst/figures/spatial_coloc_cellchat_prob_heatmap_clustered.png
#   inst/figures/spatial_coloc_integration_evidence.png
#   inst/data/spatial_coloc_cellchat_pairs.rds
#   inst/data/spatial_coloc_cellchat_lr_pairs.rds

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
  library(tidyr)
  library(tibble)
})

for (pkg in c("spatialCooccur", "CellChat", "RANN", "circlize")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Install ", pkg, " before running this script.",
      if (identical(pkg, "spatialCooccur")) {
        " remotes::install_github(\"juninamo/spatialCooccur\")"
      } else if (identical(pkg, "CellChat")) {
        " remotes::install_github(\"jinworks/CellChat\")"
      } else {
        ""
      },
      call. = FALSE
    )
  }
}

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

# Build Visium tissue coordinates for CellChat spatial mode from the CytoSPACE object.
# Each mapped cell inherits the Visium spot tissue position (imagerow/imagecol) of its
# assigned CytoSPACE location; distance constraints follow CellChat Visium defaults.
.build_cytospace_cellchat_spatial <- function(seurat, cell_ids) {
  if (!all(c("row", "col") %in% names(seurat@meta.data))) {
    stop("CytoSPACE metadata must include row/col Visium grid coordinates")
  }
  img_name <- Seurat::Images(seurat)[1]
  if (is.na(img_name) || !nzchar(img_name)) {
    stop("Seurat object has no Visium image for spatial CellChat coordinates")
  }
  spatial.locs <- Seurat::GetTissueCoordinates(
    seurat,
    scale = NULL,
    cols = c("imagerow", "imagecol")
  )[cell_ids, , drop = FALSE]
  colnames(spatial.locs) <- c("imagerow", "imagecol")
  spot.size <- 65
  spot_diameter <- seurat@images[[img_name]]@scale.factors$spot
  if (!is.finite(spot_diameter) || spot_diameter <= 0) {
    stop("Visium spot_diameter missing from image scale.factors")
  }
  spatial.factors <- data.frame(
    ratio = spot.size / spot_diameter,
    tol = spot.size / 2
  )
  d_um <- CellChat::computeCellDistance(
    coordinates = spatial.locs,
    ratio = spatial.factors$ratio,
    tol = spatial.factors$tol
  )
  min_d <- min(d_um[d_um > 0], na.rm = TRUE)
  if (!is.finite(min_d) || min_d <= 0) min_d <- spot.size * 1.5
  scale.distance <- min(1, 1.1 / min_d)
  list(
    coordinates = spatial.locs,
    spatial.factors = spatial.factors,
    scale.distance = scale.distance,
    contact.range = 100
  )
}

.download_spatial_rds(rds_cyto)
if (!file.exists(rds_cyto)) stop("CytoSPACE RDS not found: ", rds_cyto)

message("Loading CytoSPACE object ...")
seurat <- readRDS(rds_cyto)
score_col <- grep("weighted_sum_score", names(seurat@meta.data), value = TRUE)[1]
if (is.na(score_col)) stop("No weighted_sum_score column in CytoSPACE metadata")

cell_locations <- seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::select("Cell", "row", "col", "CellType", dplyr::all_of(score_col)) %>%
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

message("Running nhood_enrichment on ", nrow(scoc_df), " cells ...")
scoc_nhood <- spatialCooccur::nhood_enrichment(
  scoc_df,
  cluster_key = "ct_pg",
  neighbors.k = k_sc,
  connectivity_key = "nn",
  transformation = TRUE,
  n_perms = 100L,
  seed = 42L,
  n_jobs = 1L
)
zm <- as.matrix(scoc_nhood$zscore)
lab_clean <- function(nm) gsub("^Cluster", "", nm)
rn <- lab_clean(rownames(zm))
cn <- lab_clean(colnames(zm))
all_labs <- sort(unique(c(rn, cn)))
meta_labs <- as.data.frame(
  do.call(rbind, lapply(all_labs, function(l) {
    p <- .parse_ct_pg(l)
    data.frame(label = l, CellType = p[["ct"]], Prognostic = p[["pg"]], stringsAsFactors = FALSE)
  })),
  stringsAsFactors = FALSE
)
meta_labs <- meta_labs[!is.na(meta_labs$Prognostic), , drop = FALSE]
labs_ord <- meta_labs$label[order(meta_labs$CellType, meta_labs$Prognostic)]
labs_ord <- labs_ord[labs_ord %in% rn & labs_ord %in% cn]

source(file.path(root, "scripts", "spatial_colocal_heatmap_helpers.R"), local = TRUE)
.spatial_save_label_ncells(root, table(scoc_df$ct_pg))

mat_z <- zm[match(labs_ord, rn), match(labs_ord, cn), drop = FALSE]
dimnames(mat_z) <- list(labs_ord, labs_ord)
mat_p <- 2 * stats::pnorm(-abs(mat_z))
mat_p[!is.finite(mat_p)] <- NA_real_
mat_p_adj <- mat_p
ok_flat <- is.finite(as.vector(mat_p))
mat_p_adj[] <- NA_real_
mat_p_adj[ok_flat] <- stats::p.adjust(as.vector(mat_p)[ok_flat], method = "BH")

message("Computing local co-localization scores ...")
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
cor_min_cells <- 10L

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
  if (local_maxnsteps <= 0L) return(binary)
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

pheno_vec <- stats::setNames(as.numeric(dfb[[score_col]]), dfb$Cell)
pheno_vec <- pheno_vec[rownames(scoc_df)]
mat_nbr_pheno <- matrix(
  NA_real_,
  nrow = nrow(scoc_df),
  ncol = length(all_ref_labels),
  dimnames = list(rownames(scoc_df), all_ref_labels)
)
for (i in seq_len(nrow(scoc_df))) {
  ni <- res_local$nn.idx[i, ]
  ni <- ni[ni > 0 & ni != i]
  if (!length(ni)) next
  nl <- label_vec[ni]
  pv <- pheno_vec[ni]
  for (target_label in unique(nl)) {
    hit <- nl == target_label
    if (any(hit)) mat_nbr_pheno[i, target_label] <- mean(pv[hit], na.rm = TRUE)
  }
}

mat_cor_pheno_cooccur <- matrix(
  NA_real_,
  nrow = length(all_ref_labels),
  ncol = length(all_ref_labels),
  dimnames = list(all_ref_labels, all_ref_labels)
)
mat_cor_pheno_neighbor <- mat_cor_pheno_cooccur
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
    mat_cor_pheno_neighbor[ref_label, target_label] <- .spatial_spearman_cor(
      pheno_ref, mat_nbr_pheno[ref_idx, target_label]
    )
  }
}

.build_pair_summary <- function() {
  labs <- labs_ord
  out <- vector("list", length(labs) * length(labs))
  k <- 0L
  for (i in seq_along(labs)) {
    for (j in seq_along(labs)) {
      src <- labs[i]
      tgt <- labs[j]
      ps <- .parse_ct_pg(src)
      pt <- .parse_ct_pg(tgt)
      k <- k + 1L
      out[[k]] <- data.frame(
        source = src,
        target = tgt,
        source_ct = ps[["ct"]],
        target_ct = pt[["ct"]],
        source_pg = ps[["pg"]],
        target_pg = pt[["pg"]],
        nhood_z = mat_z[src, tgt],
        nhood_q = mat_p_adj[src, tgt],
        local_cooccur = mat_local[src, tgt],
        pheno_cooccur_rho = mat_cor_pheno_cooccur[src, tgt],
        pheno_neighbor_rho = mat_cor_pheno_neighbor[src, tgt],
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(out)
}

pair_df <- .build_pair_summary()
pair_df$prognostic_pair <- pair_df$source_pg %in% c("Adverse", "Favorable") &
  pair_df$target_pg %in% c("Adverse", "Favorable")
pair_df$cross_prognostic <- pair_df$source_pg != pair_df$target_pg &
  pair_df$source_pg %in% c("Adverse", "Favorable") &
  pair_df$target_pg %in% c("Adverse", "Favorable")
pair_df$spatial_sig <- is.finite(pair_df$nhood_z) & pair_df$nhood_z > 0 &
  is.finite(pair_df$nhood_q) & pair_df$nhood_q < 0.05
pair_df$local_high <- is.finite(pair_df$local_cooccur) & pair_df$local_cooccur >= 0.5
pair_df$pheno_sig <- (is.finite(pair_df$pheno_cooccur_rho) & abs(pair_df$pheno_cooccur_rho) >= 0.15) |
  (is.finite(pair_df$pheno_neighbor_rho) & abs(pair_df$pheno_neighbor_rho) >= 0.15)

message("Running spatial CellChat on CytoSPACE prognostic tail cells ...")
meta_all <- seurat@meta.data
meta_all$percentile <- dplyr::percent_rank(meta_all[[score_col]])
meta_all$prognostic_group <- dplyr::case_when(
  meta_all$percentile < 0.05 ~ "Favorable",
  meta_all$percentile >= 0.95 ~ "Adverse",
  TRUE ~ "Other"
)
meta_all$ct_pg <- paste0(meta_all$CellType, "_", meta_all$prognostic_group)
cells_pg <- rownames(meta_all)[meta_all$prognostic_group %in% c("Adverse", "Favorable")]
expr_counts <- Seurat::GetAssayData(seurat, assay = "Spatial", layer = "counts")[, cells_pg, drop = FALSE]
cs <- Matrix::colSums(expr_counts)
expr_norm <- Matrix::t(Matrix::t(expr_counts) / pmax(cs, 1) * 10000)
expr_norm <- log1p(expr_norm)
meta_pg <- meta_all[cells_pg, , drop = FALSE]
meta_pg$labels <- meta_pg$ct_pg
tab_pg <- table(meta_pg$labels)
keep_labs <- names(tab_pg)[tab_pg >= 10L]
cells_cc <- rownames(meta_pg)[meta_pg$labels %in% keep_labs]
expr_cc <- expr_norm[, cells_cc, drop = FALSE]
meta_cc <- meta_pg[cells_cc, , drop = FALSE]
meta_cc$samples <- factor(meta_cc$orig.ident)

spatial_inputs <- .build_cytospace_cellchat_spatial(seurat, cells_cc)
message(
  "Spatial CellChat: ", nrow(spatial_inputs$coordinates), " cells; ",
  "scale.distance = ", signif(spatial_inputs$scale.distance, 3),
  "; contact.range = ", spatial_inputs$contact.range, " um"
)

cellchat <- CellChat::createCellChat(
  object = expr_cc,
  meta = meta_cc,
  group.by = "labels",
  datatype = "spatial",
  coordinates = spatial_inputs$coordinates,
  spatial.factors = spatial_inputs$spatial.factors
)
cellchat@DB <- CellChat::CellChatDB.human
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = spatial_inputs$scale.distance,
  contact.dependent = TRUE,
  contact.range = spatial_inputs$contact.range
)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)
dir.create(file.path(root, "inst", "data"), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  cellchat,
  file.path(root, "inst", "data", "spatial_cellchat_object.rds")
)

cc_df <- CellChat::subsetCommunication(cellchat)

lr_anno <- CellChat::CellChatDB.human$interaction %>%
  dplyr::select(
    interaction_name,
    annotation
  ) %>%
  dplyr::mutate(
    interaction_mode = dplyr::case_when(
      .data$annotation %in% c("Cell-Cell Contact", "ECM-Receptor") ~ "physical",
      .data$annotation == "Secreted Signaling" ~ "secreted",
      TRUE ~ "other"
    )
  )

cc_df <- cc_df %>%
  dplyr::left_join(lr_anno, by = "interaction_name")

cc_agg <- cc_df %>%
  dplyr::group_by(source, target) %>%
  dplyr::summarise(
    cellchat_prob = sum(.data$prob, na.rm = TRUE),
    cellchat_n_lr = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    source = as.character(.data$source),
    target = as.character(.data$target)
  )

cc_agg_mode <- cc_df %>%
  dplyr::filter(.data$interaction_mode %in% c("physical", "secreted")) %>%
  dplyr::group_by(.data$source, .data$target, .data$interaction_mode) %>%
  dplyr::summarise(
    prob = sum(.data$prob, na.rm = TRUE),
    n_lr = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = .data$interaction_mode,
    values_from = c(.data$prob, .data$n_lr),
    values_fill = 0
  ) %>%
  dplyr::mutate(
    source = as.character(.data$source),
    target = as.character(.data$target),
    cellchat_prob_physical = dplyr::coalesce(.data$prob_physical, 0),
    cellchat_prob_secreted = dplyr::coalesce(.data$prob_secreted, 0),
    cellchat_n_physical = dplyr::coalesce(as.integer(.data$n_lr_physical), 0L),
    cellchat_n_secreted = dplyr::coalesce(as.integer(.data$n_lr_secreted), 0L)
  ) %>%
  dplyr::select(
    source, target,
    cellchat_prob_physical, cellchat_prob_secreted,
    cellchat_n_physical, cellchat_n_secreted
  )

pair_df <- pair_df %>%
  dplyr::left_join(cc_agg, by = c("source", "target")) %>%
  dplyr::left_join(cc_agg_mode, by = c("source", "target")) %>%
  dplyr::mutate(
    cellchat_prob = dplyr::coalesce(.data$cellchat_prob, 0),
    cellchat_n_lr = dplyr::coalesce(.data$cellchat_n_lr, 0L),
    cellchat_prob_physical = dplyr::coalesce(.data$cellchat_prob_physical, 0),
    cellchat_prob_secreted = dplyr::coalesce(.data$cellchat_prob_secreted, 0),
    cellchat_n_physical = dplyr::coalesce(.data$cellchat_n_physical, 0L),
    cellchat_n_secreted = dplyr::coalesce(.data$cellchat_n_secreted, 0L),
    cellchat_sig = .data$cellchat_prob > 0 & .data$cellchat_n_lr >= 1L,
    cellchat_sig_physical = .data$cellchat_prob_physical > 0,
    cellchat_sig_secreted = .data$cellchat_prob_secreted > 0,
    cellchat_frac_physical = dplyr::if_else(
      .data$cellchat_prob > 0,
      .data$cellchat_prob_physical / .data$cellchat_prob,
      NA_real_
    ),
    spatial_coloc_nhood = .data$spatial_sig,
    spatial_coloc_local = .data$local_high,
    spatial_coloc = .data$spatial_sig | .data$local_high,
    interaction_dominance = dplyr::case_when(
      .data$cellchat_prob <= 0 ~ "none",
      .data$cellchat_prob_physical >= 1.5 * pmax(.data$cellchat_prob_secreted, 1e-9) ~ "physical-dominant",
      .data$cellchat_prob_secreted >= 1.5 * pmax(.data$cellchat_prob_physical, 1e-9) ~ "secreted-dominant",
      TRUE ~ "mixed"
    ),
    dual_spatial_lr = .data$spatial_coloc & .data$cellchat_sig,
    dual_nhood_lr = .data$spatial_coloc_nhood & .data$cellchat_sig,
    dual_local_lr = .data$spatial_coloc_local & .data$cellchat_sig,
    integrated_score = {
      z1 <- as.numeric(scale(.data$nhood_z))
      z2 <- as.numeric(scale(log1p(pmax(.data$local_cooccur, 0))))
      z3 <- as.numeric(scale(log1p(.data$cellchat_prob)))
      z1[is.na(z1)] <- 0
      z2[is.na(z2)] <- 0
      z3[is.na(z3)] <- 0
      z1 + z2 + z3
    }
  )

pair_prognostic <- pair_df %>%
  dplyr::filter(.data$prognostic_pair) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score))

out_dir <- file.path(root, "inst", "figures")
data_dir <- file.path(root, "inst", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(pair_prognostic, file.path(data_dir, "spatial_coloc_cellchat_pairs.rds"))

lr_pair_df <- cc_df %>%
  dplyr::inner_join(
    pair_prognostic %>%
      dplyr::select(
        source, target, nhood_z, nhood_q, local_cooccur,
        pheno_cooccur_rho, pheno_neighbor_rho,
        spatial_sig, local_high, spatial_coloc_nhood, spatial_coloc_local
      ),
    by = c("source", "target")
  )
saveRDS(lr_pair_df, file.path(data_dir, "spatial_coloc_cellchat_lr_pairs.rds"))

.spatial_spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L) return(NA_real_)
  stats::cor(x[ok], y[ok], method = "spearman")
}

pair_cc <- pair_prognostic %>% dplyr::filter(.data$cellchat_prob > 0)
assoc_tbl <- tibble::tibble(
  comparison = c(
    "Physical prob ~ nhood z",
    "Physical prob ~ local co-localization",
    "Secreted prob ~ nhood z",
    "Secreted prob ~ local co-localization",
    "Physical prob ~ pheno vs co-localization rho",
    "Physical prob ~ pheno vs neighbor rho",
    "Secreted prob ~ pheno vs co-localization rho",
    "Secreted prob ~ pheno vs neighbor rho",
    "Physical fraction ~ nhood z",
    "Physical fraction ~ local co-localization",
    "Physical fraction ~ pheno vs neighbor rho"
  ),
  spearman_rho = c(
    .spatial_spearman_safe(pair_cc$nhood_z, pair_cc$cellchat_prob_physical),
    .spatial_spearman_safe(pair_cc$local_cooccur, pair_cc$cellchat_prob_physical),
    .spatial_spearman_safe(pair_cc$nhood_z, pair_cc$cellchat_prob_secreted),
    .spatial_spearman_safe(pair_cc$local_cooccur, pair_cc$cellchat_prob_secreted),
    .spatial_spearman_safe(pair_cc$pheno_cooccur_rho, pair_cc$cellchat_prob_physical),
    .spatial_spearman_safe(pair_cc$pheno_neighbor_rho, pair_cc$cellchat_prob_physical),
    .spatial_spearman_safe(pair_cc$pheno_cooccur_rho, pair_cc$cellchat_prob_secreted),
    .spatial_spearman_safe(pair_cc$pheno_neighbor_rho, pair_cc$cellchat_prob_secreted),
    .spatial_spearman_safe(pair_cc$nhood_z, pair_cc$cellchat_frac_physical),
    .spatial_spearman_safe(pair_cc$local_cooccur, pair_cc$cellchat_frac_physical),
    .spatial_spearman_safe(pair_cc$pheno_neighbor_rho, pair_cc$cellchat_frac_physical)
  )
)
saveRDS(assoc_tbl, file.path(data_dir, "spatial_coloc_cellchat_assoc.rds"))

top_dual <- pair_prognostic %>%
  dplyr::filter(.data$dual_spatial_lr) %>%
  dplyr::arrange(dplyr::desc(.data$integrated_score)) %>%
  dplyr::slice_head(n = 25)

pair_plot <- pair_prognostic %>%
  dplyr::filter(
    is.finite(.data$nhood_z),
    is.finite(.data$cellchat_prob)
  )

message("Writing spatial_coloc_cellchat_bubble.png ...")
p_bubble <- ggplot(
  pair_plot,
  aes(
    x = .data$nhood_z,
    y = log1p(.data$cellchat_prob),
    size = pmax(.data$local_cooccur, 0),
    color = .data$dual_spatial_lr
  )
) +
  geom_point(alpha = 0.65) +
  scale_color_manual(
    values = c(`TRUE` = "#B2182B", `FALSE` = "#666666"),
    name = "Spatial co-localization\n+ CellChat L-R"
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Local co-localization\n(mean non-zero score)") +
  labs(
    title = "Spatial co-localization vs distance-constrained CellChat signaling",
    subtitle = "CytoSPACE cell locations; prognostic pairs; x = neighborhood enrichment z-score",
    x = "Neighborhood enrichment z-score",
    y = "log1p(CellChat aggregated probability)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
if (nrow(top_dual) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  p_bubble <- p_bubble +
    ggrepel::geom_text_repel(
      data = top_dual,
      aes(label = paste0(.data$source, " \u2192 ", .data$target)),
      size = 2.6,
      max.overlaps = 12,
      segment.size = 0.2,
      show.legend = FALSE
    )
}
ggsave(
  file.path(out_dir, "spatial_coloc_cellchat_bubble.png"),
  p_bubble,
  width = 10,
  height = 7,
  dpi = 150
)

message("Writing spatial_coloc_cellchat_dual_heatmap.png ...")
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)) {
  labs_pg <- sort(unique(c(
    pair_prognostic$source[pair_prognostic$prognostic_pair],
    pair_prognostic$target[pair_prognostic$prognostic_pair]
  )))
  labs_pg <- labs_pg[labs_pg %in% keep_labs]
  if (length(labs_pg) >= 2L) {
    dual_mat <- matrix(0, length(labs_pg), length(labs_pg), dimnames = list(labs_pg, labs_pg))
    for (i in seq_len(nrow(pair_prognostic))) {
      r <- pair_prognostic$source[i]
      c <- pair_prognostic$target[i]
      if (r %in% labs_pg && c %in% labs_pg) {
        dual_mat[r, c] <- as.numeric(pair_prognostic$integrated_score[i])
      }
    }
    dual_mat[!is.finite(dual_mat)] <- 0
    col_fun <- circlize::colorRamp2(
      c(min(dual_mat), 0, max(dual_mat)),
      c("#2166AC", "#F7F7F7", "#B2182B")
    )
    ht_dual <- ComplexHeatmap::Heatmap(
      dual_mat,
      name = "Integrated",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_names_gp = grid::gpar(fontsize = 7),
      column_title = "Integrated spatial + spatial CellChat score (prognostic pairs)",
      heatmap_legend_param = list(title = "Score")
    )
    grDevices::png(
      file.path(out_dir, "spatial_coloc_cellchat_dual_heatmap.png"),
      width = 9,
      height = 8,
      units = "in",
      res = 150
    )
    ComplexHeatmap::draw(ht_dual)
    grDevices::dev.off()
  }
}

message("Writing spatial_coloc_cellchat_chord.png ...")
chord_groups <- unique(c(top_dual$source, top_dual$target))
if (length(chord_groups) >= 2L) {
  grDevices::png(
    file.path(out_dir, "spatial_coloc_cellchat_chord.png"),
    width = 9,
    height = 9,
    units = "in",
    res = 150
  )
  chord_ok <- FALSE
  tryCatch({
    pair_lr <- cc_df %>%
      dplyr::filter(.data$source %in% chord_groups, .data$target %in% chord_groups) %>%
      dplyr::group_by(.data$source, .data$target) %>%
      dplyr::summarise(prob = sum(.data$prob, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(.data$prob > 0)
    if (nrow(pair_lr) >= 2L) {
      mat_chord <- matrix(0, length(chord_groups), length(chord_groups), dimnames = list(chord_groups, chord_groups))
      for (i in seq_len(nrow(pair_lr))) {
        mat_chord[pair_lr$source[i], pair_lr$target[i]] <- pair_lr$prob[i]
      }
      pal_chord <- stats::setNames(rep("#888888", length(chord_groups)), chord_groups)
      pal_chord[endsWith(chord_groups, "_Adverse")] <- "#B2182B"
      pal_chord[endsWith(chord_groups, "_Favorable")] <- "#2166AC"
      circlize::chordDiagram(
        mat_chord,
        grid.col = pal_chord,
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = 0.1)
      )
      circlize::circos.clear()
      title("Spatial CellChat L-R probability: co-localizing prognostic pairs", cex.main = 0.9)
      chord_ok <- TRUE
    }
  }, error = function(e) invisible(NULL))
  if (!chord_ok) {
    plot.new()
    text(0.5, 0.5, "Chord diagram unavailable for selected prognostic pairs", cex = 1.1)
  }
  grDevices::dev.off()
}

message("Writing spatial_coloc_cellchat_lr_bubble.png ...")
if (nrow(top_dual) > 0) {
  top_pairs <- top_dual %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::mutate(pair_id = paste0(.data$source, "->", .data$target))
  lr_rows <- cc_df %>%
    dplyr::mutate(pair_id = paste0(.data$source, "->", .data$target)) %>%
    dplyr::filter(.data$pair_id %in% top_pairs$pair_id) %>%
    dplyr::group_by(.data$pair_id, .data$pathway_name) %>%
    dplyr::summarise(prob = sum(.data$prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$pair_id) %>%
    dplyr::slice_max(order_by = .data$prob, n = 3, with_ties = FALSE) %>%
    dplyr::ungroup()

  if (nrow(lr_rows) > 0) {
    p_lr <- ggplot(lr_rows, aes(x = .data$pair_id, y = .data$pathway_name, size = .data$prob, color = .data$prob)) +
      geom_point() +
      scale_color_gradient(low = "#F7F7F7", high = "#B2182B") +
      labs(
        title = "Top signaling pathways for spatially co-localizing prognostic pairs",
        x = "Sender \u2192 receiver (prognostic group)",
        y = "CellChat pathway",
        size = "Aggregated probability",
        color = "Probability"
      ) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(
      file.path(out_dir, "spatial_coloc_cellchat_lr_bubble.png"),
      p_lr,
      width = 11,
      height = 7,
      dpi = 150
    )
  }
}

message("Writing spatial_coloc_cellchat_type_spatial.png ...")
pair_type <- pair_prognostic %>%
  dplyr::filter(.data$cellchat_prob > 0) %>%
  tidyr::pivot_longer(
    cols = c(.data$cellchat_prob_physical, .data$cellchat_prob_secreted),
    names_to = "interaction_kind",
    values_to = "mode_prob"
  ) %>%
  dplyr::mutate(
    interaction_kind = dplyr::recode(
      .data$interaction_kind,
      cellchat_prob_physical = "Physical (contact / ECM)",
      cellchat_prob_secreted = "Secreted signaling"
    )
  )

p_type_spatial <- ggplot(
  pair_type,
  aes(
    x = .data$nhood_z,
    y = log1p(.data$mode_prob),
    color = .data$interaction_kind
  )
) +
  geom_point(alpha = 0.55, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c("Physical (contact / ECM)" = "#2166AC", "Secreted signaling" = "#B2182B")) +
  labs(
    title = "CellChat interaction mode vs neighborhood enrichment",
    x = "Neighborhood enrichment z-score",
    y = "log1p(mode-specific CellChat probability)",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

p_type_local <- ggplot(
  pair_type,
  aes(
    x = .data$local_cooccur,
    y = log1p(.data$mode_prob),
    color = .data$interaction_kind
  )
) +
  geom_point(alpha = 0.55, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c("Physical (contact / ECM)" = "#2166AC", "Secreted signaling" = "#B2182B")) +
  labs(
    title = "CellChat interaction mode vs local co-localization",
    x = "Local co-localization (mean non-zero score)",
    y = "log1p(mode-specific CellChat probability)",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_type_spatial_combined <- patchwork::wrap_plots(p_type_spatial, p_type_local, ncol = 2)
  ggsave(
    file.path(out_dir, "spatial_coloc_cellchat_type_spatial.png"),
    p_type_spatial_combined,
    width = 13,
    height = 6,
    dpi = 150
  )
} else {
  ggsave(
    file.path(out_dir, "spatial_coloc_cellchat_type_spatial.png"),
    p_type_spatial,
    width = 10,
    height = 6,
    dpi = 150
  )
}

message("Writing spatial_coloc_cellchat_type_pheno.png ...")
pair_pheno <- pair_prognostic %>%
  dplyr::filter(.data$cellchat_prob > 0) %>%
  tidyr::pivot_longer(
    cols = c(.data$pheno_cooccur_rho, .data$pheno_neighbor_rho),
    names_to = "pheno_metric",
    values_to = "pheno_rho"
  ) %>%
  dplyr::filter(is.finite(.data$pheno_rho)) %>%
  dplyr::mutate(
    pheno_metric = dplyr::recode(
      .data$pheno_metric,
      pheno_cooccur_rho = "PhenoMapR vs local co-localization",
      pheno_neighbor_rho = "PhenoMapR vs neighbor PhenoMapR"
    )
  )

p_type_pheno <- ggplot(
  pair_pheno,
  aes(
    x = .data$pheno_rho,
    y = .data$cellchat_frac_physical,
    color = .data$interaction_dominance
  )
) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", linewidth = 0.7) +
  facet_wrap(~ .data$pheno_metric, ncol = 2) +
  scale_color_manual(
    values = c(
      "physical-dominant" = "#2166AC",
      "secreted-dominant" = "#B2182B",
      "mixed" = "#666666",
      "none" = "#cccccc"
    ),
    name = "Dominant\ninteraction mode"
  ) +
  labs(
    title = "Physical interaction fraction vs PhenoMapR spatial correlations",
    x = "Spearman rho",
    y = "Fraction of CellChat probability from physical interactions"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "spatial_coloc_cellchat_type_pheno.png"),
  p_type_pheno,
  width = 11,
  height = 5.5,
  dpi = 150
)

message("Writing spatial_coloc_cellchat_type_assoc.png ...")
dom_summary <- pair_prognostic %>%
  dplyr::filter(.data$interaction_dominance %in% c("physical-dominant", "secreted-dominant", "mixed")) %>%
  dplyr::group_by(.data$interaction_dominance) %>%
  dplyr::summarise(
    n = dplyr::n(),
    nhood_z = mean(.data$nhood_z, na.rm = TRUE),
    local_cooccur = mean(.data$local_cooccur, na.rm = TRUE),
    pheno_neighbor_rho = mean(.data$pheno_neighbor_rho, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(.data$nhood_z, .data$local_cooccur, .data$pheno_neighbor_rho),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      .data$metric,
      nhood_z = "Neighborhood z",
      local_cooccur = "Local co-localization",
      pheno_neighbor_rho = "Pheno vs neighbor rho"
    )
  )

p_dom <- ggplot(dom_summary, aes(x = .data$interaction_dominance, y = .data$value, fill = .data$interaction_dominance)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ .data$metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "physical-dominant" = "#2166AC",
    "secreted-dominant" = "#B2182B",
    "mixed" = "#666666"
  )) +
  labs(
    title = "Spatial and PhenoMapR metrics by dominant CellChat interaction mode",
    x = NULL,
    y = "Mean value (prognostic pairs with CellChat support)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

assoc_plot <- assoc_tbl %>%
  dplyr::filter(is.finite(.data$spearman_rho)) %>%
  dplyr::mutate(
    comparison = factor(.data$comparison, levels = rev(.data$comparison))
  )

p_assoc <- ggplot(assoc_plot, aes(x = .data$spearman_rho, y = .data$comparison, fill = .data$spearman_rho > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#666666") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B")) +
  labs(
    title = "Spearman associations: interaction mode vs spatial / PhenoMapR metrics",
    subtitle = "Among prognostic pairs with any CellChat probability",
    x = "Spearman rho",
    y = NULL
  ) +
  theme_minimal(base_size = 11)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_assoc_combined <- patchwork::wrap_plots(p_dom, p_assoc, ncol = 1, heights = c(1.1, 1.4))
  ggsave(
    file.path(out_dir, "spatial_coloc_cellchat_type_assoc.png"),
    p_assoc_combined,
    width = 12,
    height = 11,
    dpi = 150
  )
} else {
  ggsave(
    file.path(out_dir, "spatial_coloc_cellchat_type_assoc.png"),
    p_assoc,
    width = 10,
    height = 7,
    dpi = 150
  )
}

integrated_hm_script <- file.path(root, "scripts", "render_spatial_integrated_heatmaps.R")
if (file.exists(integrated_hm_script)) {
  message("Rendering integrated Spearman + CellChat heatmaps ...")
  status <- system2("Rscript", c(integrated_hm_script, root), stdout = TRUE, stderr = TRUE)
  if (length(status)) cat(paste(status, collapse = "\n"), "\n")
}

expr_summary_script <- file.path(root, "scripts", "render_spatial_cellchat_expression_summary.R")
if (file.exists(expr_summary_script)) {
  message("Rendering expression-axis summary figures ...")
  status <- system2("Rscript", c(expr_summary_script, root), stdout = TRUE, stderr = TRUE)
  if (length(status)) cat(paste(status, collapse = "\n"), "\n")
}

loc_script <- file.path(root, "scripts", "render_spatial_cytospace_location_plots.R")
if (file.exists(loc_script)) {
  message("Rendering CytoSPACE location overview maps ...")
  status <- system2("Rscript", c(loc_script, root), stdout = TRUE, stderr = TRUE)
  if (length(status)) cat(paste(status, collapse = "\n"), "\n")
}

pair_maps_script <- file.path(root, "scripts", "render_spatial_pair_spatial_maps.R")
if (file.exists(pair_maps_script)) {
  message("Rendering spatial pair evidence maps ...")
  status <- system2("Rscript", c(pair_maps_script, root), stdout = TRUE, stderr = TRUE)
  if (length(status)) cat(paste(status, collapse = "\n"), "\n")
}

message(
  "Done. ", nrow(pair_prognostic), " prognostic pairs; ",
  sum(pair_prognostic$dual_spatial_lr), " with spatial co-localization + CellChat L-R; ",
  sum(pair_prognostic$interaction_dominance == "physical-dominant", na.rm = TRUE), " physical-dominant; ",
  sum(pair_prognostic$interaction_dominance == "secreted-dominant", na.rm = TRUE), " secreted-dominant."
)
message("Re-knit vignettes/spatial-transcriptomics.Rmd to pick up the figures.")
