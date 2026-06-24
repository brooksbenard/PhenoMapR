#!/usr/bin/env Rscript
# Pre-render spatial co-localization heatmap for the spatial-transcriptomics vignette.
# Uses the full CytoSPACE cell set (no subsampling) and spatialCooccur::nhood_enrichment().
#
# Run from package root:
#   Rscript scripts/render_spatial_colocalization_heatmap.R
#   Rscript scripts/render_spatial_colocalization_heatmap.R /path/to/PhenoMapR
#
# Output: inst/figures/spatial_colocalization_nhood_enrichment.png

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
  library(ComplexHeatmap)
  library(circlize)
})

if (!requireNamespace("spatialCooccur", quietly = TRUE)) {
  stop(
    "Install spatialCooccur: remotes::install_github(\"juninamo/spatialCooccur\")",
    call. = FALSE
  )
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

row_ha <- ComplexHeatmap::rowAnnotation(
  CellType = row_ct,
  `Prognostic Group` = row_pg,
  col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
  show_annotation_name = FALSE,
  simple_anno_size = grid::unit(4, "mm"),
  show_legend = c(FALSE)
)
col_ha <- ComplexHeatmap::HeatmapAnnotation(
  `# cells` = ComplexHeatmap::anno_barplot(
    col_ncells,
    gp = grid::gpar(fill = "#666666"),
    border = FALSE,
    height = grid::unit(14, "mm"),
    ylim = c(0, max(col_ncells, 1L, na.rm = TRUE))
  ),
  CellType = col_ct,
  `Prognostic Group` = col_pg,
  col = list(CellType = pal_ct, `Prognostic Group` = pal_pg),
  annotation_name_side = "right",
  simple_anno_size = grid::unit(4, "mm")
)

ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Z",
  col = col_fun,
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
  heatmap_legend_param = list(title = "Z-score"),
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
grDevices::png(out_png, width = 7, height = 6, units = "in", res = 150)
ComplexHeatmap::draw(
  ht,
  column_title = "Neighborhood Enrichment of Prognostic Cell Types",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grDevices::dev.off()

message("Done. Re-knit vignettes/spatial-transcriptomics.Rmd to pick up the figure.")
