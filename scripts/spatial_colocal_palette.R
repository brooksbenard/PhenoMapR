# Shared palettes for spatial co-localization and CellChat figures.

# Co-localization score (0 = cream/off-white, high = green; matches nhood positive arm).
.spatial_coloc_score_palette_colors <- function() {
  c("#FCFBE3", "#f7f7f7", "#b8ddd4", "#5a9e8f", "#005C55")
}

.spatial_coloc_score_col_fun <- function(lim) {
  if (!is.finite(lim) || lim <= 0) lim <- 1
  cols <- .spatial_coloc_score_palette_colors()
  breaks <- seq(0, lim, length.out = length(cols))
  circlize::colorRamp2(breaks, cols)
}

.spatial_coloc_score_scale_gg <- function(name = "Local\nco-localization", limits = c(0, NA)) {
  ggplot2::scale_color_gradientn(
    colours = .spatial_coloc_score_palette_colors(),
    limits = limits,
    name = name
  )
}

# CellChat probability (RdPu; used for L-R probability heatmaps).
.spatial_cellchat_prob_palette_colors <- function() {
  c(
    "#fff7f3", "#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1",
    "#dd3497", "#ae017e", "#7a0177", "#49006a"
  )
}

.spatial_cellchat_prob_col_fun <- function(lim) {
  if (!is.finite(lim) || lim <= 0) lim <- 1
  cols <- .spatial_cellchat_prob_palette_colors()
  breaks <- seq(0, lim, length.out = length(cols))
  circlize::colorRamp2(breaks, cols)
}

.spatial_cellchat_prob_scale_gg <- function(name = "CellChat\nprobability", limits = c(0, NA)) {
  ggplot2::scale_color_gradientn(
    colours = .spatial_cellchat_prob_palette_colors(),
    limits = limits,
    name = name
  )
}

# Backward-compatible aliases (CellChat RdPu).
.spatial_colocal_palette_colors <- .spatial_cellchat_prob_palette_colors
.spatial_colocal_col_fun <- .spatial_cellchat_prob_col_fun
.spatial_colocal_scale_gg <- .spatial_cellchat_prob_scale_gg
