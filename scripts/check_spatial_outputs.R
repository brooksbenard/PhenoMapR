#!/usr/bin/env Rscript
# Audit expected spatial pre-render outputs for vignettes/spatial-transcriptomics.Rmd

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[1L]) else normalizePath(".")

fig_dir <- file.path(root, "inst", "figures")
pair_dir <- file.path(fig_dir, "spatial_pair_maps")
data_dir <- file.path(root, "inst", "data")

expected_figs <- c(
  # CytoSPACE location overviews
  "spatial_cytospace_locations_scaled.png",
  "spatial_cytospace_locations_uniform.png",
  "spatial_cytospace_celltype_scaled.png",
  "spatial_cytospace_celltype_uniform.png",
  "spatial_cytospace_phenomapr_scaled.png",
  "spatial_cytospace_phenomapr_uniform.png",
  "spatial_cytospace_prognostic_tails_scaled.png",
  "spatial_cytospace_prognostic_tails_uniform.png",
  # Co-localization heatmaps (base)
  "spatial_colocalization_nhood_enrichment.png",
  "spatial_colocalization_nhood_enrichment_clustered.png",
  "spatial_colocalization_colocalization_scores.png",
  "spatial_colocalization_colocalization_scores_clustered.png",
  "spatial_colocalization_colocalization_scores_clustered_outlined.png",
  "spatial_colocalization_pheno_vs_cooccur_spearman.png",
  "spatial_colocalization_pheno_vs_cooccur_spearman_clustered.png",
  "spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined.png",
  "spatial_colocalization_pheno_vs_neighbor_pheno_spearman.png",
  "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered.png",
  "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined.png",
  # Integrated / CellChat
  "spatial_colocalization_pheno_vs_cooccur_spearman_clustered_outlined_integrated.png",
  "spatial_colocalization_pheno_vs_neighbor_pheno_spearman_clustered_outlined_integrated.png",
  "spatial_coloc_cellchat_bubble.png",
  "spatial_coloc_cellchat_dual_heatmap.png",
  "spatial_coloc_cellchat_chord.png",
  "spatial_coloc_cellchat_lr_bubble.png",
  "spatial_coloc_cellchat_type_spatial.png",
  "spatial_coloc_cellchat_type_pheno.png",
  "spatial_coloc_cellchat_type_assoc.png",
  "spatial_coloc_integrated_four_panel.png",
  "spatial_coloc_cellchat_prob_heatmap_clustered.png",
  "spatial_coloc_integration_evidence.png",
  # Expression-axis summaries
  "spatial_coloc_lr_topic_composition.png",
  "spatial_coloc_pathway_topic_composition.png",
  "spatial_coloc_expression_axis_heatmap.png",
  "spatial_coloc_top_lr_spatial_axes.png",
  "spatial_coloc_cellchat_spatial_network.png",
  "spatial_coloc_cellchat_spatial_feature.png"
)

expected_data <- c(
  "spatial_coloc_cellchat_pairs.rds",
  "spatial_coloc_cellchat_lr_pairs.rds",
  "spatial_coloc_cellchat_assoc.rds",
  "spatial_coloc_label_ncells.rds",
  "spatial_pair_cell_metrics.rds"
)

check <- function(paths, base) {
  data.frame(
    file = paths,
    exists = file.exists(file.path(base, paths)),
    stringsAsFactors = FALSE
  )
}

fig_chk <- check(expected_figs, fig_dir)
data_chk <- check(expected_data, data_dir)

cat("=== Spatial figure audit ===\n")
cat("Present:", sum(fig_chk$exists), "/", nrow(fig_chk), "\n")
if (any(!fig_chk$exists)) {
  cat("\nMissing figures:\n")
  writeLines(paste0("  - ", fig_chk$file[!fig_chk$exists]))
}

cat("\n=== Spatial data audit ===\n")
cat("Present:", sum(data_chk$exists), "/", nrow(data_chk), "\n")
if (any(!data_chk$exists)) {
  cat("\nMissing data:\n")
  writeLines(paste0("  - ", data_chk$file[!data_chk$exists]))
}

# Pair maps: top 5 dual-positive cross-type pairs from pairs RDS
pairs_rds <- file.path(data_dir, "spatial_coloc_cellchat_pairs.rds")
if (file.exists(pairs_rds)) {
  pair_tbl <- readRDS(pairs_rds)
  is_dual <- pair_tbl$dual_spatial_lr %in% TRUE
  is_cross <- pair_tbl$source_ct != pair_tbl$target_ct
  pairs_to_plot <- pair_tbl[is_dual & is_cross, , drop = FALSE]
  pairs_to_plot <- pairs_to_plot[order(-pairs_to_plot$integrated_score), , drop = FALSE]
  if (nrow(pairs_to_plot) > 5L) pairs_to_plot <- pairs_to_plot[seq_len(5L), , drop = FALSE]
  slug <- function(s, t) {
    clean <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
    paste0("spatial_pair_", clean(s), "__", clean(t))
  }
  pair_variants <- c("", "_uniform", "_compact", "_compact_uniform", "_tails", "_tails_uniform")
  expected_pairs <- as.character(unlist(lapply(seq_len(nrow(pairs_to_plot)), function(i) {
    pr <- pairs_to_plot[i, ]
    base <- slug(pr$source, pr$target)
    paste0(base, pair_variants, ".png")
  })))
  pair_chk <- check(expected_pairs, pair_dir)
  cat("\n=== Pair map audit (top 5 dual-positive cross-type) ===\n")
  cat("Present:", sum(pair_chk$exists), "/", nrow(pair_chk), "\n")
  if (any(!pair_chk$exists)) {
    cat("\nMissing pair maps:\n")
    writeLines(paste0("  - ", pair_chk$file[!pair_chk$exists]))
  }
}

invisible(if (any(!fig_chk$exists) || any(!data_chk$exists)) quit(status = 1L))
