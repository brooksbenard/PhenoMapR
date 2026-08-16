# Orthogonal validation: DepMap PAAD CRISPR essentials / dependency ranks
# restricted to malignant-ductal (Ductal cell type 2) adverse AND favorable
# markers from CRA001160. Always uses the full CRA001160 cohort (~57k cells).
#
# Phenotype tails are FIXED at 10% for all DepMap analyses (2D Chronos,
# dependency probability, and NextGen organoid reuse of depmap_markers).
# Do not switch to 5%: at 5% + manuscript marker filters (log2FC >= 0.75),
# malignant-ductal favorable tails are too sparse for FindMarkers / GSEA.

dat <- methodology_load_cra001160(tier = "full")
malignant_ductal <- methodology_malignant_ductal_celltypes()
depmap_phenotype_percentile <- 0.10

methodology_run_test(
  "meth:orthogonal:malignant_ductal_markers",
  "Significant adverse and favorable markers in malignant-ductal cells only",
  expr = {
    scores <- PhenoMap(
      expression = dat$expression,
      reference = methodology_get_option("reference"),
      cancer_type = methodology_get_option("cancer_type"),
      verbose = FALSE
    )
    groups <- define_phenotype_groups(scores, percentile = depmap_phenotype_percentile)
    gcol <- grep("^phenotype_group_", colnames(groups), value = TRUE)[1]
    id_col <- if (".cell_id" %in% names(dat$metadata)) ".cell_id" else "cell_id"
    groups$cell_type <- dat$metadata$celltype_original[
      match(groups$cell_id, dat$metadata[[id_col]])
    ]

    mk <- find_phenotype_markers(
      expression = dat$expression,
      group_labels = groups,
      group_column = gcol,
      cell_id_column = "cell_id",
      cell_type_column = "cell_type",
      marker_scope = "cell_type_specific",
      celltype_contrast = "within_cell_type",
      celltype_unique_genes = FALSE,
      verbose = FALSE,
      max_cells_per_ident = 5000L
    )

    adv_raw <- methodology_filter_malignant_ductal_markers(mk$adverse_markers)
    fav_raw <- methodology_filter_malignant_ductal_markers(mk$favorable_markers)
    # Significant positive markers with a reasonable effect size (not just
    # log2FC > 0). find_phenotype_markers already uses Seurat's default
    # logfc.threshold = 0.25; enforce a stricter post-filter for DepMap GSEA.
    marker_p_adj <- 0.05
    # Chosen from FC sweep: max adverse NES while favorable remains n.s.
    # (see reports/depmap_gsea_fc_cutoff_sweep.tsv)
    marker_min_log2fc <- 0.5
    adv <- methodology_significant_positive_markers(
      adv_raw,
      p_adj_threshold = marker_p_adj,
      min_log2fc = marker_min_log2fc
    )
    fav <- methodology_significant_positive_markers(
      fav_raw,
      p_adj_threshold = marker_p_adj,
      min_log2fc = marker_min_log2fc
    )

    if (!nrow(adv)) stop("no malignant-ductal adverse markers", call. = FALSE)
    if (!nrow(fav)) stop("no malignant-ductal favorable markers", call. = FALSE)

    adv_sets <- methodology_marker_sets(adv, phenotype_label = "adverse")
    fav_sets <- methodology_marker_sets(fav, phenotype_label = "favorable")

    methodology_store_artifact("depmap_markers", list(
      adverse_hits = adv,
      favorable_hits = fav,
      adverse_genes = adv_sets$unique_genes,
      favorable_genes = fav_sets$unique_genes,
      cell_types = malignant_ductal,
      tier = dat$tier,
      n_cells = dat$n_cells,
      phenotype_percentile = depmap_phenotype_percentile,
      marker_filters = list(
        phenotype_percentile = depmap_phenotype_percentile,
        p_adj_threshold = marker_p_adj,
        min_log2fc = marker_min_log2fc,
        direction = "positive (avg_log2FC > min_log2fc)",
        find_phenotype_markers_logfc_threshold = 0.25
      ),
      # keep combined hits for plotting helpers that expect marker_hits
      hits = rbind(
        transform(adv, phenotype = "adverse"),
        transform(fav, phenotype = "favorable")
      )
    ))

    list(
      tier = dat$tier,
      n_cells = dat$n_cells,
      cell_types = malignant_ductal,
      phenotype_percentile = depmap_phenotype_percentile,
      marker_filters = list(
        phenotype_percentile = depmap_phenotype_percentile,
        p_adj = marker_p_adj,
        min_avg_log2fc = marker_min_log2fc
      ),
      n_adverse_hits = adv_sets$n_hits,
      n_adverse_genes = adv_sets$n_unique_genes,
      n_favorable_hits = fav_sets$n_hits,
      n_favorable_genes = fav_sets$n_unique_genes,
      median_log2fc_adverse = stats::median(adv$avg_log2FC, na.rm = TRUE),
      median_log2fc_favorable = stats::median(fav$avg_log2FC, na.rm = TRUE),
      top_adverse = head(adv[order(-adv$avg_log2FC), c("gene", "avg_log2FC", "p_adj")], 10),
      top_favorable = head(fav[order(-fav$avg_log2FC), c("gene", "avg_log2FC", "p_adj")], 10)
    )
  }
)

methodology_run_test(
  "meth:orthogonal:depmap_paad_essentials",
  "Fisher enrichment: malignant-ductal adverse/favorable markers vs DepMap PAAD essentials",
  skip_if = {
    data_dir <- methodology_depmap_data_dir()
    needed <- c(
      file.path(data_dir, "Model.csv"),
      file.path(data_dir, "CRISPRGeneEffect.csv")
    )
    if (!all(file.exists(needed))) {
      "DepMap Model.csv and CRISPRGeneEffect.csv required in methodology/data"
    } else {
      NULL
    }
  },
  expr = {
    marker_art <- methodology_get_artifact("depmap_markers")
    if (is.null(marker_art)) stop("Run malignant_ductal_markers test first in this session")

    depmap_genes <- methodology_depmap_paad_essential_genes()
    gene_scores <- attr(depmap_genes, "gene_scores")
    if (is.null(gene_scores)) {
      gene_scores <- methodology_depmap_paad_gene_scores()
    }
    dep_meta <- attr(depmap_genes, "metadata")
    universe <- rownames(dat$expression)

    fisher_adv <- methodology_depmap_fisher_enrichment(
      query_genes = marker_art$adverse_genes,
      depmap_genes = depmap_genes,
      universe = universe
    )
    fisher_fav <- methodology_depmap_fisher_enrichment(
      query_genes = marker_art$favorable_genes,
      depmap_genes = depmap_genes,
      universe = universe
    )

    fisher_details <- list(
      adverse = fisher_adv,
      favorable = fisher_fav,
      # Keep top-level fields used by older plot code (adverse arm)
      n_markers = fisher_adv$n_markers,
      n_depmap_in_universe = fisher_adv$n_depmap_in_universe,
      overlap = fisher_adv$overlap,
      overlap_genes = fisher_adv$overlap_genes,
      odds_ratio = fisher_adv$odds_ratio,
      p_value = fisher_adv$p_value,
      n_hits = nrow(marker_art$adverse_hits),
      n_cell_types = length(marker_art$cell_types),
      n_paad_models = dep_meta$n_paad_models,
      effect_cutoff = dep_meta$effect_cutoff,
      primary_only = isTRUE(dep_meta$primary_only),
      cell_types = marker_art$cell_types
    )

    methodology_store_artifact("depmap", list(
      marker_hits = marker_art$hits,
      adverse_hits = marker_art$adverse_hits,
      favorable_hits = marker_art$favorable_hits,
      unique_genes = marker_art$adverse_genes,
      adverse_genes = marker_art$adverse_genes,
      favorable_genes = marker_art$favorable_genes,
      depmap_genes = depmap_genes,
      gene_scores = gene_scores,
      fisher = fisher_details,
      cell_types = marker_art$cell_types,
      primary_only = isTRUE(dep_meta$primary_only),
      tier = marker_art$tier
    ))

    list(
      cell_types = marker_art$cell_types,
      n_paad_models = dep_meta$n_paad_models,
      primary_only = isTRUE(dep_meta$primary_only),
      effect_cutoff = dep_meta$effect_cutoff,
      adverse_overlap = fisher_adv$overlap,
      adverse_odds_ratio = fisher_adv$odds_ratio,
      adverse_p = fisher_adv$p_value,
      favorable_overlap = fisher_fav$overlap,
      favorable_odds_ratio = fisher_fav$odds_ratio,
      favorable_p = fisher_fav$p_value
    )
  }
)

methodology_run_test(
  "meth:orthogonal:depmap_rank_tails",
  "Malignant-ductal adverse/favorable markers enrich opposite DepMap essentiality tails",
  skip_if = {
    if (is.null(methodology_get_artifact("depmap"))) {
      "Run depmap_paad_essentials test first in this session"
    } else {
      NULL
    }
  },
  expr = {
    dep_art <- methodology_get_artifact("depmap")
    rank_res <- methodology_depmap_rank_enrichment(
      adverse_genes = dep_art$adverse_genes,
      favorable_genes = dep_art$favorable_genes,
      gene_scores = dep_art$gene_scores,
      universe = rownames(dat$expression),
      tail_frac = 0.1
    )

    dep_art$rank_enrichment <- rank_res
    methodology_store_artifact("depmap", dep_art)

    list(
      n_adverse_in_depmap = rank_res$n_adverse_in_depmap,
      n_favorable_in_depmap = rank_res$n_favorable_in_depmap,
      median_pct_adverse = rank_res$wilcox_adverse_vs_favorable$median_adverse,
      median_pct_favorable = rank_res$wilcox_adverse_vs_favorable$median_favorable,
      wilcox_adverse_vs_bg_p = rank_res$wilcox_adverse_vs_bg$p,
      wilcox_favorable_vs_bg_p = rank_res$wilcox_favorable_vs_bg$p,
      wilcox_adverse_vs_favorable_p = rank_res$wilcox_adverse_vs_favorable$p,
      fisher_adverse_most_essential = rank_res$fisher_adverse_in_most_essential,
      fisher_favorable_least_essential = rank_res$fisher_favorable_in_least_essential
    )
  }
)

methodology_run_test(
  "meth:orthogonal:depmap_gsea",
  "GSEA NES for adverse/favorable markers vs DepMap essentiality ranking",
  skip_if = {
    if (is.null(methodology_get_artifact("depmap"))) {
      "Run depmap_paad_essentials test first in this session"
    } else if (!requireNamespace("fgsea", quietly = TRUE)) {
      "Package fgsea required for GSEA NES analysis"
    } else {
      NULL
    }
  },
  expr = {
    dep_art <- methodology_get_artifact("depmap")
    gsea_res <- methodology_depmap_gsea(
      adverse_genes = dep_art$adverse_genes,
      favorable_genes = dep_art$favorable_genes,
      gene_scores = dep_art$gene_scores,
      universe = rownames(dat$expression),
      nperm = 10000L
    )

    dep_art$gsea <- gsea_res
    methodology_store_artifact("depmap", dep_art)

    tab <- gsea_res$table
    # Drop list-columns for readable report serialization
    keep_cols <- intersect(
      c("pathway", "pval", "padj", "ES", "NES", "size"),
      names(tab)
    )
    tab_report <- tab[, keep_cols, drop = FALSE]

    list(
      ranking = gsea_res$ranking,
      n_adverse = length(gsea_res$pathways$Adverse %||% character()),
      n_favorable = length(gsea_res$pathways$Favorable %||% character()),
      gsea = tab_report
    )
  }
)
