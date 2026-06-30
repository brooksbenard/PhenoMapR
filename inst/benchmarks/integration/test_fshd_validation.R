# FSHD validation cache replay (bulk + snRNA)

root <- .integration_state$config$root
vdir <- integration_vignettes_dir(root)
fshd_paths <- integration_fshd_cache_paths(root)

integration_run_test(
  "fshd:cache:files",
  "FSHD validation RDS caches present",
  skip_if = if (all(is.na(fshd_paths))) "No FSHD_validation_*.rds files in vignettes/" else NULL,
  expr = {
    present <- names(fshd_paths)[!is.na(fshd_paths)]
    integration_assert(length(present) >= 3L, "expected at least 3 FSHD cache files")
    list(files = paste(present, collapse = ", "))
  }
)

bulk_expr_path <- integration_resolve_path("fshd_bulk_expr", root = root)
bulk_pheno_path <- integration_resolve_path("fshd_bulk_pheno", root = root)

if (!is.na(bulk_expr_path) && !is.na(bulk_pheno_path)) {
  integration_run_test(
    "fshd:ref:derive",
    "Derive FSHD training reference from GSE140261",
    expr = {
      expr_mat <- readRDS(bulk_expr_path)
      pheno <- readRDS(bulk_pheno_path)
      ref <- derive_reference_from_bulk(
        bulk_expression = expr_mat,
        phenotype = pheno,
        sample_id_column = "sample_id",
        phenotype_type = "binary",
        phenotype_column = "group",
        verbose = FALSE
      )
      integration_assert(nrow(ref) >= 50, "FSHD reference too small")
      list(n_genes = nrow(ref))
    }
  )

  integration_run_test(
    "fshd:ref:training_separation",
    "Training bulk scores separate FSHD vs Control",
    expr = {
      expr_mat <- readRDS(bulk_expr_path)
      pheno <- readRDS(bulk_pheno_path)
      ref <- derive_reference_from_bulk(
        bulk_expression = expr_mat,
        phenotype = pheno,
        sample_id_column = "sample_id",
        phenotype_type = "binary",
        phenotype_column = "group",
        verbose = FALSE
      )
      scores <- PhenoMap(expression = expr_mat, reference = ref, verbose = FALSE)
      sc <- scores[[1]]
      grp <- pheno$group[match(rownames(scores), pheno$sample_id)]
      integration_assert(length(unique(grp)) == 2L, "expected binary groups")
      wt <- stats::wilcox.test(sc ~ grp)
      integration_assert(wt$p.value < 0.05, sprintf("training Wilcoxon p = %.4f", wt$p.value))
      medians <- tapply(sc, grp, stats::median, na.rm = TRUE)
      integration_assert(
        length(unique(sign(diff(medians)))) >= 1L,
        "FSHD and Control training scores should differ"
      )
      list(p = wt$p.value)
    }
  )
}

if (!is.na(fshd_paths[["FSHD_validation_bulk_plot_df"]])) {
  integration_run_test(
    "fshd:cache:bulk_wilcox",
    "Cached bulk validation Wilcoxon matches recomputed from plot df",
    expr = {
      plot_df <- readRDS(fshd_paths[["FSHD_validation_bulk_plot_df"]])
      cached_w <- readRDS(fshd_paths[["FSHD_validation_wilcox"]])
      cohorts <- unique(plot_df$cohort)
      pvals <- vapply(cohorts, function(coh) {
        sub <- plot_df[plot_df$cohort == coh, , drop = FALSE]
        wt <- stats::wilcox.test(sub$score ~ sub$group)
        wt$p.value
      }, numeric(1))
      names(pvals) <- cohorts
      for (coh in cohorts) {
        if (!coh %in% names(cached_w)) next
        cached_p <- cached_w[[coh]]$p.value
        integration_assert(
          abs(pvals[[coh]] - cached_p) < 1e-4,
          sprintf("%s Wilcoxon p mismatch: %.6f vs %.6f", coh, pvals[[coh]], cached_p)
        )
      }
      list(cohorts = paste(cohorts, collapse = ", "))
    }
  )
}

if (!is.na(fshd_paths[["FSHD_validation_snRNA_scores"]])) {
  integration_run_test(
    "fshd:cache:snrna",
    "Cached snRNA Day5 FSHD2 scores exceed Control",
    expr = {
      sn <- readRDS(fshd_paths[["FSHD_validation_snRNA_scores"]])
      d5 <- sn[sn$day == "Day5", , drop = FALSE]
      integration_assert(nrow(d5) > 100, "too few Day5 nuclei")
      wt <- stats::wilcox.test(d5$score ~ d5$status)
      integration_assert(wt$p.value < 0.001, sprintf("Day5 snRNA p = %.2e", wt$p.value))
      list(p = wt$p.value, n = nrow(d5))
    }
  )
}
