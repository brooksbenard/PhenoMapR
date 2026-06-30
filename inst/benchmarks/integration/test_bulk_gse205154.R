# Bulk RNA-seq integration: GSE205154 (PAAD, n=289)

root <- .integration_state$config$root
dl <- integration_get_option("download")

expr_path <- integration_ensure_path("gse205154_expr", download = dl, root = root)
info_path <- integration_ensure_path("gse205154_info", download = dl, root = root)

integration_run_test(
  "bulk:gse205154:paths",
  "GSE205154 expression and phenotype files available",
  skip_if = if (is.na(expr_path) || is.na(info_path)) {
    "GSE205154 files not found (place under vignettes/ or enable download)"
  } else NULL,
  expr = {
    list(expr = expr_path, info = info_path)
  }
)

if (!is.na(expr_path) && !is.na(info_path)) {
  bulk_mat <- readRDS(expr_path)
  pheno <- readRDS(info_path)

  integration_run_test(
    "bulk:gse205154:phenomap",
    "PhenoMap PRECOG Pancreatic scores are finite",
    expr = {
      scores <- PhenoMap(
        expression = bulk_mat,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      integration_assert(is.data.frame(scores), "scores not a data.frame")
      sc <- scores[[grep("Pancreatic", colnames(scores), value = TRUE)[1]]]
      integration_assert(all(is.finite(sc)), "non-finite scores")
      integration_assert(length(sc) == nrow(pheno) || length(sc) == ncol(bulk_mat),
                         "score length mismatch")
      golden_path <- integration_fixture_path("gse205154_score_summary.rds")
      summary <- list(
        n = length(sc),
        mean = mean(sc, na.rm = TRUE),
        sd = stats::sd(sc, na.rm = TRUE),
        min = min(sc, na.rm = TRUE),
        max = max(sc, na.rm = TRUE)
      )
      if (integration_get_option("update_fixtures") || !file.exists(golden_path)) {
        saveRDS(summary, golden_path)
      } else {
        expected <- readRDS(golden_path)
        integration_assert(
          abs(summary$mean - expected$mean) < 1e-6,
          sprintf("mean score drift: %.6f vs %.6f", summary$mean, expected$mean)
        )
      }
      list(n_samples = summary$n, mean_score = summary$mean)
    }
  )

  integration_run_test(
    "bulk:gse205154:survival",
    "Primary PAAD score tertiles separate overall survival",
    expr = {
      if (!requireNamespace("survival", quietly = TRUE)) {
        stop("survival package required")
      }
      scores <- PhenoMap(
        expression = bulk_mat,
        reference = "precog",
        cancer_type = "Pancreatic",
        verbose = FALSE
      )
      col_pan <- grep("Pancreatic$", colnames(scores), value = TRUE)[1]
      scores_df <- data.frame(
        sample_id = rownames(scores),
        score = scores[[col_pan]],
        stringsAsFactors = FALSE
      )
      dat <- merge(pheno, scores_df, by = "sample_id")
      integration_assert(nrow(dat) > 50, "too few merged samples")
      qs <- stats::quantile(dat$score, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
      dat$tertile <- ifelse(dat$score <= qs[1], "low",
                            ifelse(dat$score >= qs[2], "high", "mid"))
      dat_km <- dat[dat$tertile %in% c("low", "high"), , drop = FALSE]
      fit <- survival::survdiff(
        survival::Surv(survival_time, survival_event) ~ tertile,
        data = dat_km
      )
      pval <- stats::pchisq(fit$chisq, df = 1, lower.tail = FALSE)
      integration_assert(pval < 0.05, sprintf("log-rank p = %.4f (expected < 0.05)", pval))
      list(logrank_p = pval, n = nrow(dat_km))
    }
  )
}
