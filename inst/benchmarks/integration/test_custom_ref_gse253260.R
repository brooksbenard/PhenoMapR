# Custom reference derivation: GSE253260 (PDAC bulk, survival)

root <- .integration_state$config$root
dl <- integration_get_option("download")

expr_path <- integration_ensure_path("gse253260_expr", download = dl, root = root)
info_path <- integration_ensure_path("gse253260_info", download = dl, root = root)

integration_run_test(
  "custom:gse253260:paths",
  "GSE253260 expression and phenotype files available",
  skip_if = if (is.na(expr_path) || is.na(info_path)) {
    "GSE253260 files not found"
  } else NULL,
  expr = list(expr = expr_path, info = info_path)
)

if (!is.na(expr_path) && !is.na(info_path)) {
  bulk_mat <- readRDS(expr_path)
  pheno <- readRDS(info_path)
  pheno <- pheno[pheno$geo_accession %in% rownames(bulk_mat), , drop = FALSE]

  integration_run_test(
    "custom:gse253260:derive",
    "derive_reference_from_bulk survival signature",
    expr = {
      ref <- suppressWarnings(derive_reference_from_bulk(
        bulk_expression = t(bulk_mat),
        phenotype = pheno,
        sample_id_column = "geo_accession",
        phenotype_type = "survival",
        survival_time = "OS_Time",
        survival_event = "OS_Censor",
        verbose = FALSE
      ))
      integration_assert(is.data.frame(ref), "reference not a data.frame")
      integration_assert(nrow(ref) >= 100, "too few signature genes")
      zcol <- grep("z_score|weighted|survival", colnames(ref), ignore.case = TRUE, value = TRUE)[1]
      integration_assert(!is.na(zcol), "no z-score column in reference")
      z <- ref[[zcol]]
      integration_assert(any(abs(z) >= 2), "expected |z| >= 2 genes in signature")
      list(n_genes = nrow(ref), n_sig = sum(abs(z) >= 2))
    }
  )

  integration_run_test(
    "custom:gse253260:score",
    "Score GSE253260 with derived custom reference",
    expr = {
      ref <- suppressWarnings(derive_reference_from_bulk(
        bulk_expression = t(bulk_mat),
        phenotype = pheno,
        sample_id_column = "geo_accession",
        phenotype_type = "survival",
        survival_time = "OS_Time",
        survival_event = "OS_Censor",
        verbose = FALSE
      ))
      scores <- PhenoMap(
        expression = t(bulk_mat),
        reference = ref,
        verbose = FALSE
      )
      sc <- scores[[1]]
      integration_assert(all(is.finite(sc)), "non-finite custom scores")
      if (requireNamespace("survival", quietly = TRUE)) {
        dat <- data.frame(
          time = pheno$OS_Time[match(rownames(scores), pheno$geo_accession)],
          event = pheno$OS_Censor[match(rownames(scores), pheno$geo_accession)],
          score = sc
        )
        dat <- dat[complete.cases(dat), , drop = FALSE]
        med <- stats::median(dat$score, na.rm = TRUE)
        dat$grp <- ifelse(dat$score >= med, "high", "low")
        fit <- survival::survdiff(
          survival::Surv(time, event) ~ grp,
          data = dat
        )
        pval <- stats::pchisq(fit$chisq, df = 1, lower.tail = FALSE)
        integration_assert(pval < 0.1, sprintf("median-split survival p = %.4f", pval))
        list(survival_p = pval)
      } else {
        list(n = length(sc))
      }
    }
  )
}
