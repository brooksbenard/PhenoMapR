# Performance / memory stress tier (--stress)

if (!isTRUE(integration_get_option("stress"))) {
  integration_record(
    "stress:skipped",
    "Performance tier not requested (pass --stress)",
    "skip",
    list(reason = "use --stress to enable")
  )
} else {
  root <- .integration_state$config$root
  baseline_path <- integration_fixture_path("stress_baselines.json")

  record_baseline <- function(id, elapsed, peak_mb = NA_real_) {
    entry <- list(elapsed_sec = elapsed, peak_mb = peak_mb)
    if (integration_get_option("update_fixtures") || !file.exists(baseline_path)) {
      baselines <- if (file.exists(baseline_path)) jsonlite::fromJSON(baseline_path) else list()
      baselines[[id]] <- entry
      jsonlite::write_json(baselines, baseline_path, auto_unbox = TRUE, pretty = TRUE)
      return(entry)
    }
    baselines <- jsonlite::fromJSON(baseline_path)
    if (!id %in% names(baselines)) {
      baselines[[id]] <- entry
      jsonlite::write_json(baselines, baseline_path, auto_unbox = TRUE, pretty = TRUE)
      return(entry)
    }
    old <- baselines[[id]]
    integration_assert(
      elapsed <= old$elapsed_sec * 1.5,
      sprintf("%s slower than baseline: %.1fs vs %.1fs", id, elapsed, old$elapsed_sec)
    )
    entry
  }

  integration_run_test(
    "stress:demo_scoring",
    "Score 5000-cell demo bundle wall time",
    skip_if = if (!requireNamespace("jsonlite", quietly = TRUE)) "jsonlite required" else NULL,
    expr = {
      helpers <- integration_source_shiny_helpers()
      demo_path <- integration_ensure_path("shiny_demo", download = FALSE, root = root)
      integration_assert(!is.na(demo_path), "demo bundle missing")
      Sys.setenv(PHENOMAPR_SHINY_DEMO_RDS = demo_path)
      helpers$clear_shiny_demo_pool_cache()
      demo <- helpers$make_shiny_demo_dataset(n_genes = 500L, n_cells = 5000L, seed = 1L)
      peak <- NA_real_
      if (requireNamespace("peakRAM", quietly = TRUE)) {
        pr <- peakRAM::peakRAM({
          t0 <- proc.time()
          scores <- PhenoMap(
            expression = demo$expression,
            reference = "precog",
            cancer_type = "Pancreatic",
            verbose = FALSE
          )
          elapsed <- (proc.time() - t0)[["elapsed"]]
        })
        peak <- pr$Peak_RAM_Used_MiB
      } else {
        t0 <- proc.time()
        scores <- PhenoMap(
          expression = demo$expression,
          reference = "precog",
          cancer_type = "Pancreatic",
          verbose = FALSE
        )
        elapsed <- (proc.time() - t0)[["elapsed"]]
      }
      integration_assert(nrow(scores) == 5000L, "score row count mismatch")
      bl <- record_baseline("demo_5000_score", elapsed, peak)
      list(elapsed = elapsed, peak_mb = peak, baseline = bl)
    }
  )

  if (isTRUE(integration_get_option("full"))) {
    h5_path <- integration_resolve_path("cra001160_h5", root = root)
    integration_run_test(
      "stress:cra_h5_load",
      "Load CRA001160 10X h5 subset (full tier)",
      skip_if = if (is.na(h5_path)) "CRA001160 h5 not found" else NULL,
      expr = {
        if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat required")
        t0 <- proc.time()
        expr <- Seurat::Read10X_h5(h5_path)
        elapsed_load <- (proc.time() - t0)[["elapsed"]]
        cells <- colnames(expr)
        set.seed(1)
        keep <- sample(cells, min(10000L, length(cells)))
        sub <- expr[, keep, drop = FALSE]
        t1 <- proc.time()
        scores <- PhenoMap(
          expression = sub,
          reference = "precog",
          cancer_type = "Pancreatic",
          verbose = FALSE
        )
        elapsed_score <- (proc.time() - t1)[["elapsed"]]
        bl <- record_baseline("cra10k_score", elapsed_score)
        list(load_sec = elapsed_load, score_sec = elapsed_score, n = ncol(sub), baseline = bl)
      }
    )
  }
}
