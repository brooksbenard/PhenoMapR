# Dataset path resolution for integration stress tests.
# Sourced by run_integration_stress.R after repo root is known.

.integration_drive_ids <- list(
  cra001160_h5   = "1PolTXggREz8XmhutCLTQJGCfKxFAzqMl",
  cra001160_meta = "17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu",
  spatial_cyto   = "1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll",
  gse205154_expr = "16P-rfXD734seGl_xxTQuSDGMStQk9G7B",
  gse205154_info = "1c0GUMVzAr6K44Z7CfVRz-dy6jIvb0nNo",
  gse253260_expr = "1NXLzKx0O6q-9Bnl_qPwA0hhijIB9jk-6",
  gse253260_info = "1j0syd8soOQNkudhGnc407TVCJuuvm4jk"
)

.integration_env_urls <- list(
  cra001160_h5   = "PHENOMAPR_CRA001160_H5_URL",
  cra001160_meta = "PHENOMAPR_CRA001160_META_URL",
  spatial_cyto   = "PHENOMAPR_SPATIAL_RDS_URL",
  gse205154_expr = "PHENOMAPR_GSE205154_MATRIX_URL",
  gse205154_info = "PHENOMAPR_GSE205154_INFO_URL",
  gse253260_expr = "PHENOMAPR_GSE253260_EXPR_URL",
  gse253260_info = "PHENOMAPR_GSE253260_INFO_URL"
)

.integration_default_filenames <- list(
  cra001160_h5   = "PAAD_CRA001160_expression.h5",
  cra001160_meta = "PAAD_CRA001160_CellMetainfo_table.tsv",
  cra001160_rds  = "PAAD_CRA001160_seurat.rds",
  spatial_cyto   = "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds",
  spatial_spot   = "HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds",
  gse205154_expr = "GSE205154_GPL20301_expression.rds",
  gse205154_info = "GSE205154_info.rds",
  gse253260_expr = "GSE253260_expression.rds",
  gse253260_info = "GSE253260_info.rds",
  shiny_demo     = "PAAD_CRA001160_demo_5000.rds",
  cosmx_h5ad     = "seurat_test_core_46.h5ad",
  fshd_bulk_expr = "FSHD_GSE140261_bulk_expression.rds",
  fshd_bulk_pheno = "FSHD_GSE140261_phenotype.rds"
)

integration_repo_root <- function() {
  root <- Sys.getenv("PHENOMAPR_REPO", unset = "")
  if (nzchar(root) && file.exists(file.path(root, "DESCRIPTION"))) {
    return(normalizePath(root, winslash = "/"))
  }
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(script_path) == 1L && nzchar(script_path)) {
    cand <- normalizePath(
      file.path(dirname(script_path), "..", "..", ".."),
      mustWork = FALSE
    )
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
  }
  wd <- normalizePath(getwd(), winslash = "/")
  if (file.exists(file.path(wd, "DESCRIPTION"))) return(wd)
  stop("Cannot locate PhenoMapR package root (set PHENOMAPR_REPO).")
}

integration_vignettes_dir <- function(root = integration_repo_root()) {
  file.path(root, "vignettes")
}

integration_first_existing <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  for (p in paths) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NA_character_
}

integration_resolve_path <- function(key, root = integration_repo_root()) {
  defaults <- .integration_default_filenames
  if (!key %in% names(defaults)) {
    stop("Unknown integration dataset key: ", key)
  }
  fname <- defaults[[key]]
  vdir <- integration_vignettes_dir(root)
  candidates <- c(
    Sys.getenv(paste0("PHENOMAPR_", toupper(gsub("[^a-z0-9]", "_", key))), unset = ""),
    file.path(vdir, fname),
    file.path(root, fname),
    file.path(root, "inst", "extdata", "shiny", fname)
  )
  if (key == "cosmx_h5ad") {
    candidates <- c(
      Sys.getenv("PHENOMAPR_COSMX_H5AD", unset = ""),
      file.path(path.expand("~/Downloads"), fname),
      candidates
    )
  }
  if (key == "shiny_demo") {
    candidates <- c(
      Sys.getenv("PHENOMAPR_SHINY_DEMO_RDS", unset = ""),
      file.path(root, "inst", "extdata", "shiny", fname),
      candidates
    )
  }
  if (grepl("^fshd_", key)) {
    candidates <- c(file.path(vdir, fname), candidates)
  }
  if (grepl("^gse", key)) {
    candidates <- c(
      file.path(root, "inst", "extdata", fname),
      candidates
    )
  }
  integration_first_existing(candidates)
}

integration_download_dataset <- function(key, root = integration_repo_root()) {
  dest <- integration_resolve_path(key, root = root)
  if (!is.na(dest)) return(dest)
  fname <- .integration_default_filenames[[key]]
  vdir <- integration_vignettes_dir(root)
  dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(vdir, fname)
  env_var <- .integration_env_urls[[key]]
  if (!is.null(env_var)) {
    url <- Sys.getenv(env_var, unset = "")
    if (nzchar(url)) {
      message("Downloading ", fname, " from ", env_var, " ...")
      utils::download.file(url, dest, mode = "wb", quiet = TRUE)
      if (file.exists(dest)) return(normalizePath(dest, winslash = "/"))
    }
  }
  gd_id <- .integration_drive_ids[[key]]
  if (is.null(gd_id)) return(NA_character_)
  if (!requireNamespace("googledrive", quietly = TRUE)) {
    warning("googledrive not installed; cannot download ", fname)
    return(NA_character_)
  }
  options(googledrive_quiet = TRUE)
  googledrive::drive_deauth()
  message("Downloading ", fname, " from Google Drive ...")
  tryCatch(
    googledrive::drive_download(googledrive::as_id(gd_id), dest, overwrite = TRUE),
    error = function(e) {
      warning("Drive download failed for ", fname, ": ", conditionMessage(e))
    }
  )
  if (file.exists(dest)) normalizePath(dest, winslash = "/") else NA_character_
}

integration_ensure_path <- function(key, download = TRUE, root = integration_repo_root()) {
  path <- integration_resolve_path(key, root = root)
  if (!is.na(path)) return(path)
  if (!isTRUE(download)) return(NA_character_)
  integration_download_dataset(key, root = root)
}

integration_load_config <- function(root = integration_repo_root()) {
  cfg <- list(
    root = root,
    vignettes = integration_vignettes_dir(root),
    fixtures = file.path(root, "inst", "benchmarks", "integration", "fixtures", "expected"),
    reports = file.path(root, "inst", "benchmarks", "integration", "reports"),
    shiny_dir = file.path(root, "inst", "shiny"),
    local_paths = list()
  )
  local_r <- file.path(root, "inst", "benchmarks", "integration", "local_paths.R")
  if (file.exists(local_r)) {
    loc <- new.env(parent = emptyenv())
    sys.source(local_r, envir = loc)
    if (exists("local_datasets", envir = loc)) {
      cfg$local_paths <- loc$local_datasets
    }
  }
  cfg
}

integration_fshd_cache_paths <- function(root = integration_repo_root()) {
  vdir <- integration_vignettes_dir(root)
  patterns <- c(
    "FSHD_validation_bulk_plot_df.rds",
    "FSHD_validation_wilcox.rds",
    "FSHD_validation_summary_table.rds",
    "FSHD_validation_snRNA_scores.rds"
  )
  setNames(
    vapply(patterns, function(p) {
      fp <- file.path(vdir, p)
      if (file.exists(fp)) normalizePath(fp, winslash = "/") else NA_character_
    }, character(1)),
    sub("\\.rds$", "", patterns)
  )
}
