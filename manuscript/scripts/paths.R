# Path helpers for manuscript/ analyses (repo root = installable package).

phenomapr_pkg_root <- function(root_arg = NULL) {
  env <- Sys.getenv("PHENOMAPR_REPO", unset = Sys.getenv("PHENOMAPR_PKG", unset = ""))
  if (nzchar(env) && file.exists(file.path(env, "DESCRIPTION"))) {
    return(normalizePath(env, winslash = "/"))
  }
  if (!is.null(root_arg) && nzchar(root_arg)) {
    ra <- normalizePath(root_arg, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(ra, "DESCRIPTION"))) return(ra)
  }
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(script_path) == 1L && nzchar(script_path)) {
    sp <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
    if (basename(sp) == "scripts" && file.exists(file.path(sp, "..", "..", "DESCRIPTION"))) {
      return(normalizePath(file.path(sp, "..", ".."), winslash = "/"))
    }
    if (grepl("/manuscript/", sp) && file.exists(file.path(sp, "..", "..", "DESCRIPTION"))) {
      return(normalizePath(file.path(sp, "..", ".."), winslash = "/"))
    }
  }
  wd <- normalizePath(getwd(), winslash = "/")
  if (file.exists(file.path(wd, "DESCRIPTION"))) return(wd)
  stop("Cannot locate PhenoMapR package root (set PHENOMAPR_REPO).")
}

manuscript_root <- function(root_arg = NULL) {
  env <- Sys.getenv("PHENOMAPR_MANUSCRIPT", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/"))
  }
  file.path(phenomapr_pkg_root(root_arg = root_arg), "manuscript")
}

integration_pkg_root <- phenomapr_pkg_root
integration_manuscript_root <- manuscript_root
integration_repo_root <- phenomapr_pkg_root

integration_vignettes_dir <- function(root_arg = NULL) {
  env <- Sys.getenv("PHENOMAPR_VIGNETTE_DATA", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/"))
  }
  file.path(phenomapr_pkg_root(root_arg = root_arg), "vignettes")
}

manuscript_data_dir <- function(root_arg = NULL) {
  file.path(manuscript_root(root_arg = root_arg), "data")
}

manuscript_benchmarks_dir <- function(root_arg = NULL) {
  file.path(manuscript_root(root_arg = root_arg), "benchmarks")
}

manuscript_results_dir <- function(root_arg = NULL) {
  file.path(manuscript_root(root_arg = root_arg), "results")
}

tcga_data_dir <- function(root_arg = NULL) {
  file.path(phenomapr_pkg_root(root_arg = root_arg), "data", "tcga")
}

load_phenomapr <- function(pkg_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_root, quiet = TRUE)
  } else {
    library(PhenoMapR)
  }
}
