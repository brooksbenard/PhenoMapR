#!/usr/bin/env Rscript
## PhenoMapR real-data integration stress test harness.
##
## Usage (from package root):
##   Rscript inst/benchmarks/integration/run_integration_stress.R
##   Rscript inst/benchmarks/integration/run_integration_stress.R --update-fixtures
##   Rscript inst/benchmarks/integration/run_integration_stress.R --stress --full
##
## Options:
##   --update-fixtures   Refresh golden marker/score fixtures
##   --stress            Run performance / memory tier
##   --full              Run heavy spatial / large CRA001160 tests
##   --no-download       Skip Google Drive downloads for missing files
##   --skip-smokes       Skip inst/benchmarks/*_smoke.R UI source tests

args <- commandArgs(trailingOnly = TRUE)
options <- list(
  update_fixtures = "--update-fixtures" %in% args,
  stress = "--stress" %in% args,
  full = "--full" %in% args,
  download = !("--no-download" %in% args),
  skip_smokes = "--skip-smokes" %in% args
)

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
int_dir <- if (length(script_path) == 1L && nzchar(script_path)) {
  normalizePath(dirname(script_path), mustWork = FALSE)
} else {
  normalizePath("inst/benchmarks/integration", mustWork = FALSE)
}

source(file.path(int_dir, "config.R"), local = TRUE)
source(file.path(int_dir, "helpers_integration.R"), local = TRUE)

root <- integration_repo_root()
cfg <- integration_init(options = options, root = root)
integration_load_package(root)

test_files <- sort(list.files(int_dir, pattern = "^test_.*\\.R$", full.names = TRUE))
for (tf in test_files) {
  message("\n=== ", basename(tf), " ===")
  source(tf, local = TRUE)
}

if (!isTRUE(options$skip_smokes)) {
  message("\n=== UI smoke scripts ===")
  integration_run_smokes()
}

summary <- integration_summarize()
integration_write_reports()

cat("\n--- Integration stress summary ---\n")
cat("Pass:", summary$pass, " Fail:", summary$fail, " Skip:", summary$skip, "\n")
if (summary$fail > 0L) quit(status = 1L)
