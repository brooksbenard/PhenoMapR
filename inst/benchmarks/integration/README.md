# PhenoMapR real-data integration stress tests

Permanent harness for manuscript-confidence validation using **real cohorts** (not synthetic fixtures). Complements `tests/testthat/` unit tests and `inst/benchmarks/*_smoke.R` UI source checks.

## Quick start

From the package root:

```bash
Rscript inst/benchmarks/integration/run_integration_stress.R
```

Refresh golden marker fixtures after intentional algorithm changes:

```bash
Rscript inst/benchmarks/integration/run_integration_stress.R --update-fixtures
```

Performance / memory tier (large CRA001160, pseudobulk):

```bash
Rscript inst/benchmarks/integration/run_integration_stress.R --stress --full
```

## Data setup

Tests look for files under `vignettes/` first, then env vars documented in [`vignettes/README.md`](../../vignettes/README.md). If missing, the runner attempts Google Drive download via **googledrive**.

| Key | Default filename |
|-----|------------------|
| CRA001160 | `PAAD_CRA001160_expression.h5`, metadata TSV, optional Seurat RDS |
| HT270 spatial | `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds` |
| GSE205154 bulk | `GSE205154_GPL20301_expression.rds`, `GSE205154_info.rds` |
| GSE253260 bulk | `GSE253260_expression.rds`, `GSE253260_info.rds` |
| CosMx h5ad | `~/Downloads/seurat_test_core_46.h5ad` or `PHENOMAPR_COSMX_H5AD` |
| FSHD validation | `vignettes/FSHD_validation_*.rds` |

Optional extra datasets: copy `local_paths.R.example` to `local_paths.R` and list paths.

## Reports

HTML and JSON reports are written to `inst/benchmarks/integration/reports/`. Golden fixtures live in `fixtures/expected/` (small RDS files, safe to commit).

## CI note

This harness is **not** wired into `R CMD check` by default (runtime and data size). Run locally or on a nightly runner with vignette data cached.
