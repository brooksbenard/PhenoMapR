# PhenoMapR manuscript analyses

This directory holds analyses, figures, and draft R Markdown used for the PhenoMapR manuscript. It is **not** part of the installable R package at the repository root.

Most large or generated content under `manuscript/` is **gitignored** so GitHub stays focused on the R package, Shiny app, vignettes, and pkgdown site. Keep local manuscript data and methodology outputs here; do not commit them.

## Layout

| Path | Contents |
|------|----------|
| `benchmarks/` | Method benchmarking, integration stress tests, FSHD/PAAD validation pipelines |
| `drafts/` | Draft analyses (Alzheimer's, FSHD, etc.) |
| `figures/` | Manuscript-only figures (reference coverage, forest plots, candidates) |
| `data/` | FSHD validation caches and other large analysis datasets (gitignored) |
| `scripts/` | Manuscript analysis scripts (TCGA/PRECOG forests, path helpers) |
| `results/` | Local analysis outputs (gitignored) |
| `docs/` | Local literature notes (gitignored) |

Vignette data, pre-rendered spatial figures, and vignette render scripts live at the **package root** (`vignettes/`, `inst/figures/`, `inst/data/`, `scripts/`).

## Prerequisites

```r
pkgload::load_all(".")   # from repo root
# or:
remotes::install_github("brooksbenard/PhenoMapR")
```

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `PHENOMAPR_REPO` | repo root | Package root (`DESCRIPTION`) |
| `PHENOMAPR_MANUSCRIPT` | `manuscript/` | Manuscript root |
| `PHENOMAPR_VIGNETTE_DATA` | `vignettes/` | Large vignette datasets |

Path helpers live in `manuscript/scripts/paths.R` (`phenomapr_pkg_root()`, `manuscript_root()`, `manuscript_benchmarks_dir()`, `manuscript_results_dir()`, `tcga_data_dir()`). Prefer these over hard-coded absolute paths; do not write analysis outputs outside this repository (legacy `~/Desktop/scIMPEL` layouts are retired).

## Common tasks

```bash
# From repo root
PHENOMAPR_REPO=$PWD Rscript manuscript/benchmarks/integration/run_integration_stress.R

# Manuscript figures / TCGA forests (defaults write to manuscript/results/)
Rscript manuscript/scripts/run_tcga_precog_forest_incremental.R

# Methodology sensitivity experiments
Rscript manuscript/benchmarks/methodology/run_methodology_experiments.R --demo

# Regenerate spatial vignette assets (from repo root)
Rscript scripts/render_spatial_colocalization_heatmap.R
```
