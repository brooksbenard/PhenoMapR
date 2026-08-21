# Vignette Data

The vignettes use data files that are too large for GitHub. Example code downloads them with **googledrive** (`drive_deauth()` — no Google login). If a file is already in the working directory it is reused.

**PhenoMapR data on Google Drive:** [Vignettes folder](https://drive.google.com/drive/folders/1unzVigogMy6XT6KwiFJ2AnKSLsMAzl2X)

| File | Drive file ID | Vignette |
|------|----------------|----------|
| `PAAD_CRA001160_expression.h5` | `1iFWJa13s5UClrP362CQtAAE7KYEo8iBc` | Single-cell — **full** CRA001160 cohort (TISCH2 log-normalized H5) |
| `PAAD_CRA001160_CellMetainfo_table.tsv` | `1yC7Vw3oQ2APB6ZlK7BUl-2DLKGJhFbGN` | Single-cell — cell metadata |
| `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test.hgnc.rds` | `1OkIr7ksAWxKVjtdlGqYHMidvHZZsySEE` | Spatial transcriptomics (spots) |
| `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds` | `1gcOyLriW9bIFNbDuQN6Vi1UsrMGKDxll` | Spatial transcriptomics (CytoSPACE) |
| `GSE205154_GPL20301_expression.rds` | `1RxKOOtsWkTkoooQld8N7xU7lY2Sh4P2b` | Bulk survival |
| `GSE205154_info.rds` | `1isX8lV9UVl_a1YLiopRlWUchPYPUTPXo` | Bulk survival |
| `GSE253260_expression.rds` | `1RaaNkdn-k1fvjWXJuipBmrf_NoqNjX-g` | Custom reference |
| `GSE253260_info.rds` | `15mw6GighBth2iZ_QQS695QRMymfYUy9W` | Custom reference |

**Single-cell vs Shiny:** the single-cell article always uses the **full** CRA001160 H5 + metadata from Drive. The Shiny quick demo uses only the small bundled subset `inst/extdata/shiny/PAAD_CRA001160_demo_5000.rds` (regenerate with `Rscript tools/make_shiny_demo_bundle.R`).

**Pre-rendered spatial figures** live under `inst/figures/` (and `inst/figures/spatial_pair_maps/`). Regenerate with `Rscript scripts/render_spatial_colocalization_heatmap.R` and related scripts.

---

## Setup

1. Install Suggests: `install.packages("PhenoMapR", dependencies = TRUE)` (includes **googledrive**).
2. Copy the load chunk from any article and run it: missing files download into your working directory.
3. If a Drive download fails, download the files in your browser from the folder above and place them in your working directory.
