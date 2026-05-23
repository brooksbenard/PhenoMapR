# PhenoMapR Shiny app

An interactive front end for the `PhenoMapR` R package. The app walks users
through the full workflow — data upload, reference choice, scoring,
phenotype-group tagging, marker discovery, and diagnostics — without writing
any R code.

The app talks to the same exported functions as the package
(`PhenoMap()`, `define_phenotype_groups()`, `find_phenotype_markers()`,
`derive_reference_from_bulk()`, `get_gene_coverage()`,
`get_top_prognostic_genes()`, `plot_score_distribution()`,
`plot_reference_signature()`), so anything you can do here you can also do
from the R console / a script with identical results.

---

## Quick start

### From an installed `PhenoMapR` package

```r
# install once
remotes::install_github("brooksbenard/PhenoMapR")

# launch
PhenoMapR::run_app()
```

A browser tab opens automatically at `http://127.0.0.1:<random port>`.

### From a clone (development / customization)

```r
# from the repo root
shiny::runApp("inst/shiny")
```

### On a remote server (workstation / HPC node)

```r
PhenoMapR::run_app(host = "0.0.0.0", port = 3838, launch.browser = FALSE)
```

Then open `http://<server>:3838` from any machine that can reach the server.
The app processes everything locally on the server — no data is uploaded
to a third party.

---

## Supported input files

| Type | Extension | Notes |
| ---- | --------- | ----- |
| R object | `.rds` | Most flexible: Seurat (v4/v5), SingleCellExperiment, SpatialExperiment, matrix, sparse `Matrix`, or AnnData (saved with `reticulate`). |
| Tabular matrix | `.tsv`, `.csv`, `.txt` | First column = gene IDs (HUGO symbols), remaining columns = samples/cells. |
| 10X HDF5 | `.h5` | Loaded via `Seurat::Read10X_h5()` or `hdf5r` fallback. |
| AnnData | `.h5ad` | Requires `reticulate` + Python `anndata`. |

You can also supply an **optional cell-metadata table** (tabular or `.rds`)
to enable cell-type-specific marker analyses and source-aware visualizations.

---

## Workflow overview (mirrors the package vignettes)

1. **Data** — upload expression, optional metadata. The app reports the data
   class, dimensions, and any input warnings (Ensembl IDs detected, axes
   flipped, etc.).
2. **Reference** — choose a built-in prognostic reference (PRECOG / TCGA /
   Pediatric PRECOG / ICI PRECOG) **or** supply a custom one by uploading
   a gene × z-score table **or** deriving one with
   `derive_reference_from_bulk()` on a separate bulk + phenotype upload
   (binary / continuous / survival outcomes). Gene-coverage stats and a
   reference-signature plot are computed automatically.
3. **Score** — run `PhenoMap()` with your choice of assay / layer, optional
   pseudobulk aggregation, and a `reference_sign` direction. The result is
   a downloadable TSV with one score column per reference.
4. **Phenotype groups** — `define_phenotype_groups()` partitions cells into
   *Most Adverse* / *Most Favorable* / *Other* by the top and bottom *N* %.
   Per-cell-type group enrichment is shown when a cell-type column is set.
5. **Markers** — `find_phenotype_markers()` with cohort-wide or
   cell-type-specific scope, exposing `min.pct`, `logfc.threshold`,
   `pval_threshold`, and `max_cells_per_ident`. Tables for the adverse and
   favorable tails are downloadable as an `.rds`.
6. **Diagnostics** — list all available cancer types per reference, inspect
   the top prognostic genes, and see how many of *your* genes overlap each
   built-in reference.

---

## Upload size limits

Single-cell `.rds` and `.h5ad` files routinely exceed several GB. The app
removes Shiny's default 5 MB upload cap entirely so you can drop in
multi-GB AnnData / Seurat objects on a workstation or HPC node — the only
real limit is your machine's RAM.

If you want to **cap** uploads (e.g. on a shared web server), pass an
explicit limit when launching:

```r
PhenoMapR::run_app(max_upload_size_mb = 500)        # cap at 500 MB
# or set the option directly before launching:
options(shiny.maxRequestSize = 5 * 1024^3)          # 5 GiB cap
PhenoMapR::run_app()
```

If you're running on a server with limited memory, prefer uploading a
preprocessed `.rds` (subset to common genes, normalized) rather than a raw
10X HDF5.

---

## Optional dependencies

* `ComplexHeatmap` + `circlize` — reference signature heatmap.
* `Seurat` — `.h5` 10X parsing and the full Seurat code path.
* `SingleCellExperiment` / `SpatialExperiment` / `SummarizedExperiment` — SCE inputs.
* `reticulate` + Python `anndata` — `.h5ad` parsing.
* `data.table` — faster tabular parsing.
* `presto` — much faster Wilcoxon for marker discovery.
* `hdf5r` — generic HDF5 fallback when Seurat is unavailable.
* `ggsignif` — Wilcoxon brackets on the score-by-cell-type plot.
* [`ggchicklet2`](https://github.com/brooksbenard/ggchicklet2) — drop-in
  rounded variants of `geom_bar()` / `geom_histogram()` / `geom_boxplot()`
  used for the bar / box / histogram visualisations in the app. When
  installed (`remotes::install_github("brooksbenard/ggchicklet2")`) the
  app automatically picks up the rounded geoms; otherwise it falls back
  to the standard ggplot2 geoms.

All of these are optional; the app degrades gracefully when they are missing.

---

## Reporting issues

Open an issue at <https://github.com/brooksbenard/PhenoMapR/issues>.
Please include the output of:

```r
sessionInfo()
packageVersion("PhenoMapR")
packageVersion("shiny")
```
