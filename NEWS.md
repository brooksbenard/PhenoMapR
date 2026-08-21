# PhenoMapR 0.1.0.9000

## Methodology / orthogonal validation

* ctPANDA (Tang et al., *Cancer Cell* 2026): cell-type-matched comparison of
  CRA001160 PhenoMapR adverse/favorable markers vs published crPRG lists lives
  under `manuscript/benchmarks/methodology/` (`--only=ctpanda`), not the
  single-cell vignette. Bundled lists remain at
  `inst/extdata/ctpanda/ctpanda_prgs.rds` (`tools/make_ctpanda_prgs.R`).

## Bug fixes

* Single-cell vignette: load the **full** CRA001160 cohort from Google Drive
  only (H5 + metadata); removed URL/secret/subset fallbacks. The Shiny quick
  demo continues to use only the small bundled
  `PAAD_CRA001160_demo_5000.rds`.
* Shiny quick demo: the bundled `PAAD_CRA001160_demo_5000.rds` is shipped
  again (a `.Rbuildignore` pattern was excluding all of `inst/extdata`), so
  the demo no longer downloads the full CRA001160 cohort at runtime.
* Shiny scoring: the app now passes `layer=`/`slot=` according to the
  installed `PhenoMap()` signature, fixing
  `PhenoMap failed: unused argument (layer = layer)` when the app is newer
  than the installed core.

## Compatibility

* Seurat **v4 and v5**: `PhenoMap()` and `find_phenotype_markers()` use
  Seurat v5 `layer=` by default (`data` / `counts` / `scale.data`) and still
  accept v4 `slot=` as an alias. Internally `GetAssayData` /
  `SetAssayData` use the argument your installed SeuratObject expects, with
  a one-shot fallback to the other name. Spatial Seurat objects honour
  the caller's `layer`/`slot` instead of forcing `counts`.

## Breaking changes

* `PhenoMap()` score columns are now named `PhenoMapR_<reference_label>` instead of
  `weighted_sum_score_<reference_label>`. Derived columns follow the same prefix
  (for example `phenotype_group_PhenoMapR_precog_BRCA`, `empirical_p_PhenoMapR_precog_BRCA`).
  Re-run scoring on saved objects or rename metadata columns manually.
