# PhenoMapR 0.1.0.9000

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
