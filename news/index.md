# Changelog

## PhenoMapR 0.1.0.9000

### Breaking changes

- [`PhenoMap()`](https://brooksbenard.github.io/PhenoMapR/reference/PhenoMap.md)
  score columns are now named `PhenoMapR_<reference_label>` instead of
  `weighted_sum_score_<reference_label>`. Derived columns follow the
  same prefix (for example `phenotype_group_PhenoMapR_precog_BRCA`,
  `empirical_p_PhenoMapR_precog_BRCA`). Re-run scoring on saved objects
  or rename metadata columns manually.
