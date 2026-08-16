# Online Methods — revisions (align with package + review)

Paste/replace the corresponding Methods subsections. Keep total Online Methods near ~3,000 words.

---

## Reference matrices shipped with PhenoMapR (disclose shrinkage)

Built-in PRECOG, TCGA, pediatric, and ICI meta-z tables are distributed in compressed form for package size. During packaging, entries with absolute z-score at or below the default scoring cutoff (|z| ≤ 2 for PRECOG/TCGA; a lower threshold for ICI) are set to missing and genes that are missing in every cancer column are dropped (`data-raw/shrink_reference_data.R`). Consequently, changing `z_score_cutoff` below these packaging thresholds cannot recover discarded associations unless users supply denser upstream tables or a custom reference via `derive_reference_from_bulk()`. This limitation is documented in the package help for `precog` / `tcga`.

---

## Phenotype extremes (tail percentile)

Unless noted, single-cell and spatial phenotype extremes used for marker discovery were the top and bottom **5%** of PhenoMapR scores (`define_phenotype_groups(percentile = 0.05)`), matching the package default and vignettes.  

**DepMap enrichment analyses** used a **10%** phenotype tail (`define_phenotype_groups(percentile = 0.10)`). A 5% tail leaves too few favorable malignant-ductal cells for within–cell-type marker discovery and GSEA under the same log2FC / FDR filters; do not mix 5% and 10% DepMap results without labeling.

---

## Therapeutic / DGIdb annotation (tone + targeted modality)

Adverse markers were annotated with DGIdb v5 flat files for **hypothesis generation only**. Genes were considered DGIdb-tractable when they had ≥1 **targeted modality** interaction (inhibitor, antibody, binder, antagonist, blocker, degrader, or equivalent drug-name cues), not merely:

- a `CLINICALLY ACTIONABLE` panel/biomarker category, or  
- an untyped association to a broad oncology compound.

Tier assignment (package helpers in `manuscript/benchmarks/methodology/dgidb_helpers.R`):

- **A:** targeted modality plus approved / antineoplastic / immunotherapy support on a targeted row  
- **B:** targeted modality without that clinical drug evidence  
- **C:** DGIdb gene-class category without targeted modality  
- **D:** curated fallback list only  
- **E:** none  

Optional Broad Drug Repurposing Hub annotations are reported in a **separate** analysis (`refine_pdac_targets_repurposing_hub.R`) and are not required for core PhenoMapR claims. DGIdb/Hub labels are not functional target validation.

---

## Benchmarking (replace accuracy wording)

Runtime and memory were measured on CRA001160 subsets after method-specific preprocessing, around each tool’s core scoring call (Extended Data Fig. 1). These comparisons speak to computational cost, not biological accuracy.  

Accuracy-oriented evaluation uses (i) **simulations** with planted adverse cells and known z signatures (`simulate_phenotype_recovery.R`; AUROC/AUPRC vs planted labels; null effect control), and (ii) **held-out / external bulk** survival discrimination (PRECOG→TCGA; GSE253260 train/test splits; multi-cohort scripts under `bulk_paad/`). Concordance with other continuous scorers is reported as secondary.

---

## Scoring (notation fix)

Retain genes with **|z| ≥ 2** (default) overlapping the expression matrix, then compute the weighted sum (R `crossprod`) between the z vector and the expression matrix. Positive bulk z indicates association with the adverse / phenotype-positive direction as defined for that reference.
