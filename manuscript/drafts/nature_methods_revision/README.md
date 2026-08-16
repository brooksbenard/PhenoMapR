# Nature Methods revision package

This folder holds **drop-in manuscript revisions** responding to an internal Nature Methods–style review of `PhenoMapR Manuscript (1).pdf`. The publisher PDF has no editable Word/TeX source in-repo; apply these texts in the authoring document before resubmission.

| File | Purpose |
|------|---------|
| [DISCUSSION_revised.md](DISCUSSION_revised.md) | Removes `ADD NEW RESULTS HERE`; adds Limitations |
| [MAIN_TEXT_TRIMS.md](MAIN_TEXT_TRIMS.md) | Length/figure budget + copy-edits |
| [RESULTS_bulk_reframed.md](RESULTS_bulk_reframed.md) | PRECOG→TCGA + held-out primary claims |
| [METHODS_revised.md](METHODS_revised.md) | Ref shrinkage, tail %, targeted DGIdb, benchmark scope |

## New / supporting analyses

```bash
# Simulation recovery (AUROC/AUPRC vs planted adverse cells)
Rscript manuscript/benchmarks/methodology/simulate_phenotype_recovery.R

# Manifest + guidance for held-out bulk tables
Rscript manuscript/benchmarks/methodology/summarize_heldout_bulk_for_manuscript.R
```

Outputs: `manuscript/benchmarks/methodology/reports/simulation_phenotype_recovery*.csv`, `bulk_heldout_*`.
