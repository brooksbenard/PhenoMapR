# Results — bulk validation (reframed)

**Action:** Replace the bulk Results subsection so non-circular designs lead, and same-cohort / TCGA→TCGA analyses are secondary (Extended Data).

---

### Pan-cancer stratification with independent signatures

We first asked whether PhenoMapR scores derived from **PRECOG** cancer meta-z signatures stratify overall survival in matched **TCGA** expression cohorts (Fig. 2). Because PRECOG meta-z values are estimated outside TCGA, this transfer is not a circular reuse of the same Cox fits. Median-split PhenoMapR scores yielded significant hazard ratios in the majority of PRECOG–TCGA matched cancers; non-significant cases were enriched for weak PRECOG–TCGA signature concordance and/or small TCGA sample size (Fig. 2b).  

As a sensitivity analysis, scoring TCGA with **built-in TCGA** meta-z signatures also stratified many cancers (Extended Data Fig. 2), but those results can partly reflect information shared with the reference construction and are not used as the primary claim.

### Custom PAAD signatures with held-out evaluation

Using clinically annotated bulk PAAD expression (GSE253260), we derived custom overall-survival signatures with `derive_reference_from_bulk()`. **Primary metrics** are from repeated train/test splits within disease stage: signatures fit on training patients only were applied to held-out patients, and discrimination was summarized by Cox models on the continuous score and by median-split log-rank tests (Methods; Extended Data Fig. 2). Held-out analyses supported stage-resolved risk separation for non-metastatic stages, with weaker performance in metastatic disease where sample size and outcome homogeneity limit discrimination.  

Full-cohort median splits are shown only as descriptive visualizations (Fig. 2c–e) and are not the basis for performance claims. Stage-agnostic scores were not strongly associated with pathologic stage, consistent with residual risk information orthogonal to stage labels in this cohort.
