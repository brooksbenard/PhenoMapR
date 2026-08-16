# Revision checklist (Nature Methods review plan)

| Plan todo | Status | Where |
|-----------|--------|--------|
| Remove Discussion placeholder + Limitations | Done | `DISCUSSION_revised.md` |
| Trim main text; move Fig. 3 target panels to ED | Done | `MAIN_TEXT_TRIMS.md`, `EXTENDED_DATA_target_nomination.md` |
| Simulation + external/held-out validation path | Done | `simulate_phenotype_recovery.R` (ran; AUROC≈1 planted, ≈0.5 null); `summarize_heldout_bulk_for_manuscript.R` + `bulk_paad/` entrypoints |
| Lead with PRECOG→TCGA + held-out metrics | Done | `RESULTS_bulk_reframed.md` |
| Align DGIdb/tail %/shrinkage; tone down actionability | Done | `METHODS_revised.md`; `R/PhenoMap.R` + Rd; `README.md`; `NEWS.md`; existing targeted DGIdb helpers |

**Authoring note:** Apply the markdown drop-ins into the Word/Google Doc source of the PDF before resubmission. Re-run `bulk_paad` held-out experiments when data are available so `reports/bulk_heldout_*` is populated.
