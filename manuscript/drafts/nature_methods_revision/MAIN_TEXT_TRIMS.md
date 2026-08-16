# Main-text trims and figure reorganization

**Goal:** Bring narrative main text to ≤ ~3,000 words (Nature Methods Article) and free display-item budget by moving clinical-target panels out of Fig. 3.

Estimated current narrative (excluding legends): ~3,200–3,300 words. Target cut: **350–600 words**.

---

## Figure plan (keep ≤6 main display items)

| Item | Keep in main? | Action |
|------|---------------|--------|
| Fig. 1 Overview | Yes | Keep; shorten legend |
| Fig. 2 Bulk / PRECOG→TCGA + held-out PAAD | Yes | Reframe (see `RESULTS_bulk_reframed.md`) |
| Fig. 3 PAAD single-cell biology | Yes | **Keep only** UMAP, score distributions, top cell-type contrasts, marker heatmap overview |
| Fig. 3 DGIdb / specificity–trajectory / drug panels | **No** | Move to **Extended Data Fig. X** (“Hypothesis-generating target nomination”) |
| Fig. 4 Spatial | Yes | Clarify spot-native vs CytoSPACE-mapped panels in legend |
| Fig. 5 AD | Yes | Compress GWAS enrichment detail to ED if needed |
| Fig. 6 FSHD | Yes | Keep; shorten legend |

Do **not** add a 7th main figure/table.

---

## Text cuts (Results)

### Development and benchmarking (cut ~120 words)
- Replace “retaining accuracy of results” with a precise claim: concordance with related continuous scorers on one cohort, plus **simulation recovery** (new ED/Supp).
- Demote runtime/memory to one sentence + Extended Data; do not lead the paper on speed alone.
- Fix count: use **eleven** methods if that matches Methods, or update Methods to twelve—keep consistent.
- Delete marketing: “only tool… Shiny App”.

**Replacement paragraph (shorter):**

> To place PhenoMapR among related approaches, we compared runtime and memory to eleven phenotype-mapping tools on matched PAAD inputs after method-specific preprocessing (Extended Data Fig. 1). PhenoMapR was among the fastest and most memory-efficient methods. Continuous scores correlated with several related tools (for example PIPET, SCIPAC, scSurv), but concordance is not a substitute for ground-truth recovery; we therefore also evaluate recovery of planted phenotype effects in simulations (Extended Data Fig. 1) and emphasize held-out bulk and independent cohort validations below.

### Bulk section (cut ~80 words; reframe)
Use `RESULTS_bulk_reframed.md` (leads with PRECOG→TCGA and held-out metrics; TCGA→TCGA only as ED).

### PAAD single-cell (cut ~100–150 words)
- Remove or move “clinically actionable” / DGIdb narrative to ED.
- Keep malignant ductal enrichment + marker examples that support biology.
- Tone: “hypothesis-generating candidate genes” not “clinically actionable targets.”

**Tone replacement:**

> Within adverse-enriched lineages we report genes with strong Fav→Adv expression trajectories and, where available, indexed drug–gene modalities in DGIdb as **hypothesis-generating** candidates (Extended Data Fig. X). These annotations do not constitute functional target validation.

### Spatial / AD / FSHD (cut ~50–80 words each if needed)
- Prefer one clear validation claim per subsection; move secondary plots to ED.
- FSHD: state CyteTypeR-derived labels explicitly once; avoid over-claiming “therapeutic targets.”

### Discussion
Use `DISCUSSION_revised.md` (already shorter; includes Limitations; no placeholder).

---

## Copy-edits to apply globally

| Find | Replace |
|------|---------|
| ADD NEW RESULTS HERE | *(delete)* |
| PhenoMaR | PhenoMapR |
| clinically-actional | clinically actionable *(or better: remove “clinically”)* |
| rak-sum | rank-sum |
| log-rank text | log-rank test |
| FSHA | FSHD |
| scPANDA | ctPANDA *(when referring to Tang et al.)* |
| \|Z ≥ 2\| | \|Z\| ≥ 2 |
| *** q≤0.001, *** q≤0.001 | *** q≤0.001 *(once)* |
| respevely | respectively |

---

## References
Target ≤50 highly cited method/data papers in the main list; move tertiary dataset citations to Methods or SI if editors request.
