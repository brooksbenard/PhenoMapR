# Vignette subsets

`PAAD_CRA001160_vignette_subset.rds` is a stratified CRA001160 subset
(~2,500 genes × 6,000 cells) with TISCH-style metadata. The single-cell
pkgdown article prefers this file when `CI` / `GITHUB_ACTIONS` is set,
because the full H5/Seurat objects on Google Drive often hit download
quotas.

Regenerate locally (requires full seurat + metadata under `vignettes/`):

```r
# See conversation history / tools: sample cells, retain high-signal genes,
# saveRDS(..., compress = "xz")
```
