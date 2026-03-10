# Vignette Data

The vignettes use data files that are too large for GitHub. You can place them locally or load them from URLs (e.g. Google Drive).

**PhenoMapR data on Google Drive:** [Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) — subfolders **Bulk_Expression**, **Single_Cell**, and **Spatial_Transcriptomics** contain the same files used by the vignettes.

---

## Option 1: Place files locally (simplest)

Download from the [Drive Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) (open each subfolder and download the files), then put them in `vignettes/`, `Vignettes/`, or the package root:

| File | Drive subfolder | Vignette |
|------|-----------------|----------|
| `PAAD_GSE111672_seurat.rds` | Single_Cell | Single-cell (GSE111672) |
| `PAAD_CRA001160_seurat.rds` | Single_Cell | Single-cell (CRA001160) |
| `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds` | Spatial_Transcriptomics | Spatial transcriptomics |
| `GSE205154.GPL20301.matrix.txt` | Bulk_Expression | Bulk survival & Custom reference |
| `GSE205154.info.txt` | Bulk_Expression | Bulk survival & Custom reference |

---

## Option 2: Load from URL (Google Drive direct links)

The single-cell and spatial vignettes can **download** files from URLs if you set environment variables. Use **direct download** URLs (see below for how to get them from Google Drive).

### Step 1: Get each file’s direct download link from Drive

1. Open the [Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) on Google Drive.
2. Open **Single_Cell**, then for each of `PAAD_GSE111672_seurat.rds` and `PAAD_CRA001160_seurat.rds`:
   - Right-click the file → **Share** → set to **Anyone with the link** (Viewer).
   - Copy the link. It looks like: `https://drive.google.com/file/d/XXXXXXXXXX/view?usp=sharing`
   - The **file ID** is the part between `/d/` and `/view` (e.g. `XXXXXXXXXX`).
   - The **direct download URL** is: `https://drive.google.com/uc?export=download&id=XXXXXXXXXX`
3. Do the same for the file in **Spatial_Transcriptomics** (`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`).
4. (Optional) For **Bulk_Expression**, get direct download URLs for `GSE205154.GPL20301.matrix.txt` and `GSE205154.info.txt` the same way.

### Step 2: Set the URLs in R and run the vignettes

In R, **before** knitting or sourcing the vignettes, run (replace the `id=...` parts with your actual file IDs):

```r
# Replace YOUR_GSE111672_ID, YOUR_CRA001160_ID, YOUR_SPATIAL_ID (and optionally bulk IDs) with the IDs from Step 1.
Sys.setenv(
  PHENOMAPR_GSE111672_RDS_URL   = "https://drive.google.com/uc?export=download&id=YOUR_GSE111672_ID",
  PHENOMAPR_CRA001160_RDS_URL   = "https://drive.google.com/uc?export=download&id=YOUR_CRA001160_ID",
  PHENOMAPR_SPATIAL_RDS_URL     = "https://drive.google.com/uc?export=download&id=YOUR_SPATIAL_ID",
  PHENOMAPR_GSE205154_MATRIX_URL = "https://drive.google.com/uc?export=download&id=YOUR_MATRIX_ID",
  PHENOMAPR_GSE205154_INFO_URL   = "https://drive.google.com/uc?export=download&id=YOUR_INFO_ID"
)

# Then knit or run the vignettes, e.g.:
# rmarkdown::render("vignettes/gse111672-single-cell.Rmd")
# rmarkdown::render("vignettes/spatial-transcriptomics.Rmd")
# rmarkdown::render("vignettes/gse205154-bulk-survival.Rmd")
# rmarkdown::render("vignettes/gse205154-custom-reference.Rmd")
```

Or set them in the shell before starting R:

```bash
export PHENOMAPR_GSE111672_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_GSE111672_ID"
export PHENOMAPR_CRA001160_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_CRA001160_ID"
export PHENOMAPR_SPATIAL_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_SPATIAL_ID"
export PHENOMAPR_GSE205154_MATRIX_URL="https://drive.google.com/uc?export=download&id=YOUR_MATRIX_ID"
export PHENOMAPR_GSE205154_INFO_URL="https://drive.google.com/uc?export=download&id=YOUR_INFO_ID"
R
```

### Large files on Google Drive

For large files, Google may show a warning page instead of downloading. If the automatic download fails:

- **Easiest:** Use Option 1 — download the file in your browser from Drive and place it in `vignettes/`.
- Or use the [gdown](https://github.com/wkentaro/gdown) Python tool with the file ID: `gdown FILE_ID`.

### GitHub Actions (CI)

Add the direct download URLs as repository secrets (e.g. `PHENOMAPR_GSE111672_RDS_URL`, etc.), then in your workflow set:

```yaml
env:
  PHENOMAPR_GSE111672_RDS_URL: ${{ secrets.PHENOMAPR_GSE111672_RDS_URL }}
  PHENOMAPR_CRA001160_RDS_URL: ${{ secrets.PHENOMAPR_CRA001160_RDS_URL }}
  PHENOMAPR_SPATIAL_RDS_URL: ${{ secrets.PHENOMAPR_SPATIAL_RDS_URL }}
  PHENOMAPR_GSE205154_MATRIX_URL: ${{ secrets.PHENOMAPR_GSE205154_MATRIX_URL }}
  PHENOMAPR_GSE205154_INFO_URL: ${{ secrets.PHENOMAPR_GSE205154_INFO_URL }}
```
