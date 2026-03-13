# Vignette Data

The vignettes use data files that are too large for GitHub. You can place them locally, let the vignettes download from Google Drive via **googledrive**, or set direct-download URLs.

**PhenoMapR data on Google Drive:** [Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) — subfolders **Bulk_Expression**, **Single_Cell**, and **Spatial_Transcriptomics** contain the same files used by the vignettes. The vignettes use these fixed file IDs when downloading with **googledrive**:

| File | Drive file ID | Vignette |
|------|----------------|----------|
| `PAAD_CRA001160_expression.h5` | `1PolTXggREz8XmhutCLTQJGCfKxFAzqMl` | Single-cell (CRA001160) |
| `PAAD_CRA001160_CellMetainfo_table.tsv` | `17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu` | Single-cell (CRA001160) |
| `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds` | `1HM0dBrQnaNsdm5mnq23aaQ2ILofJ0_vj` | Spatial transcriptomics |
| `GSE205154.GPL20301.matrix.txt` | `1Vk4KCQWF9ikpAuMsjFzVDCoy1TzDl2rN` | Bulk survival & Custom reference |
| `GSE205154.info.txt` | `1omAA2kfVn-nyyZfcc4vBhRFogC6cuoNQ` | Bulk survival & Custom reference |

---

## Option 1: Place files locally (simplest)

Download from the [Drive Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) (open each subfolder and download the files), then put them in `vignettes/`, `Vignettes/`, or the package root. If a file is present, the vignette uses it and does not download.

---

## Option 2: Automatic download from Google Drive (googledrive)

If a data file is **not** found locally, each vignette will try to download it from Google Drive using the **googledrive** package (no login: `drive_deauth()` is used for public links).

1. Install the package with Suggests: `install.packages("PhenoMapR", dependencies = TRUE)` or `install.packages("googledrive")`.
2. Run or knit the vignette from the package root (so paths like `vignettes/PAAD_CRA001160_expression.h5` resolve). The first run will download the file(s) into `vignettes/` (or the current directory).

Files must be shared **Anyone with the link**. The vignettes use the file IDs in the table above.

---

## Option 3: Load from URL (environment variables)

You can supply **direct download** URLs via environment variables. This is useful in CI (e.g. GitHub Actions) where you store the URLs in secrets; the vignettes fall back to these if the file is still missing after trying local paths and googledrive.

### Step 1: Get each file’s direct download link from Drive

1. Open the [Vignettes folder](https://drive.google.com/drive/folders/1rKGZBX7sa_Iq8AJb1wcxiRc3oD6v6B5n) on Google Drive.
2. Open **Single_Cell**, then for each of `PAAD_CRA001160_expression.h5` and `PAAD_CRA001160_CellMetainfo_table.tsv`:
   - Right-click the file → **Share** → set to **Anyone with the link** (Viewer).
   - Copy the link. It looks like: `https://drive.google.com/file/d/XXXXXXXXXX/view?usp=sharing`
   - The **file ID** is the part between `/d/` and `/view` (e.g. `XXXXXXXXXX`).
   - The **direct download URL** is: `https://drive.google.com/uc?export=download&id=XXXXXXXXXX`
3. Do the same for the file in **Spatial_Transcriptomics** (`HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds`).
4. (Optional) For **Bulk_Expression**, get direct download URLs for `GSE205154.GPL20301.matrix.txt` and `GSE205154.info.txt` the same way.

### Step 2: Set the URLs and run the vignettes

In R, **before** knitting or sourcing the vignettes, run (you can use the file IDs from the table at the top):

```r
# Use the file IDs from the table at the top of this README.
Sys.setenv(
  PHENOMAPR_CRA001160_H5_URL    = "https://drive.google.com/uc?export=download&id=1PolTXggREz8XmhutCLTQJGCfKxFAzqMl",
  PHENOMAPR_CRA001160_META_URL  = "https://drive.google.com/uc?export=download&id=17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu",
  PHENOMAPR_SPATIAL_RDS_URL     = "https://drive.google.com/uc?export=download&id=1HM0dBrQnaNsdm5mnq23aaQ2ILofJ0_vj",
  PHENOMAPR_GSE205154_MATRIX_URL = "https://drive.google.com/uc?export=download&id=1Vk4KCQWF9ikpAuMsjFzVDCoy1TzDl2rN",
  PHENOMAPR_GSE205154_INFO_URL   = "https://drive.google.com/uc?export=download&id=1omAA2kfVn-nyyZfcc4vBhRFogC6cuoNQ"
)

# Then knit or run the vignettes, e.g.:
# rmarkdown::render("vignettes/single-cell.Rmd")
# rmarkdown::render("vignettes/spatial-transcriptomics.Rmd")
# rmarkdown::render("vignettes/bulk-survival.Rmd")
# rmarkdown::render("vignettes/custom-reference.Rmd")
```

Or set them in the shell before starting R:

```bash
export PHENOMAPR_CRA001160_H5_URL="https://drive.google.com/uc?export=download&id=1PolTXggREz8XmhutCLTQJGCfKxFAzqMl"
export PHENOMAPR_CRA001160_META_URL="https://drive.google.com/uc?export=download&id=17mqxnKOZJn0jW2iD9RV0wZeWsilAIwdu"
export PHENOMAPR_SPATIAL_RDS_URL="https://drive.google.com/uc?export=download&id=1HM0dBrQnaNsdm5mnq23aaQ2ILofJ0_vj"
export PHENOMAPR_GSE205154_MATRIX_URL="https://drive.google.com/uc?export=download&id=1Vk4KCQWF9ikpAuMsjFzVDCoy1TzDl2rN"
export PHENOMAPR_GSE205154_INFO_URL="https://drive.google.com/uc?export=download&id=1omAA2kfVn-nyyZfcc4vBhRFogC6cuoNQ"
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
  PHENOMAPR_CRA001160_H5_URL: ${{ secrets.PHENOMAPR_CRA001160_H5_URL }}
  PHENOMAPR_CRA001160_META_URL: ${{ secrets.PHENOMAPR_CRA001160_META_URL }}
  PHENOMAPR_SPATIAL_RDS_URL: ${{ secrets.PHENOMAPR_SPATIAL_RDS_URL }}
  PHENOMAPR_GSE205154_MATRIX_URL: ${{ secrets.PHENOMAPR_GSE205154_MATRIX_URL }}
  PHENOMAPR_GSE205154_INFO_URL: ${{ secrets.PHENOMAPR_GSE205154_INFO_URL }}
```
