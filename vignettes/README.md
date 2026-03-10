# Vignette Data

The single-cell and spatial transcriptomics vignettes require preprocessed RDS files that are too large for GitHub. You can either place them locally or load them from a URL (e.g. Google Drive).

## Option 1: Place files locally

Place the RDS files in `vignettes/`, `Vignettes/`, or the package root:

| File | Vignette |
|------|----------|
| `PAAD_GSE111672_seurat.rds` | Single-cell (GSE111672) |
| `PAAD_CRA001160_seurat.rds` | Single-cell (CRA001160) |
| `HT270P1-S1H2Fc2U1Z1Bs1-H2Bs2-Test_processed.rds` | Spatial transcriptomics |

## Option 2: Load from URL (Google Drive, Zenodo, etc.)

Set environment variables with **direct download URLs** before knitting or running the vignettes:

```r
# In R, before running vignettes:
Sys.setenv(
  PHENOMAPR_GSE111672_RDS_URL = "https://drive.google.com/uc?export=download&id=YOUR_FILE_ID",
  PHENOMAPR_CRA001160_RDS_URL = "https://drive.google.com/uc?export=download&id=YOUR_FILE_ID",
  PHENOMAPR_SPATIAL_RDS_URL  = "https://drive.google.com/uc?export=download&id=YOUR_FILE_ID"
)
```

Or in your shell (e.g. before `R CMD build` or CI):

```bash
export PHENOMAPR_GSE111672_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_FILE_ID"
export PHENOMAPR_CRA001160_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_FILE_ID"
export PHENOMAPR_SPATIAL_RDS_URL="https://drive.google.com/uc?export=download&id=YOUR_FILE_ID"
```

### Google Drive direct links

1. Upload the RDS file to Google Drive.
2. Right-click → Share → set to "Anyone with the link".
3. Copy the link. The file ID is the long string between `/d/` and `/view`:
   - `https://drive.google.com/file/d/FILE_ID_HERE/view?usp=sharing`
4. Use: `https://drive.google.com/uc?export=download&id=FILE_ID_HERE`

### GitHub Actions (CI)

Add these as repository secrets, then reference them in your workflow:

```yaml
env:
  PHENOMAPR_GSE111672_RDS_URL: ${{ secrets.PHENOMAPR_GSE111672_RDS_URL }}
  PHENOMAPR_CRA001160_RDS_URL: ${{ secrets.PHENOMAPR_CRA001160_RDS_URL }}
  PHENOMAPR_SPATIAL_RDS_URL: ${{ secrets.PHENOMAPR_SPATIAL_RDS_URL }}
```

### Zenodo / other hosts

Use the direct download URL (not a landing page). For Zenodo, use the "Download" link for the file.
