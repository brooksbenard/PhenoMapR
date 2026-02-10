# PrognosticScorer Package - Complete Implementation Guide

## Package Overview

PrognosticScorer is an R package that calculates weighted sum prognostic scores for:
- Bulk RNA-seq data
- Single-cell RNA-seq (Seurat, SingleCellExperiment, AnnData)
- Spatial transcriptomics

Using reference z-scores from:
- PRECOG (pan-cancer meta-analysis)
- TCGA (survival analysis)
- Pediatric cancers
- ICI (immunotherapy response)

## Files Created

### Core R Files (go in `R/` directory)
1. **score_expression.R** - Main user-facing function
2. **data_handlers.R** - Input validation and format conversion
3. **weighted_sum_scoring.R** - Core scoring calculations
4. **prognostic_data.R** - Reference data documentation
5. **utils.R** - Helper utilities

### Package Configuration
6. **DESCRIPTION** - Package metadata and dependencies
7. **README.md** - User documentation
8. **introduction.Rmd** - Vignette (goes in `vignettes/`)

### Setup Scripts
9. **setup_data.R** - Prepare and save reference datasets
10. **build_package.R** - Automated build script

## Step-by-Step Setup

### 1. Create Package Structure

```bash
# Option A: Use build_package.R script
Rscript build_package.R

# Option B: Manual setup
```

```r
library(devtools)
library(usethis)

# Create package
create_package("~/PrognosticScorer")
setwd("~/PrognosticScorer")

# Setup infrastructure  
use_mit_license("Your Name")
use_readme_md()
use_roxygen_md()
use_testthat()
use_vignette("introduction")
```

### 2. Copy R Source Files

Copy all 5 R files into the `R/` directory:

```bash
cp score_expression.R ~/PrognosticScorer/R/
cp data_handlers.R ~/PrognosticScorer/R/
cp weighted_sum_scoring.R ~/PrognosticScorer/R/
cp prognostic_data.R ~/PrognosticScorer/R/
cp utils.R ~/PrognosticScorer/R/
```

### 3. Prepare Reference Data

Your reference data should be data.frames with:
- Genes as rownames
- Cancer types as column names  
- Z-scores as values

Example structure:
```r
# precog example
                 BRCA    LUAD    COAD
TP53             3.25   -2.10    1.45
MYC             -1.80    2.35   -0.95
EGFR             0.50    4.20    0.30
```

### 4. Save Data to Package

```r
# Load your reference dataframes into R environment
load("path/to/precog.RData")  # or however you load them
load("path/to/tcga.RData")
load("path/to/pediatric.RData")
load("path/to/ici.RData")

# Create datasets_info mapping
datasets_info <- data.frame(
  dataset_name = c("my_breast_study", "my_lung_study"),
  precog_label = c("BRCA", "LUAD"),
  tcga_label = c("BRCA", "LUAD"),
  pediatric_precog_label = c(NA, NA),
  ici_precog_label = c(NA, "NSCLC"),
  stringsAsFactors = FALSE
)

# Save to package
setwd("~/PrognosticScorer")
usethis::use_data(precog, overwrite = TRUE, compress = "xz")
usethis::use_data(tcga, overwrite = TRUE, compress = "xz")
usethis::use_data(pediatric, overwrite = TRUE, compress = "xz")
usethis::use_data(ici, overwrite = TRUE, compress = "xz")
usethis::use_data(datasets_info, overwrite = TRUE, compress = "xz")
```

### 5. Generate Documentation

```r
library(devtools)
setwd("~/PrognosticScorer")

# Generate documentation from roxygen comments
document()
```

### 6. Check Package

```r
# Run comprehensive checks
check()

# Common issues and fixes:
# - Missing dependencies: Add to DESCRIPTION Imports/Suggests
# - Undocumented functions: Add #' roxygen comments
# - Data documentation: Already included in prognostic_data.R
```

### 7. Build and Install

```r
# Install locally
install()

# Or build source package
build()

# Build binary package
build(binary = TRUE)
```

### 8. Test It Works

```r
library(PrognosticScorer)

# List cancer types
list_cancer_types("precog")

# Create test data
test_matrix <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
rownames(test_matrix) <- paste0("Gene", 1:1000)
colnames(test_matrix) <- paste0("Sample", 1:10)

# Score
scores <- score_expression(
  expression = test_matrix,
  reference = "precog",
  cancer_type = "BRCA"
)

head(scores)
```

## Usage Examples

### Example 1: Bulk RNA-seq

```r
library(PrognosticScorer)

# Load expression data (genes x samples)
bulk_expr <- read.csv("expression_matrix.csv", row.names = 1)

# Calculate scores
scores <- score_expression(
  expression = bulk_expr,
  reference = "tcga",
  cancer_type = "LUAD",
  z_score_cutoff = 2
)

# Normalize to z-scores
scores$zscore <- normalize_scores(scores[, 1])
```

### Example 2: Single-Cell Seurat

```r
library(Seurat)
library(PrognosticScorer)

# Load Seurat object
seurat_obj <- readRDS("seurat_object.rds")

# Score individual cells
scores <- score_expression(
  expression = seurat_obj,
  reference = "precog",
  cancer_type = "BRCA",
  assay = "RNA",
  slot = "data"
)

# Add to metadata
seurat_obj <- add_scores_to_seurat(seurat_obj, scores)

# Visualize
FeaturePlot(seurat_obj, features = "weighted_sum_score_precog_BRCA")
VlnPlot(seurat_obj, features = "weighted_sum_score_precog_BRCA", 
        group.by = "celltype")
```

### Example 3: Pseudobulk Analysis

```r
# Aggregate by patient before scoring
patient_scores <- score_expression(
  expression = seurat_obj,
  reference = "ici_precog",
  cancer_type = "MELANOMA_Metastatic",
  pseudobulk = TRUE,
  group_by = "patient_id"
)

# Each row is now a patient-level score
head(patient_scores)
```

### Example 4: Spatial Transcriptomics

```r
library(Seurat)
library(PrognosticScorer)

# Load spatial object
spatial_obj <- Load10X_Spatial("path/to/spaceranger/outs")

# Score spatial spots
scores <- score_expression(
  expression = spatial_obj,
  reference = "tcga",
  cancer_type = "BRCA",
  assay = "Spatial",
  slot = "counts"
)

# Add and visualize
spatial_obj <- add_scores_to_seurat(spatial_obj, scores)
SpatialFeaturePlot(spatial_obj, 
                   features = "weighted_sum_score_tcga_BRCA")
```

### Example 5: Custom Signatures

```r
# Create custom gene signature
my_signature <- data.frame(
  row.names = c("TP53", "MYC", "EGFR", "PTEN", "KRAS"),
  hypoxia_score = c(3.5, -2.1, 2.8, -1.9, 1.4)
)

# Score with custom signature
scores <- score_expression(
  expression = my_expression_data,
  reference = my_signature,
  z_score_cutoff = 1.5
)
```

## Key Features

### Flexible Input Formats
- Base R matrix/data.frame
- Seurat objects (single-cell and spatial)
- SingleCellExperiment
- SpatialExperiment  
- AnnData (Python, via reticulate)

### Multiple Reference Datasets
- **PRECOG**: Pan-cancer meta-analysis
- **TCGA**: TCGA survival analysis
- **Pediatric**: Pediatric cancers
- **ICI**: Immunotherapy response

### Analysis Options
- Cell-level or pseudobulk scoring
- Adjustable z-score cutoffs
- Custom signature support
- Multiple cancer type scoring

### Utility Functions
- `list_cancer_types()` - Browse available cancer types
- `get_gene_coverage()` - Check gene overlap
- `get_top_prognostic_genes()` - Extract top genes
- `normalize_scores()` - Convert to z-scores
- `plot_score_distribution()` - Visualize scores

## Package Dependencies

### Required (Imports)
- methods
- Matrix
- dplyr
- glue
- progress

### Optional (Suggests)
- Seurat
- SingleCellExperiment
- SpatialExperiment
- SummarizedExperiment
- reticulate (for AnnData)
- testthat
- knitr
- rmarkdown

## Publishing Package

### Option 1: GitHub

```bash
cd ~/PrognosticScorer
git init
git add .
git commit -m "Initial package commit"
git remote add origin https://github.com/username/PrognosticScorer.git
git push -u origin main
```

Users install with:
```r
devtools::install_github("username/PrognosticScorer")
```

### Option 2: Local Distribution

```r
# Build source package
devtools::build()

# This creates: PrognosticScorer_0.1.0.tar.gz
# Share this file, users install with:
install.packages("PrognosticScorer_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Option 3: CRAN (Advanced)

1. Ensure package passes all checks with 0 errors, 0 warnings
2. Add cran-comments.md
3. Submit via `devtools::release()`

## Troubleshooting

### Issue: "could not find function 'get_data'"
**Solution**: Make sure all R files are in the R/ directory and package is properly installed

### Issue: "object 'precog' not found"
**Solution**: Run `usethis::use_data()` to save reference datasets

### Issue: No method for Seurat objects
**Solution**: Add `Seurat` to Suggests in DESCRIPTION

### Issue: Documentation not generating
**Solution**: 
```r
roxygen2::roxygenise()
# Or
devtools::document()
```

### Issue: Namespace errors
**Solution**: Check that all functions have `@export` tag or are marked `@keywords internal`

## Best Practices

1. **Version Control**: Use Git from the start
2. **Documentation**: Write clear examples in roxygen comments
3. **Testing**: Add unit tests in tests/testthat/
4. **Vignettes**: Update introduction.Rmd with real examples
5. **NEWS.md**: Track changes between versions
6. **CITATION**: Add citation file for academic use

## Next Steps

1. Customize DESCRIPTION with your details
2. Add tests for core functions
3. Build vignette with real data examples
4. Add NEWS.md to track versions
5. Consider adding:
   - Additional visualization functions
   - Integration with pathway analysis
   - Batch effect correction options
   - Multi-reference scoring

## Support

For issues or questions:
1. Check package documentation: `?score_expression`
2. View vignette: `vignette("introduction", package = "PrognosticScorer")`
3. Report bugs on GitHub Issues
4. Email: your.email@example.com
