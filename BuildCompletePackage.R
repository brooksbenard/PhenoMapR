################################################################################
#
#  PROGNOSTICSCORER PACKAGE BUILD SCRIPT
#
#  This script reads z-score files and builds the complete package
#
################################################################################

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                       ║\n")
cat("║              PROGNOSTICSCORER PACKAGE BUILD SCRIPT                    ║\n")
cat("║                                                                       ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# CONFIGURATION - EDIT THESE PATHS
# ==============================================================================

# Where your z-score files are located
ZSCORE_DIR <- "~/Downloads"  # Change this to your z-score file location

# Where you want to create the package
PACKAGE_DIR <- "~/Desktop/scIMPEL"  # Change if desired

# Where the package R source files are located (downloaded from outputs)
SOURCE_FILES_DIR <- "~/Downloads"  # Where you downloaded the .R files

# Z-score file names - EDIT THESE to match your actual file names
ZSCORE_FILES <- list(
  precog = "precog_zscores.csv",      # Change to your actual file name
  tcga = "tcga_zscores.csv",          # Change to your actual file name
  pediatric = "pediatric_zscores.csv", # Change to your actual file name
  ici = "ici_zscores.csv"             # Change to your actual file name
)

# File format for z-score files: "csv", "rdata", "rds", or "excel"
ZSCORE_FORMAT <- "csv"  # Change if your files are different format

# ==============================================================================
# STEP 0: INSTALL REQUIRED PACKAGES
# ==============================================================================

cat("\n[Step 0/6] Checking and installing required packages...\n")

required_packages <- c(
  "devtools", "roxygen2", "usethis", "testthat",
  "Matrix", "dplyr", "glue", "progress",
  "Seurat"  # Optional but recommended
)

missing_packages <- required_packages[!sapply(required_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})]

if (length(missing_packages) > 0) {
  cat("  Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
} else {
  cat("  ✓ All required packages already installed\n")
}

library(devtools)
library(usethis)
library(dplyr)

# ==============================================================================
# STEP 1: CREATE PACKAGE STRUCTURE
# ==============================================================================

cat("\n[Step 1/6] Creating package structure...\n")

# Expand paths
PACKAGE_DIR <- path.expand(PACKAGE_DIR)
ZSCORE_DIR <- path.expand(ZSCORE_DIR)
SOURCE_FILES_DIR <- path.expand(SOURCE_FILES_DIR)

# Create package directory
if (!dir.exists(PACKAGE_DIR)) {
  suppressMessages(create_package(PACKAGE_DIR, open = FALSE))
  cat("  ✓ Created package directory:", PACKAGE_DIR, "\n")
} else {
  cat("  → Package directory already exists\n")
}

# Change to package directory
setwd(PACKAGE_DIR)

# Create subdirectories
dirs <- c("R", "data", "man", "tests/testthat", "vignettes")
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("  ✓ Created", d, "\n")
  }
}

# ==============================================================================
# STEP 2: COPY R SOURCE FILES
# ==============================================================================

cat("\n[Step 2/6] Copying R source files...\n")

r_source_files <- c(
  "score_expression.R",
  "data_handlers.R",
  "weighted_sum_scoring.R",
  "prognostic_data.R",
  "utils.R"
)

for (f in r_source_files) {
  src <- file.path(SOURCE_FILES_DIR, f)
  dst <- file.path(PACKAGE_DIR, "R", f)

  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("  ✓ Copied", f, "\n")
  } else {
    cat("  ✗ WARNING: Could not find", f, "in", SOURCE_FILES_DIR, "\n")
    cat("    Please manually copy this file to:", file.path(PACKAGE_DIR, "R"), "\n")
  }
}

# Copy package metadata files
metadata_files <- c("DESCRIPTION", "README.md")
for (f in metadata_files) {
  src <- file.path(SOURCE_FILES_DIR, f)
  dst <- file.path(PACKAGE_DIR, f)

  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("  ✓ Copied", f, "\n")
  }
}

# Copy vignette
vignette_src <- file.path(SOURCE_FILES_DIR, "introduction.Rmd")
vignette_dst <- file.path(PACKAGE_DIR, "vignettes", "introduction.Rmd")
if (file.exists(vignette_src)) {
  file.copy(vignette_src, vignette_dst, overwrite = TRUE)
  cat("  ✓ Copied vignette\n")
}

# ==============================================================================
# STEP 3: LOAD Z-SCORE DATA
# ==============================================================================

cat("\n[Step 3/6] Loading z-score data files...\n")

load_zscore_file <- function(filename, format = "csv") {
  filepath <- file.path(ZSCORE_DIR, filename)

  if (!file.exists(filepath)) {
    cat("  ✗ File not found:", filepath, "\n")
    return(NULL)
  }

  data <- switch(format,
                 "csv" = {
                   read.csv(filepath, row.names = 1, check.names = FALSE)
                 },
                 "rdata" = {
                   env <- new.env()
                   load(filepath, envir = env)
                   # Assume first object in RData file
                   env[[ls(env)[1]]]
                 },
                 "rds" = {
                   readRDS(filepath)
                 },
                 "excel" = {
                   if (!requireNamespace("readxl", quietly = TRUE)) {
                     install.packages("readxl")
                   }
                   library(readxl)
                   df <- as.data.frame(read_excel(filepath))
                   rownames(df) <- df[,1]
                   df[,-1]
                 },
                 stop("Unknown format: ", format)
  )

  # Convert to data.frame if needed
  data <- as.data.frame(data)

  # Verify structure
  if (is.null(rownames(data)) || all(rownames(data) == as.character(1:nrow(data)))) {
    warning("Data may not have gene names as rownames!")
  }

  return(data)
}

# Load each z-score dataset
precog <- load_zscore_file(ZSCORE_FILES$precog, ZSCORE_FORMAT)
tcga <- load_zscore_file(ZSCORE_FILES$tcga, ZSCORE_FORMAT)
pediatric <- load_zscore_file(ZSCORE_FILES$pediatric, ZSCORE_FORMAT)
ici <- load_zscore_file(ZSCORE_FILES$ici, ZSCORE_FORMAT)

# Check what was loaded
if (!is.null(precog)) {
  cat("  ✓ Loaded PRECOG:", nrow(precog), "genes x", ncol(precog), "cancer types\n")
  cat("    Cancer types:", paste(head(colnames(precog), 5), collapse = ", "), "...\n")
}

if (!is.null(tcga)) {
  cat("  ✓ Loaded TCGA:", nrow(tcga), "genes x", ncol(tcga), "cancer types\n")
  cat("    Cancer types:", paste(head(colnames(tcga), 5), collapse = ", "), "...\n")
}

if (!is.null(pediatric)) {
  cat("  ✓ Loaded Pediatric:", nrow(pediatric), "genes x", ncol(pediatric), "cancer types\n")
  cat("    Cancer types:", paste(head(colnames(pediatric), 5), collapse = ", "), "...\n")
}

if (!is.null(ici)) {
  cat("  ✓ Loaded ICI:", nrow(ici), "genes x", ncol(ici), "cancer types\n")
  cat("    Cancer types:", paste(head(colnames(ici), 3), collapse = ", "), "...\n")
}

# ==============================================================================
# STEP 4: CREATE DATASETS_INFO MAPPING
# ==============================================================================

cat("\n[Step 4/6] Creating datasets_info mapping...\n")

# Option 1: Create minimal empty datasets_info
datasets_info <- data.frame(
  dataset_name = character(),
  precog_label = character(),
  tcga_label = character(),
  pediatric_precog_label = character(),
  ici_precog_label = character(),
  stringsAsFactors = FALSE
)

cat("  ✓ Created empty datasets_info (you can populate this later if needed)\n")
cat("    To add dataset mappings, edit the datasets_info object before saving\n")

# Option 2: Create example mappings (uncomment and edit if you want)
# datasets_info <- data.frame(
#   dataset_name = c("my_breast_study", "my_lung_study"),
#   precog_label = c("BRCA", "LUAD"),
#   tcga_label = c("BRCA", "LUAD"),
#   pediatric_precog_label = c(NA, NA),
#   ici_precog_label = c(NA, "NSCLC"),
#   stringsAsFactors = FALSE
# )

# ==============================================================================
# STEP 5: SAVE DATA TO PACKAGE
# ==============================================================================

cat("\n[Step 5/6] Saving data to package...\n")

# Save each dataset that was successfully loaded
if (!is.null(precog)) {
  use_data(precog, overwrite = TRUE, compress = "xz")
  cat("  ✓ Saved precog.rda\n")
}

if (!is.null(tcga)) {
  use_data(tcga, overwrite = TRUE, compress = "xz")
  cat("  ✓ Saved tcga.rda\n")
}

if (!is.null(pediatric)) {
  use_data(pediatric, overwrite = TRUE, compress = "xz")
  cat("  ✓ Saved pediatric.rda\n")
}

if (!is.null(ici)) {
  use_data(ici, overwrite = TRUE, compress = "xz")
  cat("  ✓ Saved ici.rda\n")
}

use_data(datasets_info, overwrite = TRUE, compress = "xz")
cat("  ✓ Saved datasets_info.rda\n")

# Verify data files were created
data_files <- list.files("data", pattern = "\\.rda$")
cat("\n  Data files created:", paste(data_files, collapse = ", "), "\n")

# ==============================================================================
# STEP 6: BUILD PACKAGE
# ==============================================================================

cat("\n[Step 6/6] Building package...\n")

# Setup package infrastructure
cat("  Setting up package infrastructure...\n")
if (!file.exists("LICENSE")) {
  use_mit_license("Your Name")  # Edit this with your name
}
use_roxygen_md()
use_build_ignore(c("setup_data.R", "build_package.R"))

# Generate documentation
cat("  Generating documentation...\n")
document()
cat("  ✓ Documentation generated\n")

# Check package
cat("\n  Checking package (this may take a minute)...\n")
check_results <- check(quiet = TRUE)

if (length(check_results$errors) > 0) {
  cat("  ✗ ERRORS found:\n")
  print(check_results$errors)
} else {
  cat("  ✓ No errors\n")
}

if (length(check_results$warnings) > 0) {
  cat("  ⚠ WARNINGS:\n")
  print(check_results$warnings)
} else {
  cat("  ✓ No warnings\n")
}

# Install package
cat("\n  Installing package...\n")
install(upgrade = "never", quiet = TRUE)
cat("  ✓ Package installed!\n")

# ==============================================================================
# COMPLETION SUMMARY
# ==============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                       ║\n")
cat("║                     PACKAGE BUILD COMPLETE!                           ║\n")
cat("║                                                                       ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Package location:", PACKAGE_DIR, "\n\n")

cat("NEXT STEPS:\n")
cat("  1. Test the package:\n")
cat("     library(PrognosticScorer)\n")
cat("     list_cancer_types('precog')\n\n")

cat("  2. Edit DESCRIPTION file with your information:\n")
cat("     file.edit('", file.path(PACKAGE_DIR, 'DESCRIPTION'), "')\n\n", sep = "")

cat("  3. Test with your data:\n")
cat("     scores <- score_expression(\n")
cat("       expression = your_seurat_object,\n")
cat("       reference = 'tcga',\n")
cat("       cancer_type = 'BRCA'\n")
cat("     )\n\n")

cat("  4. Share on GitHub:\n")
cat("     usethis::use_git()\n")
cat("     usethis::use_github()\n\n")

# Quick test
cat("QUICK TEST:\n")
tryCatch({
  library(PrognosticScorer)

  if (!is.null(precog) && ncol(precog) > 0) {
    cancer_types <- colnames(precog)[1]
    cat("  Available cancer types in PRECOG:", paste(head(colnames(precog), 10), collapse = ", "), "\n")
    cat("\n  ✓ Package loaded successfully!\n")
  }
}, error = function(e) {
  cat("\n  ✗ Error loading package:", e$message, "\n")
})

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("\nPackage build script completed!\n\n")
