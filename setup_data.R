# Setup Script for PrognosticScorer Package Data
# Run this script to prepare your reference datasets for the package

library(usethis)

# Assuming you have your reference dataframes loaded:
# - precog
# - tcga  
# - pediatric
# - ici
# - datasets_info

# Save data to package data/ directory
# Note: This should be run from the package root directory

save_package_data <- function() {
  
  # Check if objects exist
  required_objects <- c("precog", "tcga", "pediatric", "ici", "datasets_info")
  
  for (obj_name in required_objects) {
    if (!exists(obj_name)) {
      stop(glue::glue("{obj_name} not found in environment. Please load it first."))
    }
  }
  
  # Create data directory if it doesn't exist
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  # Save each dataset
  message("Saving reference datasets...")
  
  # PRECOG
  message("  - precog")
  usethis::use_data(precog, overwrite = TRUE, compress = "xz")
  
  # TCGA
  message("  - tcga")
  usethis::use_data(tcga, overwrite = TRUE, compress = "xz")
  
  # Pediatric
  message("  - pediatric")
  usethis::use_data(pediatric, overwrite = TRUE, compress = "xz")
  
  # ICI
  message("  - ici")
  usethis::use_data(ici, overwrite = TRUE, compress = "xz")
  
  # Dataset info
  message("  - datasets_info")
  usethis::use_data(datasets_info, overwrite = TRUE, compress = "xz")
  
  message("Data saved successfully!")
  message("\nNext steps:")
  message("1. Document your data with roxygen2 comments in R/prognostic_data.R")
  message("2. Run devtools::document() to generate documentation")
  message("3. Run devtools::check() to validate package")
}


# Alternative: Save data manually
save_data_manual <- function() {
  
  # Ensure data directory exists
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  # Save as .rda files
  save(precog, file = "data/precog.rda", compress = "xz")
  save(tcga, file = "data/tcga.rda", compress = "xz")
  save(pediatric, file = "data/pediatric.rda", compress = "xz")
  save(ici, file = "data/ici.rda", compress = "xz")
  save(datasets_info, file = "data/datasets_info.rda", compress = "xz")
  
  message("Data saved to data/ directory")
}


# Example: Create dummy datasets_info if you don't have one
create_example_datasets_info <- function() {
  
  datasets_info <- data.frame(
    dataset_name = c("example_breast", "example_lung", "example_melanoma"),
    precog_label = c("BRCA", "LUAD", "SKCM"),
    tcga_label = c("BRCA", "LUAD", "SKCM"),
    pediatric_precog_label = c(NA, NA, NA),
    ici_precog_label = c(NA, "NSCLC", "MELANOMA_Metastatic"),
    stringsAsFactors = FALSE
  )
  
  return(datasets_info)
}


# Instructions for preparing your data
prepare_reference_data <- function() {
  
  cat("
=================================================================
PREPARING REFERENCE DATA FOR PROGNOSTICSCORER PACKAGE
=================================================================

1. LOAD YOUR REFERENCE DATASETS
   Make sure the following objects are in your R environment:
   
   - precog: data.frame with genes as rownames, cancer types as columns
   - tcga: data.frame with genes as rownames, cancer types as columns  
   - pediatric: data.frame with genes as rownames, cancer types as columns
   - ici: data.frame with genes as rownames, ICI columns
   - datasets_info: mapping of dataset names to cancer type labels

2. VALIDATE DATA FORMAT
   
   Each reference dataframe should have:
   - Genes as rownames
   - Cancer types as column names
   - Z-scores as values
   - No missing rownames
   
   Example:
   
                    BRCA    LUAD    COAD
   TP53             3.25   -2.10    1.45
   MYC             -1.80    2.35   -0.95
   EGFR             0.50    4.20    0.30

3. CREATE datasets_info
   
   This dataframe maps your dataset names to cancer type labels:
   
   datasets_info <- data.frame(
     dataset_name = c('my_breast_study', 'my_lung_study'),
     precog_label = c('BRCA', 'LUAD'),
     tcga_label = c('BRCA', 'LUAD'),
     pediatric_precog_label = c(NA, NA),
     ici_precog_label = c(NA, 'NSCLC')
   )

4. SAVE DATA TO PACKAGE
   
   From the package root directory:
   
   source('setup_data.R')
   save_package_data()
   
   Or manually:
   
   save_data_manual()

5. BUILD PACKAGE
   
   devtools::document()
   devtools::check()
   devtools::install()

=================================================================
")
  
}

# Run this to see instructions
prepare_reference_data()
