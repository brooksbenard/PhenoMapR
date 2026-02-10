#!/usr/bin/env Rscript

# Complete Package Build Script for PrognosticScorer
# This script guides you through building the entire package

cat("
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║           PROGNOSTICSCORER PACKAGE BUILD SCRIPT              ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝

")

# Check required packages
required_pkgs <- c("devtools", "roxygen2", "usethis", "testthat")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("\n[!] Missing required packages:", paste(missing_pkgs, collapse = ", "), "\n")
  cat("    Install with: install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))\n\n")
  stop("Please install missing packages first")
}

library(devtools)
library(usethis)
library(roxygen2)

# ==============================================================================
# STEP 1: CREATE PACKAGE STRUCTURE
# ==============================================================================

create_package_structure <- function(package_path = "~/PrognosticScorer") {
  
  cat("\n[1/7] Creating package structure...\n")
  
  if (!dir.exists(package_path)) {
    usethis::create_package(package_path, open = FALSE)
    cat("    ✓ Package directory created at:", package_path, "\n")
  } else {
    cat("    → Package directory already exists\n")
  }
  
  setwd(package_path)
  
  # Create subdirectories
  dirs <- c("R", "data", "man", "tests/testthat", "vignettes")
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      cat("    ✓ Created", d, "\n")
    }
  }
  
  return(package_path)
}

# ==============================================================================
# STEP 2: COPY R FILES
# ==============================================================================

copy_r_files <- function() {
  
  cat("\n[2/7] Setting up R source files...\n")
  
  r_files <- c(
    "score_expression.R",
    "data_handlers.R", 
    "weighted_sum_scoring.R",
    "prognostic_data.R",
    "utils.R"
  )
  
  cat("    Please ensure these files are in the R/ directory:\n")
  for (f in r_files) {
    cat("      -", f, "\n")
  }
  
  cat("    Copy them manually or run:\n")
  cat("    cp /path/to/files/*.R R/\n")
}

# ==============================================================================
# STEP 3: SETUP PACKAGE METADATA
# ==============================================================================

setup_metadata <- function() {
  
  cat("\n[3/7] Setting up package metadata...\n")
  
  # License
  if (!file.exists("LICENSE")) {
    usethis::use_mit_license("Your Name")
    cat("    ✓ MIT license added\n")
  }
  
  # README
  if (!file.exists("README.md")) {
    usethis::use_readme_md(open = FALSE)
    cat("    ✓ README.md template created\n")
  }
  
  # .Rbuildignore entries
  usethis::use_build_ignore(c("setup_data.R", "build_package.R"))
  
  # Roxygen
  usethis::use_roxygen_md()
  cat("    ✓ Configured roxygen2\n")
  
  # Testing
  if (!file.exists("tests/testthat.R")) {
    usethis::use_testthat()
    cat("    ✓ Test infrastructure added\n")
  }
}

# ==============================================================================
# STEP 4: PREPARE DATA
# ==============================================================================

prepare_data <- function() {
  
  cat("\n[4/7] Preparing package data...\n")
  cat("    
    IMPORTANT: You need to have the following objects in your environment:
    - precog
    - tcga
    - pediatric
    - ici
    - datasets_info
    
    If you have these loaded, they will be saved to the package.
    If not, you'll need to load them and run save_package_data() manually.
    \n")
  
  # Check if data objects exist
  data_objects <- c("precog", "tcga", "pediatric", "ici", "datasets_info")
  existing <- sapply(data_objects, exists)
  
  if (all(existing)) {
    cat("    ✓ All data objects found in environment\n")
    
    response <- readline("    Save data to package now? (y/n): ")
    if (tolower(response) == "y") {
      for (obj in data_objects) {
        usethis::use_data(get(obj), overwrite = TRUE, compress = "xz")
        cat("      ✓ Saved", obj, "\n")
      }
    }
  } else {
    cat("    ✗ Missing data objects:", paste(data_objects[!existing], collapse = ", "), "\n")
    cat("    → You'll need to run setup_data.R later\n")
  }
}

# ==============================================================================
# STEP 5: GENERATE DOCUMENTATION
# ==============================================================================

generate_docs <- function() {
  
  cat("\n[5/7] Generating documentation...\n")
  
  tryCatch({
    roxygen2::roxygenise()
    cat("    ✓ Documentation generated\n")
  }, error = function(e) {
    cat("    ✗ Error generating docs:", e$message, "\n")
    cat("    → Run devtools::document() manually after fixing R files\n")
  })
}

# ==============================================================================
# STEP 6: CHECK PACKAGE
# ==============================================================================

check_package <- function() {
  
  cat("\n[6/7] Checking package...\n")
  
  response <- readline("    Run R CMD check? This may take a few minutes (y/n): ")
  
  if (tolower(response) == "y") {
    tryCatch({
      check_results <- devtools::check()
      cat("    ✓ Package check complete\n")
      
      if (length(check_results$errors) > 0) {
        cat("    ✗ Errors found:\n")
        print(check_results$errors)
      }
      
      if (length(check_results$warnings) > 0) {
        cat("    ⚠ Warnings:\n")
        print(check_results$warnings)
      }
      
    }, error = function(e) {
      cat("    ✗ Check failed:", e$message, "\n")
    })
  } else {
    cat("    → Skipped package check\n")
  }
}

# ==============================================================================
# STEP 7: BUILD & INSTALL
# ==============================================================================

build_and_install <- function() {
  
  cat("\n[7/7] Building and installing package...\n")
  
  response <- readline("    Build and install package? (y/n): ")
  
  if (tolower(response) == "y") {
    tryCatch({
      devtools::install()
      cat("    ✓ Package installed successfully\n")
      cat("    → You can now use: library(PrognosticScorer)\n")
    }, error = function(e) {
      cat("    ✗ Installation failed:", e$message, "\n")
    })
  } else {
    cat("    → Skipped installation\n")
  }
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main <- function() {
  
  cat("\nThis script will help you build the PrognosticScorer package.\n")
  cat("Make sure you have all R source files ready.\n\n")
  
  response <- readline("Continue? (y/n): ")
  
  if (tolower(response) != "y") {
    cat("\nExiting...\n")
    return(invisible())
  }
  
  # Get package path
  pkg_path <- readline("\nEnter package path [~/PrognosticScorer]: ")
  if (pkg_path == "") pkg_path <- "~/PrognosticScorer"
  pkg_path <- path.expand(pkg_path)
  
  # Run build steps
  tryCatch({
    create_package_structure(pkg_path)
    copy_r_files()
    setup_metadata()
    prepare_data()
    generate_docs()
    check_package()
    build_and_install()
    
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════════╗\n")
    cat("║                     BUILD COMPLETE!                           ║\n")
    cat("╚═══════════════════════════════════════════════════════════════╝\n")
    cat("\n")
    cat("Next steps:\n")
    cat("  1. Edit DESCRIPTION file with your information\n")
    cat("  2. Update README.md with your examples\n")
    cat("  3. Add tests to tests/testthat/\n")
    cat("  4. Build vignette: devtools::build_vignettes()\n")
    cat("  5. Share on GitHub!\n\n")
    
  }, error = function(e) {
    cat("\n✗ Build process failed:", e$message, "\n")
    cat("  Check the errors above and try again.\n\n")
  })
}

# Run if called as script
if (!interactive()) {
  main()
} else {
  cat("\nTo build the package, run: source('build_package.R'); main()\n\n")
}
