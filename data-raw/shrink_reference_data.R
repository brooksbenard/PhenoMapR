# Shrink reference z-score RDA files for smaller package size and faster load.
# Scoring uses only genes with |z-score| > z_score_cutoff (default 2).
# This script:
#   1. Sets z-scores with |z| <= cutoff to NA (smaller when compressed).
#   2. Drops genes that are NA across all columns.
# Run from package root: Rscript data-raw/shrink_reference_data.R
# Requires devtools or usethis for use_data (or use save() as below).

Z_CUTOFF <- 2  # match default z_score_cutoff in score_expression()

shrink_zscore_table <- function(x, cutoff = Z_CUTOFF) {
  x <- as.matrix(x)
  x[!is.na(x) & abs(x) <= cutoff] <- NA
  all_na <- rowSums(!is.na(x)) == 0
  x[!all_na, , drop = FALSE]
}

run_shrink <- function() {
  data_dir <- "data"
  if (!dir.exists(data_dir)) stop("Run from package root (directory 'data' not found).")

  # Reference datasets only (not datasets_info)
  ref_names <- c("precog", "tcga", "pediatric", "ici")

  for (nm in ref_names) {
    path <- file.path(data_dir, paste0(nm, ".rda"))
    if (!file.exists(path)) {
      message("  Skip ", nm, " (file not found)")
      next
    }
    load(path)
    obj <- get(nm)
    dim_before <- dim(obj)
    size_before <- file.size(path)

    obj <- shrink_zscore_table(obj, Z_CUTOFF)
    dim_after <- dim(obj)

    # Save with same object name so load() still works
    assign(nm, obj)
    save(list = nm, file = path, compress = "xz")

    size_after <- file.size(path)
    message(sprintf(
      "  %s: %d x %d -> %d x %d rows; %.1f MB -> %.1f MB",
      nm, dim_before[1], dim_before[2], dim_after[1], dim_after[2],
      size_before / 1e6, size_after / 1e6
    ))
  }
  message("Done.")
}

run_shrink()
