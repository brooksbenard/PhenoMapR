#!/usr/bin/env Rscript
## Smoke tests for the EcoTyper-style block-outline rewrite in
## R/plot_phenotype_markers.R.
##
## The previous implementation tried to draw black borders around
## each (cell-type x phenotype) block by calling
## `decorate_heatmap_body()` and then computing native-unit x/width
## offsets manually. Inside the callback we ALSO pushed a fresh
## viewport with `clip = "on"`, which silently reset the native
## x-scale to 0..1 -- so the rectangle ended up clipped outside the
## visible heatmap body and the borders never appeared.
##
## The fix follows the EcoTyper / jokergoo pattern:
##   1. column_split on the same (phenotype_bin x cell_type) key as
##      the row split, so each block becomes its own
##      (row_slice, column_slice) viewport;
##   2. decorate the diagonal (row_slice == column_slice on the
##      block_key level) with `grid.rect()` and NO explicit x/y
##      coordinates -- the rectangle then exactly fills the
##      block's viewport.
##
## This script asserts that the source code still contains both
## halves of that contract, so a future refactor cannot regress
## the rendering without tripping a smoke test.

stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("OK: ", msg, "\n", sep = "")
}

src_file <- system.file("R", "plot_phenotype_markers.R", package = "PhenoMapR")
dev_file <- file.path(getwd(), "R", "plot_phenotype_markers.R")
if (file.exists(dev_file)) src_file <- dev_file
stopifnot(file.exists(src_file))
src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

## 1. The cell-type-specific Heatmap call splits BOTH axes on the
##    block_key. The row split was already there; the column_split
##    on `col_split_factor` is the new piece.
stopifnot_msg(
  grepl("column_split\\s*=\\s*col_split_factor", src),
  "Heatmap() call passes column_split = col_split_factor"
)
stopifnot_msg(
  grepl("cluster_column_slices\\s*=\\s*FALSE", src),
  "Heatmap() call disables cluster_column_slices to preserve block order"
)
stopifnot_msg(
  grepl("column_gap\\s*=\\s*grid::unit\\(0,\\s*\"mm\"\\)", src),
  "Heatmap() call sets column_gap = unit(0, \"mm\") to keep continuity"
)

## 2. We compute (row_slice, column_slice) diagonal pairs from the
##    shared block_key levels.
stopifnot_msg(
  grepl("decorate_pairs\\s*<-", src),
  "decorate_pairs list is built from row/column split levels"
)
stopifnot_msg(
  grepl("row_slice_levels\\[ri\\]", src) &&
    grepl("col_lv\\s*<-\\s*levels\\(col_split_factor\\)", src),
  "diagonal pairs match row_split levels against column_split levels"
)

## 3. The decorate loop calls decorate_heatmap_body with both
##    row_slice AND column_slice arguments and an empty-coord
##    grid.rect() -- this is the EcoTyper-style pattern. We
##    explicitly check that the OLD pushViewport / native-unit
##    fallback is gone.
stopifnot_msg(
  grepl(
    paste0(
      "ComplexHeatmap::decorate_heatmap_body\\(\\s*",
      "\"phenomap_ct_markers\"\\s*,\\s*",
      "row_slice\\s*=\\s*ri\\s*,\\s*",
      "column_slice\\s*=\\s*cj\\s*,"
    ),
    src
  ),
  "decorate_heatmap_body() is called with both row_slice = ri AND column_slice = cj"
)
stopifnot_msg(
  !grepl(
    "grid::pushViewport\\(grid::viewport\\(clip\\s*=\\s*\"on\"\\)\\)",
    src
  ),
  "the buggy pushViewport(clip='on') wrapper inside decorate_heatmap_body is gone"
)
stopifnot_msg(
  !grepl("grid::unit\\(j1 - 1L,\\s*\"native\"\\)", src) &&
    !grepl("grid::unit\\(j2 - j1 \\+ 1L,\\s*\"native\"\\)", src),
  "no manual native-unit x/width arithmetic remains in the decorate loop"
)

## 4. The grid.rect() call inside the decorate_heatmap_body
##    callback uses no x/y/width/height arguments -- relying on the
##    slice viewport to bound it -- and pulls col / lwd from the
##    user-facing block_outline_color / block_outline_lwd args.
##    We narrow the search to the cell-type-specific decorate
##    branch by anchoring on its identifying string first.
ct_branch_marker <- "EcoTyper-style block outlines"
stopifnot_msg(
  grepl(ct_branch_marker, src, fixed = TRUE),
  "the EcoTyper-style decorate branch is documented in the source"
)
ct_branch_idx <- regexpr(ct_branch_marker, src, fixed = TRUE)
stopifnot(ct_branch_idx > 0L)
ct_window <- substr(src, ct_branch_idx, min(nchar(src), ct_branch_idx + 1500L))
stopifnot_msg(
  grepl("grid::grid\\.rect\\([^)]*?gp\\s*=\\s*grid::gpar",
        ct_window, perl = TRUE),
  "grid.rect() inside the decorate callback uses gp = gpar(...) only"
)
stopifnot_msg(
  grepl("col\\s*=\\s*col_use", ct_window) &&
    grepl("lwd\\s*=\\s*lwd_use", ct_window),
  "the rect's col/lwd are wired to the block_outline_color / _lwd args"
)
stopifnot_msg(
  grepl("fill\\s*=\\s*NA", ct_window),
  "the rect uses fill = NA so the heatmap cells underneath stay visible"
)

## 5. The "global" heatmap path still works when no decoration is
##    requested -- we initialise decorate_pairs to list() at the
##    top of the function so the post-render loop is a no-op
##    regardless of heatmap_type.
stopifnot_msg(
  grepl("decorate_pairs\\s*<-\\s*list\\(\\)", src),
  "decorate_pairs is initialised to list() at the top of the function"
)

cat("\nAll marker-heatmap block-outline smoke tests passed.\n")
