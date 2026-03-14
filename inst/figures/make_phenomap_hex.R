# Build the hex sticker for PhenoMapR. Run from package root:
#   Rscript inst/figures/make_phenomap_hex.R
# Output: inst/figures/PhenoMapR_logo.png (hex shape + "PhenoMapR" text + map image).
# Use that file (or man/figures/logo.png copied from it) for README and pkgdown.

if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker", repos = "https://cloud.r-project.org")
}

# Ensure we run from package root (parent of inst/)
if (basename(getwd()) == "figures" && dir.exists("../../inst")) {
  setwd("../..")
} else if (!dir.exists("inst/figures")) {
  if (file.exists("make_phenomap_hex.R")) setwd(dirname(dirname(normalizePath("make_phenomap_hex.R"))))
}

library(hexSticker)

hex_logo_path <- file.path("inst", "figures", "PhenoMapR_map_logo_2.png")
if (!file.exists(hex_logo_path)) stop("Source image not found: ", hex_logo_path)

sticker(
  hex_logo_path,
  package = "PhenoMapR",
  p_size = 20,
  p_color = "#36454f",
  s_x = 1,
  s_y = 0.72,
  s_width = 0.6,
  s_height = 0.6,
  h_fill = "#fff",
  h_color = "#536878",
  filename = file.path("inst", "figures", "PhenoMapR_logo.png")
)
message("Hex sticker written to inst/figures/PhenoMapR_logo.png. Copy to man/figures/logo.png for README/pkgdown.")
