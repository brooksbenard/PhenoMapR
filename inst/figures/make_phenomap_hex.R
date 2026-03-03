if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}

library(hexSticker)

hex_logo_path <- file.path("inst", "figures", "PhenoMap_waddington.png")

sticker(
  hex_logo_path,
  package = "PhenoMap",
  p_size = 20,
  p_color = "#36454f",
  s_x = 1,
  s_y = 0.85,
  s_width = 0.75,
  s_height = 0.75,
  h_fill = "#fff",
  h_color = "#536878",
  filename = file.path("inst", "figures", "PhenoMap_logo.png")
)
