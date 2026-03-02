if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}

library(hexSticker)

hex_logo_path <- path.expand("~/Desktop/scIMPEL/inst/figures/hex_logo.png")

sticker(
  hex_logo_path,
  package = "scIMPEL",
  p_size = 20,
  p_color = "#36454f",
  s_x = 1,
  s_y = 0.85,
  s_width = 0.75,
  s_height = 0.75,
  h_fill = "#fff",
  h_color = "#536878",
  filename = path.expand("~/Desktop/scIMPEL/inst/figures/SCIMPEL_logo.png")
)
