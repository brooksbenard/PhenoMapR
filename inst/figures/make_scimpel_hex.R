if (!requireNamespace("hexSticker", quietly = TRUE)) {
  install.packages("hexSticker")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(hexSticker)
library(ggplot2)

# Simple hex logo for scIMPEL
p <- ggplot() +
  geom_text(
    aes(0, 0, label = "scIMPEL"),
    size = 8,
    fontface = "bold"
  ) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  theme_void()

sticker(
  p,
  package = "scIMPEL",
  p_size = 8,
  s_x = 1,
  s_y = 0.8,
  s_width = 1.4,
  s_height = 1.2,
  h_fill = "#1b2838",
  h_color = "#66c2a5",
  filename = "inst/figures/SCIMPEL_logo.png"
)

