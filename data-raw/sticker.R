library(hexSticker)
library(showtext)
font_add_google("Gochi Hand", "gochi")
showtext_auto()
img <- "sticker.png"
sticker(
  img, package = "ASP", p_size = 8, s_x = 1, s_y = .75, s_width = .6, s_height = .5,
  h_size = .7, h_fill = "#FFE4E1", p_color = "#FA8072", h_color = "#F08080", p_family = "gochi",
  filename = "sticker.svg"
)
