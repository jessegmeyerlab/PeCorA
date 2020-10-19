library("hexSticker")
library("magick")
library(here)

p_size <- 7
p_y <- 1.7
p_color <- "#FFFFFFDD"
h_color <- "black"
h_fill <-  "seagreen"
s_x <- 0.94
s_y <- 0.87
s_width <- 0.7
s_height <- 0.7

sticker(here("Sheep.png"), package="PeCorA",
        p_color = p_color,
        p_size = p_size,
        p_y = p_y,
        h_fill = h_fill,
        h_color = h_color,
        p_family = "RobotoCondensed-Regular",
        s_width = s_width, s_height = s_height,
        s_x = s_x, s_y = s_y,
        filename="PECORA_hex.png")

