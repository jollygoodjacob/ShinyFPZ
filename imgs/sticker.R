library(ggplot2)
library(dplyr)
library(hexSticker)
library(grid)
library(png)

# -----------------------------
# 1. Stylized river segments
# -----------------------------
river <- tibble::tribble(
  ~x,   ~y,   ~seg,
  0.18, 0.18, "FPZ 1",
  0.24, 0.24, "FPZ 1",
  0.30, 0.29, "FPZ 1",
  0.36, 0.34, "FPZ 2",
  0.42, 0.39, "FPZ 2",
  0.48, 0.45, "FPZ 3",
  0.54, 0.51, "FPZ 3",
  0.60, 0.58, "FPZ 4",
  0.66, 0.65, "FPZ 4",
  0.72, 0.72, "FPZ 5",
  0.78, 0.79, "FPZ 5"
)

# Split into short path chunks so each FPZ can have its own color
river_segments <- river %>%
  mutate(xend = lead(x), yend = lead(y), seg_end = lead(seg)) %>%
  filter(!is.na(xend)) %>%
  mutate(seg_draw = seg)

# -----------------------------
# 2. Simplified dendrogram
# -----------------------------
dendro_segments <- tibble::tribble(
  ~x,   ~y,   ~xend, ~yend,
  0.50, 0.86, 0.50, 0.70,  # trunk
  0.50, 0.70, 0.35, 0.60,
  0.50, 0.70, 0.65, 0.60,
  0.35, 0.60, 0.26, 0.50,
  0.35, 0.60, 0.42, 0.50,
  0.65, 0.60, 0.58, 0.50,
  0.65, 0.60, 0.74, 0.50
)

# Connector from dendrogram to river
connector <- tibble::tribble(
  ~x,   ~y,   ~xend, ~yend,
  0.50, 0.50, 0.48, 0.45
)

# -----------------------------
# 3. Palette
# -----------------------------
fpz_cols <- c(
  "FPZ 1" = "#1B9E77",
  "FPZ 2" = "#2C7FB8",
  "FPZ 3" = "#D95F02",
  "FPZ 4" = "#7570B3",
  "FPZ 5" = "#E6AB02"
)

# -----------------------------
# 4. Build icon graphic
# -----------------------------
icon_plot <- ggplot() +
  # Dendrogram
  geom_segment(
    data = dendro_segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 1.8,
    color = "white",
    lineend = "round"
  ) +
  # Connector
  geom_segment(
    data = connector,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 2.1,
    color = "white",
    lineend = "round"
  ) +
  # River segments
  geom_segment(
    data = river_segments,
    aes(x = x, y = y, xend = xend, yend = yend, color = seg_draw),
    linewidth = 4.2,
    lineend = "round"
  ) +
  # Thin white halo under river for contrast
  geom_segment(
    data = river_segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 5.6,
    color = alpha("white", 0.18),
    lineend = "round"
  ) +
  scale_color_manual(values = fpz_cols) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save transparent icon
ggsave(
  filename = "ShinyFPZ_icon.png",
  plot = icon_plot,
  width = 5,
  height = 5,
  dpi = 300,
  bg = "transparent"
)

# -----------------------------
# 5. Read icon and place in hex sticker
# -----------------------------
img <- png::readPNG("ShinyFPZ_icon.png")

sticker(
  img,
  package = "ShinyFPZ",
  p_size = 7.2,
  p_color = "white",
  p_y = 0.13,
  s_x = 1,
  s_y = 0.58,
  s_width = 0.92,
  s_height = 0.92,
  h_fill = "#062047",
  h_color = "#062047",
  filename = "ShinyFPZ_hex.png"
)