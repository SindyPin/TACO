library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Users/sindypinero/Downloads/")

df <- read.csv("DCE_results_1.csv", header = TRUE, stringsAsFactors = FALSE)

df <- df %>%
  filter(predictor != "#N/A", outcome != "#N/A")

df_top <- df %>%
  mutate(
    abs_dce    = abs(DCE),
    path_label = paste(predictor, outcome, sep = "->")
  ) %>%
  arrange(desc(abs_dce))

dot_plot <- ggplot(df_top, aes(x = p_value, y = abs_dce)) +
  geom_point(aes(color = DCE), size = 2) +
  geom_text_repel(
    data = subset(df_top, adj_p_value < 0.05 & abs_dce > 0.41),
    aes(label = path_label),
    size = 3,
    max.overlaps = Inf,
    force = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.size = 0.3,
    segment.color = "grey50",
    seed = 42
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "DCE") +
  scale_x_continuous(trans = "log10") +
  labs(x = "p-value (log scale)", y = "Absolute DCE") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 10))

print(dot_plot)

ggsave(
  filename = "Gene_pairs_DCE_dotplot.png",
  plot     = dot_plot,
  device   = "png",
  width    = 10,
  height   = 8,
  dpi      = 300
)

message("Saved to: ", getwd())