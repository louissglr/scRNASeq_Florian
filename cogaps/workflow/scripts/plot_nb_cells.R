suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tibble)
  library(forcats)
  library(RColorBrewer)
  library(scales)
})

contrib_file <- snakemake@input[["contrib"]] 
metadata_file <- snakemake@input[["metadata"]]
output_png  <- snakemake@output[["png"]]
sample_name <- snakemake@wildcards[["sample"]]

contrib <- read.table(contrib_file, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

metadata <- metadata |>
  mutate(group = paste(orig_tumor, tnt, Sting_status, sep = "_"))


contrib <- contrib |>
  filter(barcode %in% metadata$cell_id)

df <- contrib |>
  inner_join(
    metadata |>
      select(cell_id, group),
    by = c("barcode" = "cell_id")
  )

df_prop <- df |>
  group_by(group, dominant_pattern) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(group) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

df_prop$dominant_pattern <- factor(df_prop$dominant_pattern)

df_prop <- df_prop |>
  mutate(group = factor(group, levels = sort(unique(group))))

df_counts <- df |>
  count(group, name = "n_cells") |>
  mutate(group = factor(group, levels = levels(df_prop$group)))

p <- ggplot(df_prop, aes(x = group, y = prop, fill = dominant_pattern)) +
  geom_col(width = 0.8, color = "white") +
  
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1),
    expand = c(0, 0)
  ) +
  
  geom_text(
    data = df_counts,
    aes(
      x = group,
      y = 1.02,
      label = paste0("n = ", n_cells)
    ),
    inherit.aes = FALSE,
    vjust = 0,  
    size = 2.5
  ) +
  
  labs(
    x = "Group (orig_tumor_tnt_Sting_status)",
    y = "Proportion of cells",
    fill = "Dominant pattern"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "right"
  )

png(output_png, width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

