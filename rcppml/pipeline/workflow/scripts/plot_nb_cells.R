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
output_pdf  <- snakemake@output[["pdf"]]
sample_name <- snakemake@wildcards[["sample"]]

# Lecture des fichiers
contrib <- read.table(contrib_file, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Créer colonne combinée dans metadata
metadata <- metadata |>
  mutate(group = paste(orig_tumor, tnt, Sting_status, sep = "_"))

# Filtrer contrib pour ne garder que les barcodes présents dans metadata
contrib <- contrib |>
  filter(barcode %in% metadata$cell_id)

# Fusion contrib + metadata
df <- contrib |>
  inner_join(
    metadata |>
      select(cell_id, group),
    by = c("barcode" = "cell_id")
  )

# Calcul des proportions de dominant_program par groupe
df_prop <- df |>
  group_by(group, dominant_program) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(group) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

# Ordre des patterns
pattern_levels <- unique(df_prop$dominant_program)
df_prop <- df_prop |>
  mutate(dominant_program = factor(dominant_program, levels = pattern_levels))

# Ordre des groupes
df_prop <- df_prop |>
  mutate(group = factor(group, levels = sort(unique(group))))

# Nombre de cellules par groupe pour annotation
df_counts <- df |>
  count(group, name = "n_cells") |>
  mutate(group = factor(group, levels = levels(df_prop$group)))

# --------------------------------------------------
# Plot (vertical + couleurs automatiques)
# --------------------------------------------------
p <- ggplot(df_prop, aes(x = group, y = prop, fill = dominant_program)) +
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
      y = 1.05,
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

pdf(output_pdf, width = 10, height = 8)
print(p)
dev.off()