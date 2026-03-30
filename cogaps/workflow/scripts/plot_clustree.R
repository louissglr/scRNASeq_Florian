library(clustree)
library(readr)
library(dplyr)
library(ggplot2)

input_file <- snakemake@input[["table"]]
output_file <- snakemake@output[["pdf"]]

clusters_df <- read_tsv(input_file, show_col_types = FALSE)

clusters_df <- clusters_df |>
  mutate(across(starts_with("k"), as.factor))

# plot
p <- clustree(clusters_df, prefix = "k",node_label = "size")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = output_file,
  plot = p,
  width = 16,
  height = 8
)
