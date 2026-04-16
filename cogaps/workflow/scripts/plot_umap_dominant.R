suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

# --------------------------------------------------
# Inputs from snakemake
# --------------------------------------------------

contrib_file <- snakemake@input[["contributions"]]
umap_file    <- snakemake@input[["umap"]]
output_png   <- snakemake@output[["png"]]

sample_name  <- snakemake@wildcards[["sample"]]
dataset      <- snakemake@params[["dataset"]]


contrib_df <- read_tsv(contrib_file, guess_max = 1e6)

colnames(contrib_df) <- tolower(colnames(contrib_df))
colnames(contrib_df) <- gsub("[^a-z0-9_]", "_", colnames(contrib_df))

barcode_candidates <- c("barcode", "barcodes", "cell", "cells", "cell_id", "cellid")
barcode_col <- intersect(colnames(contrib_df), barcode_candidates)[1]

exclude_cols <- c("sample_name", "dominant_pattern", barcode_col)
pattern_cols <- setdiff(colnames(contrib_df), exclude_cols)

clean_umap <- function(file) {

  df <- read_tsv(file, guess_max = 1e6)
  colnames(df) <- tolower(colnames(df))
  colnames(df) <- gsub("[^a-z0-9_]", "_", colnames(df))

  barcode_candidates <- c(
    "barcode", "barcodes",
    "cell", "cells",
    "cell_id", "cellid"
  )

  barcode_col <- intersect(colnames(df), barcode_candidates)[1]
  umap_candidates <- colnames(df)[grepl("umap|dim|coord|x$|y$", colnames(df))]

  umap_x <- umap_candidates[1]
  umap_y <- umap_candidates[2]

  df_clean <- df |>
    rename(
      barcode = all_of(barcode_col),
      UMAP_1  = all_of(umap_x),
      UMAP_2  = all_of(umap_y)
    ) |>
    select(barcode, UMAP_1, UMAP_2, everything())

  return(df_clean)
}

umap_df <- clean_umap(umap_file)


plot_df <- umap_df |>
  left_join(
    contrib_df |> select(all_of(barcode_col), dominant_pattern),
    by = c("barcode" = barcode_col)
  ) |>
  filter(!is.na(dominant_pattern))


plot_df$dominant_pattern <- factor(plot_df$dominant_pattern)


x_min <- min(plot_df$UMAP_1, na.rm = TRUE)
x_max <- max(plot_df$UMAP_1, na.rm = TRUE)
y_min <- min(plot_df$UMAP_2, na.rm = TRUE)
y_max <- max(plot_df$UMAP_2, na.rm = TRUE)

x_offset <- 0.05 * (x_max - x_min)
y_offset <- 0.05 * (y_max - y_min)
arrow_length <- 1.5

p <- ggplot(plot_df,
            aes(x = UMAP_1,
                y = UMAP_2,
                color = dominant_pattern)) +
  geom_point(size = 0.8, alpha = 0.9) +
  theme_void() +
  coord_fixed() +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  ) +
  geom_segment(x = x_min - x_offset, y = y_min - y_offset,
               xend = x_min - x_offset + arrow_length, yend = y_min - y_offset,
               colour = "black", size = 0.5,
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = x_min - x_offset, y = y_min - y_offset,
               xend = x_min - x_offset, yend = y_min - y_offset + arrow_length,
               colour = "black", size = 0.5,
               arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text",
           x = x_min - x_offset + arrow_length/2,
           y = y_min - y_offset - 0.5,
           label = "UMAP1",
           size = 3) +
  annotate("text",
           x = x_min - x_offset - 0.5,
           y = y_min - y_offset + arrow_length/2,
           label = "UMAP2",
           size = 3,
           angle = 90) +
  labs(
    title = paste0("UMAP - Dominant Pattern - ", sample_name, " (", dataset, ")"),
    color = "Dominant pattern"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

png(output_png, width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

message("Finished: ", output_png)