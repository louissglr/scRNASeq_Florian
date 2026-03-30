suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(grid)
})

# --------------------------------------------------
# INPUTS
# --------------------------------------------------

contrib_file <- snakemake@input[["contributions"]]
umap_file    <- snakemake@input[["umap"]]
output_pdf   <- snakemake@output[["pdf"]]
sample_name  <- snakemake@wildcards[["sample"]]
dataset      <- snakemake@params[["dataset"]]

if (!dir.exists(dirname(output_pdf)))
  dir.create(dirname(output_pdf), recursive = TRUE)

# --------------------------------------------------
# 1. Read Contributions
# --------------------------------------------------

contrib_df <- read_tsv(contrib_file, guess_max = 1e6)

colnames(contrib_df) <- tolower(colnames(contrib_df))
colnames(contrib_df) <- gsub("[^a-z0-9_]", "_", colnames(contrib_df))

barcode_candidates <- c("barcode","barcodes","cell","cells","cell_id","cellid")
barcode_col <- intersect(colnames(contrib_df), barcode_candidates)[1]

if (is.na(barcode_col))
  stop("No barcode column detected")

exclude_cols <- c("sample_name","dominant_program",barcode_col)
pattern_cols <- setdiff(colnames(contrib_df), exclude_cols)

if (length(pattern_cols) == 0)
  stop("No pattern columns detected")

# --------------------------------------------------
# 2. Read UMAP
# --------------------------------------------------

clean_umap <- function(file){

  df <- read_tsv(file, guess_max = 1e6)

  colnames(df) <- tolower(colnames(df))
  colnames(df) <- gsub("[^a-z0-9_]", "_", colnames(df))

  barcode_candidates <- c("barcode","barcodes","cell","cells","cell_id","cellid")
  barcode_col <- intersect(colnames(df), barcode_candidates)[1]

  if (is.na(barcode_col))
    stop("No barcode-like column in UMAP")

  umap_candidates <- colnames(df)[grepl("umap|dim|coord|x$|y$", colnames(df))]

  if (length(umap_candidates) < 2)
    stop("No UMAP coordinates detected")

  umap_x <- umap_candidates[1]
  umap_y <- umap_candidates[2]

  df |> rename(
    barcode = all_of(barcode_col),
    UMAP1   = all_of(umap_x),
    UMAP2   = all_of(umap_y)
  ) |> select(barcode, UMAP1, UMAP2, everything())
}

umap_df <- clean_umap(umap_file)

# --------------------------------------------------
# 3. Merge + Long Format (Weights)
# --------------------------------------------------

plot_df <- umap_df |>
  left_join(
    contrib_df |> select(all_of(barcode_col), all_of(pattern_cols)),
    by = c("barcode" = barcode_col)
  ) |>
  filter(if_any(all_of(pattern_cols), ~ !is.na(.)))

plot_long <- plot_df |>
  pivot_longer(
    cols = all_of(pattern_cols),
    names_to = "Pattern",
    values_to = "weight"
  )

pattern_names <- unique(plot_long$Pattern)
n_patterns <- length(pattern_names)

# --------------------------------------------------
# 4. FIXED PALETTE (Same as your dominant script)
# --------------------------------------------------

pattern_colors <- colorRampPalette(
  brewer.pal(min(n_patterns,12), "Set3")
)(n_patterns)

names(pattern_colors) <- pattern_names

# --------------------------------------------------
# 5. Function Theme With Axes Arrows
# --------------------------------------------------

add_axes_arrows <- function(p, df, xcol, ycol){

  x_min <- min(df[[xcol]], na.rm = TRUE)
  x_max <- max(df[[xcol]], na.rm = TRUE)
  y_min <- min(df[[ycol]], na.rm = TRUE)
  y_max <- max(df[[ycol]], na.rm = TRUE)

  x_offset <- 0.05 * (x_max - x_min)
  y_offset <- 0.05 * (y_max - y_min)
  arrow_length <- 1.5

  p +
    geom_segment(
      x = x_min - x_offset,
      y = y_min - y_offset,
      xend = x_min - x_offset + arrow_length,
      yend = y_min - y_offset,
      colour = "black",
      size = 0.5,
      arrow = arrow(length = unit(0.2,"cm"))
    ) +
    geom_segment(
      x = x_min - x_offset,
      y = y_min - y_offset,
      xend = x_min - x_offset,
      yend = y_min - y_offset + arrow_length,
      colour = "black",
      size = 0.5,
      arrow = arrow(length = unit(0.2,"cm"))
    ) +
    annotate("text",
             x = x_min - x_offset + arrow_length/2,
             y = y_min - y_offset - 0.5,
             label = xcol,
             size = 3) +
    annotate("text",
             x = x_min - x_offset - 0.5,
             y = y_min - y_offset + arrow_length/2,
             label = ycol,
             size = 3,
             angle = 90)
}

# --------------------------------------------------
# 6. PLOT WEIGHTS WITH FACET
# --------------------------------------------------

p_weights <- ggplot(plot_long,
                    aes(x = UMAP1,
                        y = UMAP2,
                        color = weight)) +
  geom_point(size = 0.8, alpha = 1) +
  facet_wrap(~ Pattern, ncol = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "right"
  ) +
  guides(color = guide_colorbar(barwidth = 0.8, barheight = 5)) +
  labs(
    title = paste0("UMAP Weights — ", sample_name,
                   " (", n_patterns, " patterns)"),
    color = "Weight"
  )

# Ajouter les flèches axes (calcul global sur données UMAP)
p_weights <- add_axes_arrows(p_weights,
                             plot_long,
                             "UMAP1",
                             "UMAP2")

# --------------------------------------------------
# 7. EXPORT
# --------------------------------------------------

pdf(output_pdf, width = 14, height = 10)
print(p_weights)
dev.off()

message("Finished: ", output_pdf)
