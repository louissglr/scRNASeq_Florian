suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
})


umap_file    <- snakemake@input[["umap"]]
metadata_file <- snakemake@input[["metadata"]]
output_pdf   <- snakemake@output[["pdf"]]
sample_name  <- snakemake@wildcards[["sample"]]
dataset      <- snakemake@params[["dataset"]]


umap_df   <- read_tsv(umap_file)
metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

colnames(umap_df) <- tolower(colnames(umap_df))
colnames(metadata) <- tolower(colnames(metadata))

metadata <- metadata |>
  mutate(group = paste(orig_tumor, tnt, sting_status, sep = "_"))

umap_cols <- colnames(umap_df)[grepl("umap|dim|coord|x$|y$", colnames(umap_df))]
barcode_col <- intersect(colnames(umap_df), c("barcode","barcodes","cell","cells","cell_id","cellid"))[1]

umap_df <- umap_df |>
  rename(
    barcode = all_of(barcode_col),
    UMAP_1  = all_of(umap_cols[1]),
    UMAP_2  = all_of(umap_cols[2])
  )


plot_df <- umap_df |>
  inner_join(metadata |>
               select(cell_id, group),
             by = c("barcode" = "cell_id"))

plot_df$group <- factor(plot_df$group)


x_min <- min(plot_df$UMAP_1)
x_max <- max(plot_df$UMAP_1)
y_min <- min(plot_df$UMAP_2)
y_max <- max(plot_df$UMAP_2)
x_offset <- 0.05 * (x_max - x_min)
y_offset <- 0.05 * (y_max - y_min)
arrow_length <- 1.5

p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
  geom_point(size = 0.8, alpha = 0.9) +
  theme_void() +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  geom_segment(x = x_min - x_offset, y = y_min - y_offset,
               xend = x_min - x_offset + arrow_length, yend = y_min - y_offset,
               colour = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(x = x_min - x_offset, y = y_min - y_offset,
               xend = x_min - x_offset, yend = y_min - y_offset + arrow_length,
               colour = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = x_min - x_offset + arrow_length/2, y = y_min - y_offset - 0.5,
           label = "UMAP1", size = 3) +
  annotate("text", x = x_min - x_offset - 0.5, y = y_min - y_offset + arrow_length/2,
           label = "UMAP2", size = 3, angle = 90) +
  labs(title = paste0("UMAP - Metadata groups - ", sample_name, " (", dataset, ")"),
       color = "Group") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")

pdf(output_pdf, width = 10, height = 8)
print(p)
dev.off()