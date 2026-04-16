library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggtext)

contrib_file <- snakemake@input[["contrib"]]
seurat_file <- snakemake@input[["tumor_rna_seuratobj"]]
markers_file <- snakemake@input[["markers_excel"]]
output_png  <- snakemake@output[["png"]]

contrib <- read.table(contrib_file, header = TRUE)
seurat <- readRDS(seurat_file)
markers <- read_xlsx(markers_file)

rownames(contrib) <- contrib$barcode
seurat$dominant_pattern <- contrib[colnames(seurat), "dominant_pattern"]

features <- markers %>%
  slice(1:5) %>%
  select(starts_with("Pattern_")) %>%
  unlist() %>%
  unique() %>%
  na.omit() %>%
  .[!grepl("^Gm|Rik", .)]

seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)

#https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

pattern_levels <- sort(unique(seurat$dominant_pattern))
color_list <- ggplotColours(n = length(pattern_levels))
pattern_colors <- setNames(color_list, pattern_levels)

colored_labels <- sapply(pattern_levels, function(x) {
  color <- pattern_colors[x]
  paste0("<span style='color:", color, "'>", x, "</span>")
})
names(colored_labels) <- pattern_levels

labels_x <- sapply(features, function(gene) {
  if (gene %in% c("Krt14", "TdTomato")) {
    paste0("<span style='color:red; font-weight:bold'>", gene, "</span>")
  } else {
    gene
  }
})

# --- DotPlot ---
p <- DotPlot(seurat,
        features = features,
        group.by = "dominant_pattern") +
  RotatedAxis() +
  scale_y_discrete(labels = colored_labels) +
  scale_x_discrete(labels = labels_x) +
  theme(
    axis.text.y = element_markdown(size = 12, face = "bold"),
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10, face = "bold")
  )

png(output_png, width = 10, height = 8)
print(p)
dev.off()