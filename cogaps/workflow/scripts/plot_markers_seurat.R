suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})

# Inputs / Outputs
contrib_file <- snakemake@input[["contrib"]]
seurat_file <- snakemake@input[["tumor_rna_seuratobj"]]
markers_file <- snakemake@input[["markers_excel"]]
output_png <- snakemake@output[["png"]]

# Read data
seurat <- readRDS(seurat_file)
contrib <- read_tsv(contrib_file, show_col_types = FALSE)

# Keep only common cells
common_cells <- intersect(contrib$barcode, colnames(seurat))
contrib <- contrib %>% filter(barcode %in% common_cells)
seurat <- subset(seurat, cells = common_cells)

# Add dominant pattern / cluster to Seurat object
seurat$dominant_pattern <- contrib$dominant_pattern[match(colnames(seurat), contrib$barcode)]
Idents(seurat) <- "dominant_pattern"

# Normalize & scale
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)

# --- Read all sheets and get top genes per pattern ---
sheets <- excel_sheets(markers_file)
top_features <- lapply(sheets, function(sheet) {
  df <- read_xlsx(markers_file, sheet = sheet)
  df %>%
    slice(1:5) %>%  # top 5 gènes par pattern
    pull(gene)
}) |> unlist() |> unique() |> na.omit()

# Filter out uninformative genes
top_features <- top_features[!grepl("^Gm|Rik", top_features)]

# --- Color setup ---
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

labels_x <- sapply(top_features, function(gene) {
  if (gene %in% c("Krt14", "TdTomato")) {
    paste0("<span style='color:red; font-weight:bold'>", gene, "</span>")
  } else {
    gene
  }
})

# --- DotPlot ---
p <- DotPlot(seurat,
             features = top_features,
             group.by = "dominant_pattern") +
  RotatedAxis() +
  scale_y_discrete(labels = colored_labels) +
  scale_x_discrete(labels = labels_x) +
  theme(
    axis.text.y = element_markdown(size = 12, face = "bold"),
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10, face = "bold")
  )

# Save png
png(output_png, width = 10, height = 8)
print(p)
dev.off()