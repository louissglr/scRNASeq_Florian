suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(uwot)
  library(grid) 
})

# -------------------
# Inputs/outputs Snakemake
# -------------------
contrib_file <- snakemake@input[["tsv"]]
seurat_file  <- snakemake@input[["seurat_rds"]]  # pour récupérer barcodes si nécessaire
output_pdf   <- snakemake@output[["pdf"]]
sample_name  <- snakemake@wildcards[["sample"]]

# Créer dossier si nécessaire
if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

# -------------------
# Charger les contributions
# -------------------
cell_contrib <- read.delim(contrib_file)

# Extraire uniquement les colonnes nmf*
pattern_cols <- grep("^programs", colnames(cell_contrib), value = TRUE)
mat <- as.matrix(cell_contrib[, pattern_cols])
pattern_names <- colnames(mat)
n_patterns <- ncol(mat)

# -------------------
# Fonction de plot
# -------------------
plot_patterns <- function(df_coords, coord_names = c("UMAP1", "UMAP2"), title_prefix = "UMAP Patterns") {
  
  df_coords$cell <- rownames(df_coords)
  df_coords <- cbind(df_coords, mat)
  
  df_long <- df_coords %>%
    pivot_longer(cols = all_of(pattern_names),
                 names_to = "Pattern",
                 values_to = "weight") %>%
    arrange(weight) 
  
  x_min <- min(df_long[[coord_names[1]]], na.rm = TRUE)
  x_max <- max(df_long[[coord_names[1]]], na.rm = TRUE)
  y_min <- min(df_long[[coord_names[2]]], na.rm = TRUE)
  y_max <- max(df_long[[coord_names[2]]], na.rm = TRUE)
  x_offset <- 0.05 * (x_max - x_min)
  y_offset <- 0.05 * (y_max - y_min)
  arrow_length <- 1.5
  
  p <- ggplot(df_long, aes_string(x = coord_names[1], y = coord_names[2], color = "weight")) +
    geom_point(size = 0.8, alpha = 1) +
    scale_color_viridis_c(option = "plasma") +
    facet_wrap(~Pattern, ncol = 3) +
    theme_void() +
    theme(strip.text = element_text(size = 10),
          legend.position = "right") +
    guides(color = guide_colorbar(barwidth = 0.8, barheight = 5)) +
    geom_segment(x = x_min - x_offset, y = y_min - y_offset,
                 xend = x_min - x_offset + arrow_length, yend = y_min - y_offset,
                 colour = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
    geom_segment(x = x_min - x_offset, y = y_min - y_offset,
                 xend = x_min - x_offset, yend = y_min - y_offset + arrow_length,
                 colour = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", x = x_min - x_offset + arrow_length/2, y = y_min - y_offset - 0.5,
             label = coord_names[1], color = "black", size = 3) +
    annotate("text", x = x_min - x_offset - 0.5, y = y_min - y_offset + arrow_length/2,
             label = coord_names[2], color = "black", size = 3, angle = 90) +
    labs(title = paste0(title_prefix, " — ", sample_name, " (", n_patterns, " patterns)"),
         color = "weight")
  
  return(p)
}

# -------------------
# UMAP
# -------------------
set.seed(123)
umap_res <- uwot::umap(mat, n_neighbors = 30, min_dist = 0.3, metric = "cosine")
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
p_umap <- plot_patterns(umap_df, c("UMAP1", "UMAP2"), "UMAP Patterns")

# -------------------
# PCA
# -------------------
pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
p_pca <- plot_patterns(pca_df, c("PC1", "PC2"), "PCA Patterns")

# -------------------
# Sauvegarde PDF
# -------------------
pdf(output_pdf, width = 14, height = 10)
print(p_umap)
print(p_pca)
dev.off()

message("End: ", output_pdf)

