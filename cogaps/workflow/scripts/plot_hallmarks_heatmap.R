library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(grid)  # pour gpar

# Inputs/outputs Snakemake
input_csv_enrichment  <- snakemake@input[["csv_enrich"]]
input_csv_overrep     <- snakemake@input[["csv_overr"]]

output_png_enrich     <- snakemake@output[["png_heatmap_enrich"]]
output_png_overrep    <- snakemake@output[["png_heatmap_overr"]]

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- as.integer(snakemake@wildcards[["npatterns"]])

# Paramètres
max_color_val <- 50
sig_threshold <- -10*log10(0.05)

generate_heatmap <- function(input_csv, output_png, title_suffix) {
  
  # Lecture et nettoyage des noms de hallmarks
  hallmarks_df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
    mutate(gene.set = str_replace(gene.set, "^HALLMARK[_ ]?", ""))
  
  # Pivot et filtrage des hallmarks significatifs
  heatmap_df <- hallmarks_df %>%
    select(gene.set, Pattern, neg.log.padj) %>%
    mutate(Pattern = factor(Pattern, levels = sort(unique(Pattern)))) %>%
    pivot_wider(
      names_from  = Pattern,
      values_from = neg.log.padj,
      values_fill = 0
    ) %>%
    (\(df) df[rowSums(df[, -1] >= sig_threshold) > 0, ])()
  
  # Conversion en matrice
  heatmap_mat <- heatmap_df %>%
    column_to_rownames("gene.set") %>%
    as.matrix()
  
  colnames(heatmap_mat) <- str_replace(colnames(heatmap_mat), "^Pattern_", "")
  
  # Cap valeurs max et mise à zéro si non significatif
  heatmap_mat[heatmap_mat > max_color_val] <- max_color_val
  heatmap_mat[heatmap_mat < sig_threshold] <- 0
  
  # Gradient couleur
  col_fun <- colorRamp2(
    c(0, sig_threshold, sig_threshold + 0.001, max_color_val),
    c("white", "white", "#ffe5e5", "darkred")
  )
  
  # Création de la heatmap
  ht <- Heatmap(
    heatmap_mat,
    name = "-10*log10(padj)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 0,
    rect_gp = gpar(col = "black", lwd = 0.5),
    column_title = paste0("Hallmark ", title_suffix, " (", sample_name, ")"),
    
    heatmap_legend_param = list(
      title = "-10*log10(padj)",
      at = c(0, sig_threshold, seq(10, max_color_val, by = 10)),
      labels = c(
        "NS (padj > 0.05)",
        "S (padj ≤ 0.05)",
        as.character(seq(10, max_color_val - 10, by = 10)),
        paste0("≥ ", max_color_val)
      ),
      border = "black",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Sauvegarde au format png
  png(output_png, width = 10, height = 8, units = "in", res = 300)
  draw(ht)
  dev.off()
}

# Génération heatmap enrichment
generate_heatmap(
  input_csv_enrichment,
  output_png_enrich,
  "enrichment"
)

# Génération heatmap overrepresentation
generate_heatmap(
  input_csv_overrep,
  output_png_overrep,
  "overrepresentation"
)
