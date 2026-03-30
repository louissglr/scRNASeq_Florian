library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(grid)

# Inputs/outputs Snakemake
input_csv_enrichment  <- snakemake@input[["csv_enrich"]]
input_csv_overrep     <- snakemake@input[["csv_overr"]]

output_pdf_enrich     <- snakemake@output[["pdf_heatmap_enrich"]]
output_pdf_overrep    <- snakemake@output[["pdf_heatmap_overr"]]

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- as.integer(snakemake@wildcards[["npatterns"]])

generate_heatmap_continuous <- function(input_csv, output_pdf, title_suffix, max_val_cutoff = 10) {
  
  # Lecture et nettoyage des noms de hallmarks
  hallmarks_df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
    mutate(
      gene.set = str_replace(gene.set, "^HALLMARK[_ ]?", ""),
      neg_log10_padj = -log10(padj)   # recalcul correct
    )
  
  # Pivot table
  heatmap_df <- hallmarks_df %>%
    select(gene.set, Pattern, neg_log10_padj) %>%
    mutate(Pattern = factor(Pattern, levels = sort(unique(Pattern)))) %>%
    pivot_wider(
      names_from  = Pattern,
      values_from = neg_log10_padj,
      values_fill = 0
    )
  
  # Conversion en matrice
  heatmap_mat <- heatmap_df %>%
    column_to_rownames("gene.set") %>%
    as.matrix()
  
  colnames(heatmap_mat) <- str_replace(colnames(heatmap_mat), "^Pattern_", "")
  
  # Appliquer cutoff max
  heatmap_mat[heatmap_mat >= max_val_cutoff] <- max_val_cutoff
  
  # Définition du gradient continu
  col_fun <- colorRamp2(
    c(0, max_val_cutoff/2, max_val_cutoff),
    c("white", "#ffcccc", "darkred")
  )
  
  # Création de la heatmap
  ht <- Heatmap(
    heatmap_mat,
    name = "-log10(padj)",
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
      title = "-log10(padj)",
      border = "black",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, col = "black")
    )
  )
  
  # Sauvegarde en PDF
  pdf(output_pdf, width = 8, height = 10)
  draw(ht)
  dev.off()
}

# Génération heatmap enrichment
generate_heatmap_continuous(
  input_csv_enrichment,
  output_pdf_enrich,
  "enrichment",
  max_val_cutoff = 10
)

# Génération heatmap overrepresentation
generate_heatmap_continuous(
  input_csv_overrep,
  output_pdf_overrep,
  "overrepresentation",
  max_val_cutoff = 10
)


