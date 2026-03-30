suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  library(grid)
})

# Inputs/outputs Snakemake
input_csv_enrichment  <- snakemake@input[["csv_enrich"]]
input_csv_overrep     <- snakemake@input[["csv_overr"]]

output_pdf_enrich     <- snakemake@output[["pdf_heatmap_enrich"]]
output_pdf_overrep    <- snakemake@output[["pdf_heatmap_overr"]]

sample_name <- snakemake@wildcards[["sample"]]

# Paramètres
max_color_val <- 50
sig_threshold <- -10 * log10(0.05)

generate_heatmap <- function(input_csv, output_pdf, title_suffix) {
  
  enrich_df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
    mutate(
      neg.log.padj = -10 * log10(p.adjust),
      gene.set = str_replace(ID, "^HALLMARK[_ ]?", "")
    )
  
  # Pivot et filtrage
  heatmap_df <- enrich_df %>%
    select(gene.set, program, neg.log.padj) %>%
    pivot_wider(
      names_from  = program,
      values_from = neg.log.padj,
      values_fill = 0
    ) %>%
    (\(df) df[rowSums(df[, -1] >= sig_threshold) > 0, ])()
  
  # Matrice
  heatmap_mat <- heatmap_df %>%
    column_to_rownames("gene.set") %>%
    as.matrix()
  
  colnames(heatmap_mat) <- str_replace(colnames(heatmap_mat), "^Pattern_", "")
  
  # Nettoyage
  heatmap_mat[heatmap_mat > max_color_val] <- max_color_val
  heatmap_mat[heatmap_mat < sig_threshold] <- 0
  
  # Couleurs
  col_fun <- colorRamp2(
    c(0, sig_threshold, sig_threshold + 0.001, max_color_val),
    c("white", "white", "#ffe5e5", "darkred")
  )
  
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
    column_title = paste0("Hallmark enrichment (", sample_name, ")"),
    
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
      labels_gp = gpar(fontsize = 8)
    )
  )
  
  pdf(output_pdf, width = 8, height = 10)
  draw(ht)
  dev.off()
}

# Appels
generate_heatmap(
  input_csv_enrichment,
  output_pdf_enrich,
  "enrichment"
)

generate_heatmap(
  input_csv_overrep,
  output_pdf_overrep,
  "overrepresentation"
)