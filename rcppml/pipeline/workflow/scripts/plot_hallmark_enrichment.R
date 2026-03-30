suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
})

input_csv <- snakemake@input[["csv"]]        # CSV produit par hallmark_enrichment.R
output_pdf <- snakemake@output[["pdf"]]     # PDF à générer

padj_cutoff <- 0.05
threshold <- -10 * log10(padj_cutoff)

# ============================
# Lecture du CSV d'enrichissement
# ============================
hallmarks_df <- read.csv(input_csv, stringsAsFactors = FALSE)

if (nrow(hallmarks_df) == 0) {
  message("Aucun enrichissement à afficher, PDF vide créé.")
  pdf(output_pdf, width = 11, height = 8)
  plot.new()
  text(0.5, 0.5, "Aucun enrichissement détecté", cex = 1.5)
  dev.off()
  quit(save = "no")
}

# Calcul du -log10(p.adjust)
hallmarks_df <- hallmarks_df %>%
  mutate(neg.log.padj = -log10(p.adjust))

# ============================
# Boucle par programme
# ============================
program_levels <- unique(hallmarks_df$program)

plot_list <- lapply(program_levels, function(prog) {
  df <- hallmarks_df %>% filter(program == prog)
  df <- df[order(-df$neg.log.padj), ]
  df$ID <- factor(df$ID, levels = rev(df$ID))
  
  ggplot(df, aes(x = neg.log.padj, y = ID, fill = neg.log.padj)) +
    geom_col() +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_fill_viridis_c(option = "C") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    ) +
    labs(
      title = paste0("Hallmark Enrichment - Program ", prog),
      x = expression(-10 * log[10](padj)),
      y = "Hallmark gene set"
    )
})

# ============================
# Export PDF
# ============================
pdf(output_pdf, width = 11, height = 8)
for (p in plot_list) {
  print(p)
}
dev.off()

message("PDF généré : ", output_pdf)
