library(ggplot2)

input_csv_enrichment  <- snakemake@input[["csv_enrich"]]
output_png_enrichment <- snakemake@output[["png_enrich"]]
input_csv_overrep  <- snakemake@input[["csv_overr"]]
output_png_overrep <- snakemake@output[["png_overr"]]



padj_cutoff <- 0.05
threshold <- -10 * log10(padj_cutoff)

# ============================
# Enrichment
# ============================

hallmarks_df <- read.csv(input_csv_enrichment, stringsAsFactors = FALSE)

pattern_levels <- unique(hallmarks_df$Pattern)
pattern_df_list <- lapply(pattern_levels, function(pat) {
  subset(hallmarks_df, Pattern == pat)
})

pattern_plot_list <- lapply(pattern_df_list, function(df) {
  df <- df[order(-df$neg.log.padj), ]
  df$gene.set <- factor(df$gene.set, levels = rev(df$gene.set))
  
  ggplot(df, aes(x = neg.log.padj, y = gene.set, fill = NES)) +
    geom_col() +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_fill_viridis_c(option = "C") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Hallmark Enrichment - ", unique(df$Pattern)),
      x = expression(-10 * log[10](padj)),
      y = "Hallmark gene set"
    )
})

png(output_png_enrichment, width = 11, height = 8)
for (p in pattern_plot_list) {
  print(p)
}
dev.off()
message("png généré : ", output_png_enrichment)

# ============================
# Overrepresentation
# ============================

hallmarks_df <- read.csv(input_csv_overrep, stringsAsFactors = FALSE)
names(hallmarks_df)[names(hallmarks_df) == "k/K"] <- "k_K"
pattern_levels <- unique(hallmarks_df$Pattern)
pattern_df_list <- lapply(pattern_levels, function(pat) {
  subset(hallmarks_df, Pattern == pat)
})

pattern_plot_list <- lapply(pattern_df_list, function(df) {
  df <- df[order(-df$neg.log.padj), ]
  df$gene.set <- factor(df$gene.set, levels = rev(df$gene.set))
  
  ggplot(df, aes(x = neg.log.padj, y = gene.set, fill = neg.log.padj)) +
    geom_col() +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_fill_viridis_c(option = "C") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Hallmark Overrepresentation -  ", unique(df$Pattern)),
      x = expression(-10 * log[10](padj)),
      y = "Hallmark gene set"
    )
})

png(output_png, width = 10, height = 8, units = "in", res = 300)
for (p in pattern_plot_list) {
  print(p)
}
dev.off()
message("png généré : ", output_png_overrep)

