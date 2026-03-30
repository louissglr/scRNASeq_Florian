suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

contributions_file <- snakemake@input[["contributions"]]
meta_file <- snakemake@input[["meta"]]
output_pdf <- snakemake@output[["pdf"]]
sample_name <- snakemake@wildcards[["sample"]]
n_patterns <- as.integer(snakemake@wildcards[["npatterns"]])

contrib <- read.table(
  contributions_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

meta <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

contrib_meta <- contrib |>
  inner_join(meta, by = c("barcode" = "cell_id"))


contrib_meta <- contrib_meta |>
  mutate(
    group = paste(orig_tumor, tnt, Sting_status, sep = "_")
  )


contrib_long <- contrib_meta |>
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "Pattern",
    values_to = "Weight"
  ) |>
  mutate(
    Pattern = factor(
      Pattern,
      levels = paste0("Pattern_", seq_len(n_patterns))
    ),
    group = factor(group)
  )


p <- ggplot(contrib_long, aes(x = Pattern, y = Weight, fill = Pattern)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Pattern",
    y = "Weight",
    title = paste(
      "Distribution des poids par pattern\nSample:",
      sample_name,
      "- n_patterns:",
      n_patterns
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


p_stat <- ggboxplot(
  contrib_long,
  x = "Pattern",
  y = "Weight",
  color = "Pattern",
  palette = "npg",
  outlier.size = 0.5,
  legend = "none"
) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format"
  ) +
  geom_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "{p.adj.signif}",
    hide.ns = TRUE
  )


p_group_pattern <- ggplot(contrib_long, aes(x = group, y = Weight, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ Pattern, scales = "free_y") +
  labs(
    title = paste(
      "Distribution des poids par groupe dans chaque pattern\nSample:",
      sample_name
    ),
    x = "Groupe",
    y = "Weight"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


p_group_pattern_stat <- p_group_pattern +
  stat_compare_means(
    aes(group = group),
    method = "kruskal.test",
    label = "p.format"
  )


pdf(output_pdf, width = 14, height = 10)

print(p)
print(p_stat)
print(p_group_pattern)
print(p_group_pattern_stat)

dev.off()

message("End: ", output_pdf)