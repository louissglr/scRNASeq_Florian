library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)


contrib_file <- snakemake@input[["contrib"]]
seurat_file <- snakemake@input[["tumor_rna_seuratobj"]]
signature_file <- snakemake@input[["signature"]]
coo_file <- snakemake@input[["coo"]]
png_file <- snakemake@output[["png"]]
sample_name <- snakemake@wildcards[["sample"]]

contrib <- read.table(contrib_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
seurat <- readRDS(seurat_file)
signature <- read_excel(signature_file)
coo <- read.table(coo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(contrib) <- contrib$barcode
seurat$dominant_pattern <- contrib[colnames(seurat), "dominant_pattern"]

sig_p3 <- signature |>
  dplyr::pull(Pattern_3) |>
  na.omit() |>
  as.character()

sig_p3 <- intersect(sig_p3, rownames(seurat))

DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat)

seurat <- AddModuleScore(
  object = seurat,
  features = list(sig_p3),
  assay = "RNA",
  slot = "data",
  name = "SigP3"
)

score_name <- "SigP31"

df_plot <- data.frame(
  cell_id = colnames(seurat),
  score = seurat[[score_name]][,1],
  dominant_pattern = seurat$dominant_pattern 
)

df_plot <- df_plot |>
  left_join(coo, by = "cell_id") |>
  mutate(dominant_pattern = str_remove(dominant_pattern, "Pattern_"))


p_umap <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = score)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "viridis") +
  theme_classic() +
  labs(color = "SigP3 score")

p_violin <- ggplot(df_plot, aes(x = as.factor(dominant_pattern), y = score, fill = as.factor(dominant_pattern))) +
  geom_violin(trim = FALSE, alpha = 1) +      
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1) +  
  theme_classic() +
  labs(
    x = "Dominant Pattern",
    y = "SigP3 score",
    title = "Distribution des scores SigP3 par pattern",
    fill = "Dominant Pattern"
  ) +
  theme(legend.position = "none")

  png(png_file, width = 10, height = 8, units = "in", res = 300)
p_umap + p_violin
dev.off()
