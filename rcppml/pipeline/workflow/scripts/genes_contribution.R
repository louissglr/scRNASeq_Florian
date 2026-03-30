suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(readr)
  library(grid)
  library(openxlsx)
})

# -------------------
# Inputs/outputs Snakemake
# -------------------
rds_path <- snakemake@input[["rcppml_rds"]]
seurat_rds <- snakemake@input[["seurat_rds"]]
output_path <- snakemake@output[["tsv"]]  
sample_name <- snakemake@wildcards[["sample"]]

# -------------------
# Lecture des objets
# -------------------
model <- readRDS(rds_path)
seurat_obj <- readRDS(seurat_rds)

gene_names <- rownames(seurat_obj)

# -------------------
# Préparer les contributions
# -------------------
gene_contrib <- as.data.frame(model@w)
rownames(gene_contrib) <- gene_names
colnames(gene_contrib) <- paste0("programs_", 1:ncol(gene_contrib))

# Ajouter une colonne "gene" pour le nom des gènes
gene_contrib_export <- gene_contrib %>%
  mutate(gene = rownames(gene_contrib)) %>%
  select(gene, everything())  

# -------------------
# Export TXT tab-delimited
# -------------------
write.table(
  gene_contrib_export,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
message("TXT exporté avec toutes les contributions : ", output_path)