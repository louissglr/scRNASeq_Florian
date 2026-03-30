suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -------------------
# Inputs/outputs Snakemake
# -------------------
contrib_path <- snakemake@input[["contrib_tsv"]]  # nouveau fichier en entrée
output_path <- snakemake@output[["tsv"]]
sample_name <- snakemake@wildcards[["sample"]]

top_n <- 50  # nombre de gènes top par programme

# -------------------
# Lecture du fichier de contributions
# -------------------
gene_contrib <- read_tsv(contrib_path)


# -------------------
# Extraction des top gènes
# -------------------
top_genes_list <- lapply(
  gene_contrib %>% select(-gene),
  function(x) {
    gene_contrib$gene[order(x, decreasing = TRUE)][1:top_n]
  }
)

# Convertir en data.frame
top_genes_df <- as.data.frame(top_genes_list)

# -------------------
# Export
# -------------------
write.csv(top_genes_df, output_path, row.names = FALSE)