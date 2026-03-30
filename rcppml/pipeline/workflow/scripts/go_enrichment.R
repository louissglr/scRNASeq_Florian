suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)   # Adapté à la souris
  library(enrichplot)
  library(ggplot2)
})

# Inputs/outputs via Snakemake
markers_file <- snakemake@input[["tsv"]]
pdf_file     <- snakemake@output[["pdf"]]

# Lecture des top gènes
markers <- read.csv(markers_file, stringsAsFactors = FALSE)

# Si les gènes sont en colonnes pour chaque pattern (comme dans top_genes)
gene_list <- lapply(seq_len(ncol(markers)), function(i) {
  genes <- as.character(na.omit(markers[[i]]))
  genes
})
names(gene_list) <- colnames(markers)

# GO enrichment pour chaque pattern
ego_list <- list()
for(pattern_name in names(gene_list)) {
  genes <- gene_list[[pattern_name]]
  
  if(length(genes) > 0){
    ego <- enrichGO(
      gene         = genes,
      OrgDb        = org.Mm.eg.db,   # souris
      keyType      = "SYMBOL",
      ont          = "BP",            # Biological Process
      pAdjustMethod= "BH",
      qvalueCutoff = 0.05,
      readable     = TRUE
    )
    ego_list[[pattern_name]] <- ego
  }
}

# Fusionner tous les résultats
merged <- merge_result(ego_list)

# Générer le PDF
pdf(pdf_file, width = 12, height = 10)
dotplot(
  merged,
  showCategory = 10,
  font.size = 8,
  label_format = 40,
  title = "GO Enrichment (BP) - Top 10 categories shown"
)
dev.off()
