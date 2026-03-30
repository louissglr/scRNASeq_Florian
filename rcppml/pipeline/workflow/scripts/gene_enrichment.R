suppressPackageStartupMessages({
  #library(clusterProfiler)
  library(msigdbr)
  library(openxlsx)
  library(dplyr)
})

tsv_file <- snakemake@input[["tsv"]]
output_csv_enrich <- snakemake@output[["csv_enrich"]]

top_genes_df <- read.xlsx(tsv_file)
gene_names <- rownames(top_genes_df)
if (is.null(gene_names)) {
  gene_names <- top_genes_df[[1]]
}

detect_species <- function(gene_names) {
  human_prop <- mean(gene_names == toupper(gene_names), na.rm = TRUE)
  mouse_prop <- mean(grepl("^[A-Z][a-z0-9]+$", gene_names), na.rm = TRUE)
  
  if (human_prop > 0.8) {
    return("Homo sapiens")
  } else if (mouse_prop > 0.8) {
    return("Mus musculus")
  } else {
    return("Homo sapiens")
  }
}

species <- detect_species(gene_names)
message("Espèce détectée : ", species)

msigdbr_df <- msigdbr(species = species, category = "H")
hallmark_list <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)
hallmark_list <- lapply(hallmark_list, unique)

# Enrichissement pour chaque pattern (programme NMF)
enrich_list <- lapply(seq_len(ncol(top_genes_df)), function(i) {
  program_genes <- as.character(na.omit(top_genes_df[[i]]))
  lapply(names(hallmark_list), function(gs_name) {
    geneset <- hallmark_list[[gs_name]]
    overlap <- intersect(geneset, program_genes)
    p_value <- phyper(length(overlap) - 1,
                      length(geneset),
                      length(unique(msigdbr_df$gene_symbol)) - length(geneset),
                      length(program_genes),
                      lower.tail = FALSE)
    data.frame(
      gene.set = gs_name,
      overlap = paste(overlap, collapse = ";"),
      n_overlap = length(overlap),
      n_geneset = length(geneset),
      n_input = length(program_genes),
      p_value = p_value,
      pattern = paste0("Pattern_", i),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Trier par pattern puis p-value
enrich_list <- enrich_list %>%
  mutate(pattern_num = as.numeric(gsub("Pattern_", "", pattern))) %>%
  arrange(pattern_num, p_value) %>%
  select(-pattern_num)

write.csv(enrich_list, output_csv_enrich, row.names = FALSE)

