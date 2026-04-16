library(CoGAPS)
library(msigdbr)

set.seed(123)

detect_species <- function(gene_names) {
  gene_names <- unique(gene_names)
  human_prop <- mean(gene_names == toupper(gene_names), na.rm = TRUE)
  mouse_prop <- mean(grepl("^[A-Z][a-z0-9]+$", gene_names), na.rm = TRUE)
  
  if (human_prop > 0.8) {
    return("Homo sapiens")
  } else if (mouse_prop > 0.8) {
    return("Mus musculus")
  } else {
    message("Espèce ni humain ni souris")
    return("Homo sapiens")
  }
}

cogaps_file <- snakemake@input[["cogaps_rds"]]
output_csv_enrichment  <- snakemake@output[["csv_enrich"]]
output_csv_overrepresentation <- snakemake@output[["csv_overr"]]

cogaps <- readRDS(cogaps_file)

fl <- cogaps@featureLoadings

dup_tab <- table(rownames(fl))
dup_genes <- names(dup_tab[dup_tab > 1])

if (length(dup_genes) > 0) {
  keep <- !(rownames(fl) %in% dup_genes)
  cogaps@featureLoadings <- fl[keep, , drop = FALSE]
}

gene_names <- rownames(cogaps@featureLoadings)
species_detected <- detect_species(gene_names)


msigdbr_df <- msigdbr(species = species_detected, category = "C8")

hallmark_ls <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)
hallmark_ls <- lapply(hallmark_ls, unique)


# --------------------------
# Bloc 1 : enrichment
# --------------------------
c2_ora_enrichment <- getPatternGeneSet(
  cogaps,
  gene.sets = hallmark_ls,
  method = "enrichment"
)

df_list_enrichment <- lapply(seq_along(c2_ora_enrichment), function(i) {
  df <- c2_ora_enrichment[[i]]
  df$Pattern <- paste0("Pattern_", i)
  df
})

c2_df_enrichment <- do.call(rbind, df_list_enrichment)
c2_df_enrichment$gene.set <- as.character(c2_df_enrichment$gene.set)

c2_df_enrichment <- c2_df_enrichment[, -7]
c2_df_enrichment <- as.data.frame(c2_df_enrichment)

write.csv(
  c2_df_enrichment,
  output_csv_enrichment,
  row.names = FALSE
)

message("Fichier produit : ", output_csv_enrichment)

# --------------------------
# Bloc 2 : overrepresentation
# --------------------------
c2_ora_overrep <- getPatternGeneSet(
  cogaps,
  gene.sets = hallmark_ls,
  method = "overrepresentation"
)

df_list_overrep <- lapply(seq_along(c2_ora_overrep), function(i) {
  df <- c2_ora_overrep[[i]]
  df$Pattern <- paste0("Pattern_", i)
  df
})

c2_df_overrep <- do.call(rbind, df_list_overrep)
c2_df_overrep$gene.set <- as.character(c2_df_overrep$gene.set)
c2_df_overrep <- c2_df_overrep[, -7]
c2_df_overrep <- as.data.frame(c2_df_overrep)

write.csv(
  c2_df_overrep,
  output_csv_overrepresentation,
  row.names = FALSE
)

message("Fichier produit : ", output_csv_overrepresentation)

