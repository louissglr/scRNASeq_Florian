suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readr)
})

set.seed(123)
# ----------------------------
# Détection automatique de l'espèce
# ----------------------------
detect_species <- function(gene_names) {
  gene_names <- unique(gene_names)
  human_prop <- mean(gene_names == toupper(gene_names), na.rm = TRUE)
  mouse_prop <- mean(grepl("^[A-Z][a-z0-9]+$", gene_names), na.rm = TRUE)
  
  if (human_prop > 0.8) {
    return("Homo sapiens")
  } else if (mouse_prop > 0.8) {
    return("Mus musculus")
  } else {
    message("Espèce ni humain ni souris, par défaut : Homo sapiens")
    return("Homo sapiens")
  }
}

# ----------------------------
# Snakemake inputs/outputs
# ----------------------------
contrib_path <- snakemake@input[["contrib_txt"]]  # fichier TXT tab-delimited
output_csv   <- snakemake@output[["csv"]]
sample_id    <- snakemake@wildcards[["sample"]]

# ----------------------------
# Lecture du fichier contrib
# ----------------------------
contrib_df <- read.delim(contrib_path, check.names = FALSE)  # TXT tab-delimited

# La première colonne doit être "gene"
rownames(contrib_df) <- contrib_df$gene
contrib_df$gene <- NULL

# Renommer les colonnes pour avoir program_1, program_2, ...
colnames(contrib_df) <- paste0("program_", seq_len(ncol(contrib_df)))

# Récupérer tous les gènes uniques pour détection d'espèce
gene_names <- rownames(contrib_df)
species_detected <- detect_species(gene_names)
message("Espèce détectée : ", species_detected)

# ----------------------------
# Récupération Hallmark MSigDB
# ----------------------------
msigdbr_df <- msigdbr(
  species = species_detected,
  collection = "H"
)

hallmark_t2g <- msigdbr_df %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

# ----------------------------
# Boucle sur chaque programme pour GSEA
# ----------------------------
res_list <- list()

for (prog in colnames(contrib_df)) {

  # Création du vecteur nommé pour GSEA
  gene_scores <- contrib_df[[prog]]
  names(gene_scores) <- rownames(contrib_df)

  # Supprimer les NA et trier décroissant
  gene_scores <- gene_scores[!is.na(gene_scores)]
  gene_scores <- sort(gene_scores, decreasing = TRUE)

  # Lancer GSEA
  gsea_res <- tryCatch(
    GSEA(
      geneList    = gene_scores,
      TERM2GENE   = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,  # on récupère tout pour filtrer plus tard
      verbose       = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    df <- gsea_res@result
    df$sample  <- sample_id
    df$program <- prog
    res_list[[prog]] <- df
  }
}

# ----------------------------
# Création CSV final
# ----------------------------
if (length(res_list) == 0) {
  empty_df <- data.frame(
    sample = character(),
    program = character(),
    ID = character(),
    Description = character(),
    setSize = integer(),
    enrichmentScore = numeric(),
    NES = numeric(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalues = numeric(),
    rank = character()
  )
  write.csv(empty_df, output_csv, row.names = FALSE)
} else {
  final_df <- bind_rows(res_list)
  final_df <- final_df[order(final_df$program, final_df$p.adjust), ]
  write.csv(final_df, output_csv, row.names = FALSE)
}

message("GSEA terminée et résultats exportés dans : ", output_csv)