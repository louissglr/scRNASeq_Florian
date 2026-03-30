suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(readr)
})

set.seed(123)

# ----------------------------
# PARAMETRES
# ----------------------------
top_n <- 100

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
contrib_path   <- snakemake@input[["contrib_txt"]]
output_gsea    <- snakemake@output[["gsea_csv"]]
output_ora     <- snakemake@output[["ora_csv"]]
sample_id      <- snakemake@wildcards[["sample"]]

# ----------------------------
# Lecture du fichier contrib
# ----------------------------
contrib_df <- read.delim(contrib_path, check.names = FALSE)

rownames(contrib_df) <- contrib_df$gene
contrib_df$gene <- NULL

gene_names <- rownames(contrib_df)
species_detected <- detect_species(gene_names)
message("Espèce détectée : ", species_detected)

# ----------------------------
# MSigDB Hallmark
# ----------------------------
msigdbr_df <- msigdbr(
  species = species_detected,
  collection = "H"
)

hallmark_t2g <- msigdbr_df %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

# ----------------------------
# LISTES RESULTATS
# ----------------------------
gsea_list <- list()
ora_list  <- list()

# ----------------------------
# BOUCLE PRINCIPALE
# ----------------------------
for (prog in colnames(contrib_df)) {

  gene_scores <- contrib_df[[prog]]
  names(gene_scores) <- rownames(contrib_df)

  gene_scores <- gene_scores[!is.na(gene_scores)]
  gene_scores <- sort(gene_scores, decreasing = TRUE)

  # ============================
  # GSEA
  # ============================
  gsea_res <- tryCatch(
    GSEA(
      geneList      = gene_scores,
      TERM2GENE     = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    df <- gsea_res@result
    df$sample  <- sample_id
    df$program <- prog
    gsea_list[[prog]] <- df
  }

  # ============================
  # ORA (top N genes)
  # ============================
  top_genes <- names(gene_scores)[1:min(top_n, length(gene_scores))]

  ora_res <- tryCatch(
    enricher(
      gene          = top_genes,
      universe      = names(gene_scores),  
      TERM2GENE     = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1
    ),
    error = function(e) NULL
  )

  if (!is.null(ora_res) && nrow(ora_res@result) > 0) {
    df <- ora_res@result
    df$sample  <- sample_id
    df$program <- prog
    ora_list[[prog]] <- df
  }
}

# ----------------------------
# EXPORT GSEA
# ----------------------------
if (length(gsea_list) == 0) {
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
  write.csv(empty_df, output_gsea, row.names = FALSE)
} else {
  gsea_df <- bind_rows(gsea_list)
  gsea_df <- gsea_df[order(gsea_df$program, gsea_df$p.adjust), ]
  write.csv(gsea_df, output_gsea, row.names = FALSE)
}

# ----------------------------
# EXPORT ORA
# ----------------------------
if (length(ora_list) == 0) {
  empty_df <- data.frame(
    sample = character(),
    program = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
  write.csv(empty_df, output_ora, row.names = FALSE)
} else {
  ora_df <- bind_rows(ora_list)
  ora_df <- ora_df[order(ora_df$program, ora_df$p.adjust), ]
  write.csv(ora_df, output_ora, row.names = FALSE)
}

# ----------------------------
# LOG FINAL
# ----------------------------
message("GSEA exportée : ", output_gsea)
message("ORA exportée : ", output_ora)