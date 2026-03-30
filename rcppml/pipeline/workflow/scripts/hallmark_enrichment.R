suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
})

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
rds_path    <- snakemake@input[["rcppml_rds"]]
seurat_rds <- snakemake@input[["seurat_rds"]]
output_csv <- snakemake@output[["csv"]]
sample_id  <- snakemake@wildcards[["sample"]]

top_n <- 50

# ----------------------------
# Lecture des objets
# ----------------------------
model <- readRDS(rds_path)
seurat_obj <- readRDS(seurat_rds)

gene_names <- rownames(seurat_obj)
species_detected <- detect_species(gene_names)
message("Espèce détectée : ", species_detected)

gene_contrib <- as.data.frame(model@w)
rownames(gene_contrib) <- gene_names

# ----------------------------
# Récupération Hallmark MSigDB
# ----------------------------
msigdbr_df <- msigdbr(
  species = species_detected,
  collection = "H"
)

hallmark_t2g <- msigdbr_df |>
  dplyr::select(gs_name, gene_symbol) |>
  distinct()

# ----------------------------
# Boucle sur chaque programme
# ----------------------------
res_list <- list()

for (prog in colnames(gene_contrib)) {

  scores <- gene_contrib[, prog]
  names(scores) <- rownames(gene_contrib)

  top_genes <- names(sort(scores, decreasing = TRUE))[1:top_n]
  top_genes <- intersect(top_genes, gene_names)

  if (length(top_genes) < 5) next

  enrich <- tryCatch(
    enricher(
      gene          = top_genes,
      universe      = gene_names,
      TERM2GENE     = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1
    ),
    error = function(e) NULL
  )

  if (is.null(enrich) || nrow(enrich@result) == 0) next

  df <- enrich@result
  df$sample  <- sample_id
  df$program <- prog

  res_list[[prog]] <- df
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
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
  write.csv(empty_df, output_csv, row.names = FALSE)
} else {
  final_df <- bind_rows(res_list)
  final_df <- final_df[, c(
    "sample","program","ID","Description",
    "GeneRatio","BgRatio","pvalue","p.adjust",
    "qvalue","geneID","Count"
  )]
  final_df <- final_df[order(final_df$program, final_df$p.adjust), ]
  write.csv(final_df, output_csv, row.names = FALSE)
}
