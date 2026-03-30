library(CoGAPS)
library(openxlsx)

# Inputs/outputs Snakemake
cogaps_file  <- snakemake@input[["cogaps_rds"]]
output_file  <- snakemake@output[["excel"]]
sample_name  <- snakemake@wildcards[["sample"]]
n_patterns   <- snakemake@wildcards[["npatterns"]]

message("Traitement de ", cogaps_file)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Lecture de l'objet CoGAPS
cogaps <- readRDS(cogaps_file)

# =========================
# PatternMarkers (top gènes par pattern)
# =========================
patternMarkerResults <- patternMarkers(cogaps, threshold = "cut")

max_len <- max(lengths(patternMarkerResults$PatternMarkers))
df_markers <- as.data.frame(
  lapply(patternMarkerResults$PatternMarkers, function(x) {
    c(x, rep(NA, max_len - length(x)))
  }),
  stringsAsFactors = FALSE
)
colnames(df_markers) <- paste0("Pattern_", seq_along(df_markers))

# Ajout de la colonne Rank
df_markers <- cbind(Rank = seq_len(nrow(df_markers)), df_markers)

# =========================
# PatternScores (score de chaque gène pour chaque pattern)
# =========================
df_scores <- as.data.frame(patternMarkerResults$PatternScores)
df_scores <- cbind(Gene = rownames(df_scores), df_scores)
df_scores <- df_scores[order(apply(df_scores[,-1], 1, min)), ]


# =========================
# Création fichier Excel avec deux onglets
# =========================
wb <- createWorkbook()

# Onglet 1 : Markers
sheet_name_markers <- substr(sprintf("S_%s_P%s_markers", sample_name, n_patterns), 1, 31)
addWorksheet(wb, sheet_name_markers)
writeData(wb, sheet = sheet_name_markers, x = df_markers)

# Onglet 2 : Scores
sheet_name_scores <- substr(sprintf("S_%s_P%s_scores", sample_name, n_patterns), 1, 31)
addWorksheet(wb, sheet_name_scores)
writeData(wb, sheet = sheet_name_scores, x = df_scores)

# Sauvegarde
saveWorkbook(wb, output_file, overwrite = TRUE)
message("Fichier Excel sauvegardé : ", output_file)
