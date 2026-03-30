suppressPackageStartupMessages({
  library(uwot)
})

contrib_file <- snakemake@input[["contrib_tsv"]]
output_tsv   <- snakemake@output[["tsv"]]

if (!dir.exists(dirname(output_tsv)))
  dir.create(dirname(output_tsv), recursive = TRUE)

# --------------------------------------------------
# Lecture
# --------------------------------------------------

contrib <- read.table(
  contrib_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# --------------------------------------------------
# Identifier les colonnes patterns
# --------------------------------------------------

pattern_cols <- setdiff(
  colnames(contrib),
  c("sample_name", "barcode", "dominant_program")
)

# Matrice numérique des patterns uniquement
pattern_df <- contrib[, pattern_cols, drop = FALSE]

# Conversion sécurisée en numeric
pattern_df[] <- lapply(pattern_df, function(x) as.numeric(as.character(x)))

if (any(is.na(pattern_df))) {
  stop("Non-numeric values detected in pattern matrix.")
}

mat <- as.matrix(pattern_df)

# Utiliser barcode comme ID cellule
cell_ids <- contrib$barcode

if (any(duplicated(cell_ids))) {
  warning("Duplicate cell IDs detected — making them unique")
  cell_ids <- make.unique(cell_ids)
}

rownames(mat) <- cell_ids

# --------------------------------------------------
# UMAP
# --------------------------------------------------

set.seed(123)

umap_res <- uwot::umap(
  mat,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine"
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df$cell_id <- rownames(mat)
umap_df <- umap_df[, c("cell_id", "UMAP1", "UMAP2")]

# --------------------------------------------------
# Export
# --------------------------------------------------

write.table(
  umap_df,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("TSV exported to: ", output_tsv)
