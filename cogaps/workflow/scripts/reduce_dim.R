suppressPackageStartupMessages({
  library(uwot)
})

contrib_file <- snakemake@input[["contrib_tsv"]]
output_tsv   <- snakemake@output[["tsv"]]

if (!dir.exists(dirname(output_tsv)))
  dir.create(dirname(output_tsv), recursive = TRUE)

contrib <- read.table(
  contrib_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pattern_cols <- setdiff(
  colnames(contrib),
  c("sample_name", "barcode", "dominant_pattern")
)

pattern_df <- contrib[, pattern_cols, drop = FALSE]
pattern_df[] <- lapply(pattern_df, function(x) as.numeric(as.character(x)))

mat <- as.matrix(pattern_df)
cell_ids <- contrib$barcode

if (any(duplicated(cell_ids))) {
  warning("Duplicate cell IDs detected — making them unique")
  cell_ids <- make.unique(cell_ids)
}

rownames(mat) <- cell_ids

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

write.table(
  umap_df,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("TSV exported to: ", output_tsv)
