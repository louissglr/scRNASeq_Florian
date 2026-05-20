library(dplyr)
library(data.table)
library(pheatmap)

contrib_file <- snakemake@input[["contrib"]]
cluster_file <- snakemake@input[["clusters"]]
output_png   <- snakemake@output[["png"]]

# ---------------------------
# Load data
# ---------------------------
contrib <- fread(contrib_file)
clusters <- fread(cluster_file)

# clusters: cell | cluster
setnames(clusters, c("cell", "cluster"))

# ---------------------------
# Merge contributions + clusters
# ---------------------------
df <- contrib %>%
  inner_join(
    clusters,
    by = c("barcode" = "cell")
  )

# ---------------------------
# Pattern columns
# ---------------------------
pattern_cols <- setdiff(
  colnames(contrib),
  c("barcode", "sample_name", "dominant_pattern")
)

# ---------------------------
# Average per NMF cluster
# ---------------------------
avg_scores <- df %>%
  group_by(cluster) %>%
  summarize(
    across(
      all_of(pattern_cols),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# ---------------------------
# Matrix for heatmap
# ---------------------------
mat <- as.matrix(avg_scores[, pattern_cols])
rownames(mat) <- paste0("cluster_", avg_scores$cluster)

colnames(mat) <- gsub("Pattern[_]?", "", colnames(mat))

# ---------------------------
# Heatmap
# ---------------------------
png(
  output_png,
  width = 1800,
  height = 1000,
  res = 200
)

pheatmap(
  mat,
  scale = "row",
  clustering_method = "ward.D2",
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 14,
  angle_col = 0
)

dev.off()