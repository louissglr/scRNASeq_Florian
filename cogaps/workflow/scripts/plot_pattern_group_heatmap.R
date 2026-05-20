library(dplyr)
library(data.table)
library(pheatmap)

contrib_file <- snakemake@input[["contrib"]]
metadata_file <- snakemake@input[["metadata"]]
output_png <- snakemake@output[["png"]]

contrib <- fread(contrib_file)
metadata <- fread(metadata_file)


metadata <- metadata %>%
  mutate(
    group = paste(orig_tumor, tnt, Sting_status, sep = "_")
  )


contrib <- contrib %>%
  filter(barcode %in% metadata$cell_id)


df <- contrib %>%
  inner_join(
    metadata %>%
      select(cell_id, group),
    by = c("barcode" = "cell_id")
  )

pattern_cols <- setdiff(
  colnames(contrib),
  c("barcode", "sample_name", "dominant_pattern")
)


avg_scores <- df %>%
  group_by(group) %>%
  summarize(
    across(
      all_of(pattern_cols),
      \(x) mean(x, na.rm = TRUE)
    )
  )


mat <- as.matrix(avg_scores[, pattern_cols])

rownames(mat) <- avg_scores$group

# Rename columns:
# Pattern_1 -> 1
colnames(mat) <- gsub("Pattern[_]?", "", colnames(mat))


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