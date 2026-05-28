suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

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
      select(cell_id, group, orig_tumor, tnt, Sting_status),
    by = c("barcode" = "cell_id")
  )


pattern_cols <- setdiff(
  colnames(contrib),
  c("barcode", "sample_name", "dominant_pattern")
)



avg_scores <- df %>%
  group_by(group, orig_tumor, tnt, Sting_status) %>%
  summarize(
    across(all_of(pattern_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )


avg_scores <- avg_scores %>%
  mutate(
    orig_tumor = trimws(tolower(orig_tumor)),
    tnt = trimws(tnt),
    tnt = gsub("\\.", "", tnt)
  )

avg_scores$orig_tumor <- as.character(avg_scores$orig_tumor)
avg_scores$tnt <- as.character(avg_scores$tnt)

# Matrix

mat <- as.matrix(avg_scores[, pattern_cols])
rownames(mat) <- avg_scores$group

colnames(mat) <- gsub("Pattern[_]?", "", colnames(mat))

# Scale rows

mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

# Clustering

row_clust <- hclust(dist(mat_scaled), method = "ward.D2")
col_clust <- hclust(dist(t(mat_scaled)), method = "ward.D2")

# Colors

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

# Annotation colors

annot_colors <- list(
  Tumor = c(
    meta = "#D55E00",
    early = "#0072B2",
    late = "#009E73"
  ),
  Treatment = c(
    NT = "#999999",
    Tax = "#CC79A7"
  )
)

# Row annotation

left_annot <- rowAnnotation(
  
  Tumor = avg_scores$orig_tumor,
  Treatment = avg_scores$tnt,
  
  col = annot_colors,
  
  show_annotation_name = TRUE,
  
  annotation_name_gp = grid::gpar(
    fontsize = 12,
    fontface = "bold"
  )
)

# ----------------------------
# Heatmap
# ----------------------------

ht <- Heatmap(
  
  mat_scaled,
  
  name = "Z-score",
  
  col = col_fun,
  
  cluster_rows = as.dendrogram(row_clust),
  cluster_columns = as.dendrogram(col_clust),
  
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  
  column_names_rot = 0,  
  column_names_gp = grid::gpar(fontsize = 13),
  
  border = FALSE,
  
  heatmap_legend_param = list(
    title = "Scaled\nscore"
  )
)

# Save

png(
  output_png,
  width = 1800,
  height = 1000,
  res = 200
)

draw(
  left_annot + ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()
