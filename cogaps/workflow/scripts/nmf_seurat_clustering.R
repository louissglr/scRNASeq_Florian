library(Seurat)
library(dplyr)
library(readr)
library(tibble)

contrib_path <- snakemake@input[["contrib"]]
output_path  <- snakemake@output[["clusters"]]
sample_name  <- snakemake@wildcards[["sample"]]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load contributions matrix
# ---------------------------
df <- read_tsv(contrib_path, show_col_types = FALSE)

# barcode = rownames
df <- df %>%
  column_to_rownames("barcode")

# remove metadata columns (IMPORTANT)
df <- df %>%
  select(-sample_name, -dominant_pattern)

# NMF matrix (cells x patterns)
A <- as.matrix(df)
k <- ncol(A)

# ---------------------------
# Seurat object + embedding
# ---------------------------
seu <- CreateSeuratObject(counts = t(A))

embeddings <- A
colnames(embeddings) <- paste0("nmf_f", 1:k)

seu[["nmf"]] <- CreateDimReducObject(
  embeddings = embeddings,
  key = "nmf_"
)

# ---------------------------
# Graph + UMAP + clustering
# ---------------------------
seu <- FindNeighbors(
  seu,
  reduction = "nmf",
  dims = 1:k,
  k.param = 20,
  annoy.metric = "cosine"
)

seu <- RunUMAP(
  seu,
  reduction = "nmf",
  dims = 1:k,
  n.neighbors = 50
)

seu <- FindClusters(
  seu,
  resolution = 0.7
)

# ---------------------------
# Output
# ---------------------------
out <- data.frame(
  cell = rownames(seu@meta.data),
  cluster = seu$seurat_clusters
)

write.table(
  out,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done: ", output_path)
