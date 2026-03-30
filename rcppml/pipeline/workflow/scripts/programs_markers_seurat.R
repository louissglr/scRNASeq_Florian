suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(writexl)
})

# -------------------
# Inputs Snakemake
# -------------------
seurat_path  <- snakemake@input[["seurat_rds"]]
contrib_path <- snakemake@input[["contrib_tsv"]]
output_path  <- snakemake@output[["xlsx"]]  # fichier Excel

top_n <- 50

# -------------------
# Load data
# -------------------
seurat_obj <- readRDS(seurat_path)
contrib <- read_tsv(contrib_path)

# -------------------
# Vérification critique ⚠️
# -------------------
# barcodes doivent matcher les cellules Seurat
common_cells <- intersect(contrib$barcode, colnames(seurat_obj))

if (length(common_cells) == 0) {
  stop("Aucune correspondance entre barcodes et Seurat object")
}

# subset + ordre cohérent
contrib <- contrib %>%
  filter(barcode %in% common_cells)

seurat_obj <- subset(seurat_obj, cells = common_cells)

# remettre dans le même ordre
contrib <- contrib[match(colnames(seurat_obj), contrib$barcode), ]

# -------------------
# Assignation identités
# -------------------
seurat_obj$nmf_program <- contrib$dominant_program
Idents(seurat_obj) <- "nmf_program"

seurat_obj <- NormalizeData(seurat_obj)

# -------------------
# Find markers
# -------------------
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# -------------------
# Top genes par programme
# -------------------
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n)

# -------------------
# Export Excel multi-feuilles
# -------------------
# Créer une liste avec 1 élément par programme
split_markers <- split(top_markers, top_markers$cluster)

# Optionnel : noms de feuilles plus lisibles
names(split_markers) <- paste0("Program_", names(split_markers))

# Écrire le fichier Excel
write_xlsx(split_markers, path = output_path)