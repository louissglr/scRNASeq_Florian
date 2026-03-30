suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(writexl)
})

seurat_path  <- snakemake@input[["seurat_rds"]]
contrib_path <- snakemake@input[["contrib_tsv"]]
output_path  <- snakemake@output[["excel"]]  


seurat_obj <- readRDS(seurat_path)
contrib <- read_tsv(contrib_path, show_col_types = FALSE)


common_cells <- intersect(contrib$barcode, colnames(seurat_obj))

contrib <- contrib |>
  filter(barcode %in% common_cells)

seurat_obj <- subset(seurat_obj, cells = common_cells)
contrib <- contrib[match(colnames(seurat_obj), contrib$barcode), ]

seurat_obj$nmf_program <- contrib$dominant_pattern
Idents(seurat_obj) <- "nmf_program"

seurat_obj <- NormalizeData(seurat_obj)

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)


top_markers <- markers |>
  filter(p_val_adj < 0.05) |>                  
  group_by(cluster) |>
  filter(n() > 0) |>                           
  arrange(p_val_adj, .by_group = TRUE) |>      
  ungroup() |>
  select(gene, everything())                    


split_markers <- split(top_markers, top_markers$cluster)


#names(split_markers) <- paste0("Program_", names(split_markers))

write_xlsx(split_markers, path = output_path)