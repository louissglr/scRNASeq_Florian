library(dplyr)
library(tidyr)

rds_path <- snakemake@input[["rcppml_rds"]]
seurat_rds <- snakemake@input[["seurat_rds"]]
output_path <- snakemake@output[["tsv"]]
sample_name <- snakemake@wildcards[["sample"]]

model <- readRDS(rds_path)

seurat_obj <- readRDS(seurat_rds)
barcodes <- colnames(seurat_obj)

cell_contrib <- as.data.frame(t(model@h))

colnames(cell_contrib) <- paste0("programs_", 1:ncol(cell_contrib))

cell_contrib$barcode <- barcodes
cell_contrib$sample_name <- sample_name

cell_contrib <- cell_contrib %>%
  relocate(any_of(c("sample_name", "barcode")), .before = everything())

nmf_cols <- setdiff(
  colnames(cell_contrib),
  c("sample_name", "barcode")
)

cell_contrib$dominant_program <- apply(
  cell_contrib[, nmf_cols, drop = FALSE],
  1,
  function(x) nmf_cols[which.max(x)]
)

cell_contrib <- cell_contrib %>%
  relocate(dominant_program, .after = barcode)

if (!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)

write.table(cell_contrib, file = output_path, sep = "\t", row.names = FALSE)

message("End: ", output_path)