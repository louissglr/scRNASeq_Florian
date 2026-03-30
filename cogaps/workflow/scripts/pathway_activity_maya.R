## ============================================================
## MAYA Pathway Activity Workflow - Adapted for Metadata
## Snakemake Compatible
## ============================================================

if (!requireNamespace("MAYA", quietly = TRUE)) {
  remotes::install_github("one-biosciences/maya")
}
library(MAYA)
library(Seurat)
library(msigdbr)
library(dplyr)
library(tibble)

## ============================================================
## Inputs / Outputs
## ============================================================

seurat_file    <- snakemake@input$seurat
metadata_file  <- snakemake@input$metadata  # Nouveau fichier metadata

heatmap1   <- snakemake@output$heatmap1
heatmap2   <- snakemake@output$heatmap2
umap_file  <- snakemake@output$umap
top_pdf    <- snakemake@output$top_genes
spec_pdf   <- snakemake@output$specificity

## ============================================================
## Load Seurat Object
## ============================================================

message("Loading Seurat object")

seurat_obj <- tryCatch(
  readRDS(seurat_file),
  error = function(e) stop("Failed to load Seurat object: ", e$message)
)

count_mat <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer = "counts"
)

message(
  "Matrix dimension: ",
  nrow(count_mat), " genes x ",
  ncol(count_mat), " cells"
)

## ============================================================
## Gene sets
## ============================================================

message("Loading MSigDB gene sets")

genesets_mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

modules_list_mouse <- split(
  genesets_mouse$gene_symbol,
  genesets_mouse$gs_name
)

rownames(count_mat) <- toupper(rownames(count_mat))
modules_list_mouse  <- lapply(modules_list_mouse, toupper)

common_genes <- intersect(
  rownames(count_mat),
  unlist(modules_list_mouse)
)

message("Common genes between matrix and gene sets: ", length(common_genes))

## ============================================================
## Load Metadata
## ============================================================

message("Loading metadata")

metadata_raw <- read.table(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

message("Metadata rows: ", nrow(metadata_raw))

## check barcode overlap
seurat_cells <- colnames(seurat_obj)

common_cells <- intersect(
  metadata_raw$cell_id,
  seurat_cells
)

message("Cells in Seurat: ", length(seurat_cells))
message("Cells in metadata: ", nrow(metadata_raw))
message("Matching cells: ", length(common_cells))

if (length(common_cells) == 0) {
  stop("No matching barcodes between metadata and Seurat object")
}

metadata <- metadata_raw |>
  filter(cell_id %in% common_cells) |>
  mutate(
    group = paste(orig_tumor, tnt, Sting_status, sep = "_")
  ) |>
  column_to_rownames("cell_id")

metadata <- metadata[colnames(seurat_obj), , drop = FALSE]

message("Metadata aligned with Seurat object")

## ============================================================
## Pathway Activity Analysis
## ============================================================

message("Running pathway analysis")

activity_summary <- MAYA_pathway_analysis(
  expr_mat = count_mat,
  modules_list = modules_list_mouse,
  is_logcpm = FALSE
)

activity_mat_scale <- scale_0_1(
  activity_summary$activity_matrix
)

modules <- names(activity_summary$PCA_obj)

## ============================================================
## Heatmap 1
## ============================================================

message("Generating heatmap 1")

png(
  heatmap1,
  width = 1200,
  height = 1000,
  res = 150
)

plot_heatmap_activity_mat(
  activity_mat = activity_mat_scale,
  meta = metadata,
  annot_name = "group",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  clustering_distance = "euclidean",
  clustering_method = "ward.D2"
)

dev.off()

## ============================================================
## Heatmap 2
## ============================================================

message("Generating heatmap 2")

png(
  heatmap2,
  width = 1200,
  height = 1000,
  res = 150
)

plot_heatmap_activity_mat(
  activity_mat = activity_mat_scale,
  meta = metadata,
  annot_name = "group",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance = "euclidean",
  clustering_method = "ward.D2"
)

dev.off()

## ============================================================
## UMAP
## ============================================================

message("Generating UMAP")

png(
  umap_file,
  width = 1000,
  height = 800,
  res = 150
)

plot_umap_annot(
  umap = activity_summary$umap,
  labels = as.factor(metadata$group),
  title = "Combined annotation - HALLMARK"
)

dev.off()

## ============================================================
## LogCPM normalization
## ============================================================

message("Running logCPM normalization")

options(future.globals.maxSize = 8 * 1024^3)

logcpm <- logcpmNormalization(count_mat)

## ============================================================
## Top contributing genes
## ============================================================

message("Plotting top contributing genes")

pdf(
  top_pdf,
  width = 10,
  height = 8
)

for (mod in modules) {
  
  try({
    
    plot_heatmap_pathway_top_contrib_genes(
      expr_mat = logcpm,
      PCA_object = activity_summary$PCA_obj,
      module = mod,
      n = 20,
      meta = metadata,
      annot_name = "group"
    )
    
  }, silent = TRUE)
  
}

dev.off()

## ============================================================
## Pathway specificity
## ============================================================

message("Plotting pathway specificity")

pdf(
  spec_pdf,
  width = 8,
  height = 6
)

for (mod in modules) {
  
  try({
    
    plot_pathway_specificity(
      PCA_object = activity_summary$PCA_obj,
      module = mod,
      meta = metadata,
      annot_name = "group"
    )
    
  }, silent = TRUE)
  
}

dev.off()

message("Workflow completed successfully")
