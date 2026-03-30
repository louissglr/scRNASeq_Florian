# Logging
if (exists("snakemake")) {
  log_file <- snakemake@log[[1]]           
  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  mainlog <- file(log_file, open = "wt")
  sink(mainlog, append = FALSE, type = "output")
  sink(mainlog, append = FALSE, type = "message")

  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(close(mainlog), add = TRUE)
}

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(Matrix)
    library(CoGAPS)
})

# Retrieve inputs and parameters from Snakemake
input_file <- snakemake@input[["tumor_rna_seuratobj"]]
output_file <- snakemake@output[["cogaps_rds"]]
sample_name <- snakemake@wildcards[["sample"]]
n_patterns <- as.integer(snakemake@wildcards[["npatterns"]])
n_iterations <- as.integer(snakemake@params[["nIterations"]]) 

message(paste(">>> Sample:", sample_name, "| nPatterns:", n_patterns))

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load Seurat object and extract expression matrix
seurat_obj <- readRDS(input_file)
expr_data <- as.matrix(seurat_obj@assays$RNA$counts)
norm_expr_data <- log1p(expr_data)

# Configure CoGAPS parameters
cogaps_params <- CogapsParams(
    nIterations = n_iterations,     
    seed = 42,
    nPatterns = n_patterns,
    sparseOptimization = TRUE,
    distributed = "genome-wide"
)

# Distributed runs
cogaps_params <- setDistributedParams(cogaps_params, nSets = 7)

# Run CoGAPS
message(">>> Starting CoGAPS...")
cogaps_result <- CoGAPS(norm_expr_data, cogaps_params)

# Save result
saveRDS(cogaps_result, output_file)
