library(readr)
library(dplyr)
library(tidyr)
library(tibble)

rds_path    <- snakemake@input[["cogaps_rds"]]
output_path <- snakemake@output[["tsv"]]
sample_name <- snakemake@wildcards[["sample"]]

if (!dir.exists(dirname(output_path))) {
  dir.create(dirname(output_path), recursive = TRUE)
}

cogaps_object <- readRDS(rds_path)

sample_factors <- as.data.frame(cogaps_object@sampleFactors) |>
  rownames_to_column("barcode") |>
  mutate(sample_name = sample_name) |>
  relocate(sample_name, barcode) 

pattern_cols <- setdiff(
  colnames(sample_factors),
  c("sample_name", "barcode")
)

sample_factors <- sample_factors |>
  mutate(
    dominant_pattern = pattern_cols[
      max.col(across(all_of(pattern_cols)), ties.method = "first")
    ]
  ) |>
  relocate(dominant_pattern, .after = barcode)

write.table(
  sample_factors,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("End: ", output_path)