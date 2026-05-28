library(readr)
library(dplyr)
library(tidyr)

rds_path    <- snakemake@input[["cogaps_rds"]]
output_path  <- snakemake@output[["tsv"]]
sample_name  <- snakemake@wildcards[["sample"]]

if (!dir.exists(dirname(output_path))) {
  dir.create(dirname(output_path), recursive = TRUE)
}

cogaps_object <- readRDS(rds_path)


sample_factors <- as.data.frame(cogaps_object@sampleFactors) %>%
  mutate(
    barcode = rownames(.),
    sample_name = sample_name
  ) %>%
  relocate(any_of(c("sample_name", "barcode")), .before = everything())

pattern_cols <- setdiff(
  colnames(sample_factors),
  c("sample_name", "barcode")
)

sample_factors$dominant_pattern <- apply(
  sample_factors[, pattern_cols, drop = FALSE],
  1,
  function(x) pattern_cols[which.max(x)]
)

sample_factors <- sample_factors %>%
  relocate(dominant_pattern, .after = last_col())

write.table(
  sample_factors,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("End: ", output_path)