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


feature_loadings <- as.data.frame(cogaps_object@featureLoadings) %>%
  mutate(
    gene = rownames(.)
  ) %>%
  relocate(any_of(c("gene")), .before = everything())


write.table(
  feature_loadings,
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("End: ", output_path)
