library(dplyr)
library(readr)
library(stringr)

files <- snakemake@input
output_path <- snakemake@output[["tsv"]]

files <- unlist(files)

all_tables <- lapply(files, function(f) {

  npattern <- str_extract(f, "npatterns-\\d+") |>
    str_remove("npatterns-")

  df <- read_tsv(f, show_col_types = FALSE)

  df <- df |>
    select(barcode, dominant_pattern) |>
    mutate(
      !!paste0("k", npattern) := str_remove(dominant_pattern, "Pattern_")
    ) |>
    select(barcode, starts_with("k"))

  return(df)
})

# Merge
cluster_df <- Reduce(function(x, y) full_join(x, y, by = "barcode"), all_tables)

cluster_df[is.na(cluster_df)] <- "NA"

k_cols <- grep("^k\\d+", colnames(cluster_df), value = TRUE)
k_order <- as.numeric(str_remove(k_cols, "k"))
k_cols <- k_cols[order(k_order)]

cluster_df <- cluster_df[, c("barcode", k_cols), drop = FALSE]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

write_tsv(cluster_df, output_path)
