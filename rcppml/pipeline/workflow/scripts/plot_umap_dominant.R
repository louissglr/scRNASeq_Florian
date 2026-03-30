suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --------------------------------------------------
# Inputs from snakemake
# --------------------------------------------------

contrib_file <- snakemake@input[["contributions"]]
umap_file    <- snakemake@input[["umap"]]
output_pdf   <- snakemake@output[["pdf"]]

sample_name  <- snakemake@wildcards[["sample"]]
dataset      <- snakemake@params[["dataset"]]

# --------------------------------------------------
# 1. Read contributions file
# --------------------------------------------------

contrib_df <- read_tsv(contrib_file, guess_max = 1e6)

# Normalisation des noms de colonnes
colnames(contrib_df) <- tolower(colnames(contrib_df))
colnames(contrib_df) <- gsub("[^a-z0-9_]", "_", colnames(contrib_df))

# Vérification colonne barcode
barcode_candidates <- c("barcode", "barcodes", "cell", "cells", "cell_id", "cellid")
barcode_col <- intersect(colnames(contrib_df), barcode_candidates)[1]

if (is.na(barcode_col)) {
  stop("No barcode column detected in contributions file")
}

# Colonnes patterns = tout sauf sample_name / barcode / dominant_pattern si présent
exclude_cols <- c("sample_name", "dominant_pattern", barcode_col)
pattern_cols <- setdiff(colnames(contrib_df), exclude_cols)

if (length(pattern_cols) == 0) {
  stop("No pattern columns detected in contributions file")
}

# --------------------------------------------------
# 2. Ensure dominant_pattern exists (if not already computed)
# --------------------------------------------------

if (!"dominant_pattern" %in% colnames(contrib_df)) {

  contrib_df$dominant_pattern <- apply(
    contrib_df[, pattern_cols, drop = FALSE],
    1,
    function(x) pattern_cols[which.max(x)]
  )
}

# --------------------------------------------------
# 3. Clean UMAP file
# --------------------------------------------------

clean_umap <- function(file) {

  df <- read_tsv(file, guess_max = 1e6)

  colnames(df) <- tolower(colnames(df))
  colnames(df) <- gsub("[^a-z0-9_]", "_", colnames(df))

  # Barcode detection
  barcode_candidates <- c(
    "barcode", "barcodes",
    "cell", "cells",
    "cell_id", "cellid"
  )

  barcode_col <- intersect(colnames(df), barcode_candidates)[1]

  if (is.na(barcode_col)) {
    stop("No barcode-like column detected in UMAP file")
  }

  # UMAP coordinate detection
  umap_candidates <- colnames(df)[grepl("umap|dim|coord|x$|y$", colnames(df))]

  if (length(umap_candidates) < 2) {
    stop("UMAP-like columns not detected")
  }

  umap_x <- umap_candidates[1]
  umap_y <- umap_candidates[2]

  df_clean <- df |>
    rename(
      barcode = all_of(barcode_col),
      UMAP_1  = all_of(umap_x),
      UMAP_2  = all_of(umap_y)
    ) |>
    select(barcode, UMAP_1, UMAP_2, everything())

  return(df_clean)
}

umap_df <- clean_umap(umap_file)

# --------------------------------------------------
# 4. Merge UMAP + dominant pattern
# --------------------------------------------------

plot_df <- umap_df |>
  left_join(
    contrib_df |> select(all_of(barcode_col), dominant_pattern),
    by = c("barcode" = barcode_col)
  ) |>
  filter(!is.na(dominant_pattern))

# --------------------------------------------------
# 5. Plot
# --------------------------------------------------

p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = dominant_pattern)) +
  geom_point(size = 1) +
  theme_classic() +
  labs(
    title = paste0("UMAP - ", sample_name, " (", dataset, ")"),
    color = "Dominant Pattern"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# --------------------------------------------------
# 6. Export PDF
# --------------------------------------------------

pdf(output_pdf, width = 10, height = 8)
print(p)
dev.off()

message("Finished: ", output_pdf)
