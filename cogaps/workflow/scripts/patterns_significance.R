suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

# ---- Inputs ----
contributions_file <- snakemake@input[["contributions_file"]]
kw_output_csv <- snakemake@output[["kw_csv"]]
pairwise_output_csv <- snakemake@output[["pairwise_csv"]]

# ---- Metadata ----
sample_name <- snakemake@wildcards[["sample"]]
n_patterns <- as.integer(snakemake@wildcards[["npatterns"]])

# ---- Read contributions ----
contrib <- read.table(
  contributions_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) |>
  mutate(barcode = sub("^.*#", "", barcode))

# ---- Long format (patterns only) ----
contrib_long <- contrib |>
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "Pattern",
    values_to = "Weight"
  )

# ---- Kruskal-Wallis : comparaison des patterns ----
kw_results <- contrib_long |>
  kruskal_test(Weight ~ Pattern) |>
  adjust_pvalue(method = "BH") |>
  add_significance("p.adj") |>
  mutate(
    sample = sample_name,
    n_patterns = n_patterns
  )

write.csv(kw_results, kw_output_csv, row.names = FALSE)
message("Kruskal-Wallis results written to: ", kw_output_csv)

# ---- Pairwise Wilcoxon entre patterns ----
pairwise_results <- pairwise.wilcox.test(
  x = contrib_long$Weight,
  g = contrib_long$Pattern,
  p.adjust.method = "BH"
)

pairwise_df <- as.data.frame(as.table(pairwise_results$p.value)) |>
  setNames(c("group1", "group2", "p.adj")) |>
  drop_na(p.adj) |>
  mutate(
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    sample = sample_name,
    n_patterns = n_patterns
  )

write.csv(pairwise_df, pairwise_output_csv, row.names = FALSE)
message("Pairwise Wilcoxon results written to: ", pairwise_output_csv)
