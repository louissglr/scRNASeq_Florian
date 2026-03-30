# -----------------------------
# Inputs Snakemake
# -----------------------------
cogaps_files <- snakemake@input[["cogaps_files"]]
output_csv <- snakemake@output[["csv"]]

# -----------------------------
# Libraries
# -----------------------------
library(CoGAPS)
library(dplyr)
library(stringr)

# -----------------------------
# Loop over all files and extract info
# -----------------------------
recap_list <- lapply(cogaps_files, function(file) {
  
  # Extract sample name and npatterns from filename
  fname <- basename(file)
  # Exemple: ALL_t_variable.npatterns-10.cogaps-object.Rds
  parts <- str_match(fname, "(.*)\\.npatterns-(\\d+)\\.cogaps-object\\.Rds")
  sample_name <- parts[2]
  expected_patterns <- as.integer(parts[3])
  
  # Load CoGAPS object
  cogaps <- readRDS(file)
  obtained_patterns <- ncol(cogaps@sampleFactors)
  
  mean_chisq <- getMeanChiSq(cogaps)

  data.frame(
    sample = sample_name,
    expected_patterns = expected_patterns,
    obtained_patterns = obtained_patterns,
    meanChisqValue = mean_chisq
  )
})

# -----------------------------
# Combine all runs
# -----------------------------
recap_table <- bind_rows(recap_list)

# -----------------------------
# Save
# -----------------------------
write.csv(recap_table, output_csv, row.names = FALSE)
