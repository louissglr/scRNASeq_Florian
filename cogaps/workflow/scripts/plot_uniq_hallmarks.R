suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(viridis)
  library(tidytext)
  library(ggplot2)
})

# =========================================================
# INPUT / OUTPUT
# =========================================================
input_csv_enrichment <- snakemake@input[["csv_enrich"]]
input_csv_overr      <- snakemake@input[["csv_overr"]]
output_pdf_enrich    <- snakemake@output[["pdf_enrich"]]
output_pdf_overr     <- snakemake@output[["pdf_overr"]]

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- as.integer(snakemake@wildcards[["npatterns"]])

padj_cutoff <- 0.05
threshold   <- -10 * log10(padj_cutoff)

# =========================================================
# SAFE READER
# =========================================================
safe_read <- function(path){
  if(!file.exists(path)){
    stop(paste("Input file does not exist:", path))
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if(nrow(df) == 0){
    warning(paste("Empty dataframe:", path))
  }
  return(df)
}

# =========================================================
# GENERIC PLOT FUNCTION
# =========================================================
generate_hallmark_plots <- function(df, title_prefix){

  # Basic column checks
  required_cols <- c("padj", "gene.set", "Pattern")
  if(!all(required_cols %in% colnames(df))){
    stop("Missing required columns in dataframe.")
  }

  # Remove NA padj
  df <- df %>% filter(!is.na(padj))

  if(nrow(df) == 0){
    warning(paste0("No rows after NA removal for ", title_prefix))
    return(list(NULL, NULL, NULL))
  }

  # Filter significant
  hallmarks_sig <- df %>%
    filter(padj < padj_cutoff)

  if(nrow(hallmarks_sig) == 0){
    warning(paste0("No significant hallmarks for ", title_prefix))
    return(list(NULL, NULL, NULL))
  }

  # Clean + prepare
  hallmarks_sig <- hallmarks_sig %>%
    mutate(
      neg_log_padj = -10 * log10(padj),
      gene.set_short = str_remove(gene.set, "^HALLMARK_"),
      Pattern_num = suppressWarnings(as.integer(str_extract(Pattern, "\\d+")))
    )

  # Remove rows where Pattern_num extraction failed
  hallmarks_sig <- hallmarks_sig %>% filter(!is.na(Pattern_num))

  if(nrow(hallmarks_sig) == 0){
    warning("No valid Pattern numbers extracted.")
    return(list(NULL, NULL, NULL))
  }

  hallmarks_sig <- hallmarks_sig %>%
    mutate(
      Pattern = factor(
        Pattern,
        levels = paste0("Pattern_", sort(unique(Pattern_num)))
      )
    ) %>%
    select(-Pattern_num)

  # =====================================================
  # 1) UNIQUE
  # =====================================================
  hallmarks_unique <- hallmarks_sig %>%
    group_by(gene.set) %>%
    filter(n_distinct(Pattern) == 1) %>%
    ungroup()

  p_unique <- NULL
  if(nrow(hallmarks_unique) > 0){
    p_unique <- ggplot(
      hallmarks_unique,
      aes(
        x = reorder_within(gene.set_short, neg_log_padj, Pattern),
        y = neg_log_padj,
        fill = neg_log_padj
      )
    ) +
      geom_col() +
      geom_hline(
        yintercept = threshold,
        color = "red",
        linetype = "dashed",
        linewidth = 0.8
      ) +
      scale_fill_viridis_c(option = "C") +
      scale_x_reordered() +
      coord_flip() +
      facet_grid(Pattern ~ ., scales = "free_y", space = "free_y") +
      labs(
        title = paste0(
          "Unique Significant Hallmarks per Pattern - ",
          title_prefix, " ",
          sample_name, " (", n_patterns, " patterns)"
        ),
        x = "",
        y = expression(-10*log[10](padj))
      ) +
      theme_classic() +
      theme(
        strip.text.y = element_text(angle = 0, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.position = "right"
      )
  }

  # =====================================================
  # 2) NON-UNIQUE
  # =====================================================
  hallmarks_nonunique <- hallmarks_sig %>%
    group_by(gene.set) %>%
    filter(n_distinct(Pattern) > 1) %>%
    ungroup()

  p_nonunique <- NULL
  if(nrow(hallmarks_nonunique) > 0){
    p_nonunique <- ggplot(
      hallmarks_nonunique,
      aes(
        x = reorder_within(gene.set_short, neg_log_padj, Pattern),
        y = neg_log_padj,
        fill = neg_log_padj
      )
    ) +
      geom_col() +
      geom_hline(
        yintercept = threshold,
        color = "red",
        linetype = "dashed",
        linewidth = 0.8
      ) +
      scale_fill_viridis_c(option = "C") +
      scale_x_reordered() +
      coord_flip() +
      facet_grid(Pattern ~ ., scales = "free_y", space = "free_y") +
      labs(
        title = paste0(
          "Non-Unique Significant Hallmarks per Pattern - ",
          title_prefix, " ",
          sample_name, " (", n_patterns, " patterns)"
        ),
        x = "",
        y = expression(-10*log[10](padj))
      ) +
      theme_classic() +
      theme(
        strip.text.y = element_text(angle = 0, face = "bold"),
        axis.text.y = element_text(size = 6),
        legend.position = "right"
      )
  }

  # =====================================================
  # 3) COMBINED
  # =====================================================
  hallmarks_combined <- hallmarks_sig %>%
    group_by(gene.set) %>%
    mutate(
      uniqueness = ifelse(
        n_distinct(Pattern) == 1,
        "Unique",
        "Non-unique"
      )
    ) %>%
    ungroup()

  p_combined <- NULL
  if(nrow(hallmarks_combined) > 0){
    p_combined <- ggplot(
      hallmarks_combined,
      aes(
        x = reorder_within(gene.set_short, neg_log_padj, Pattern),
        y = neg_log_padj,
        fill = uniqueness
      )
    ) +
      geom_col() +
      geom_hline(
        yintercept = threshold,
        color = "red",
        linetype = "dashed",
        linewidth = 0.6
      ) +
      scale_x_reordered() +
      coord_flip() +
      facet_grid(Pattern ~ ., scales = "free_y", space = "free_y") +
      scale_fill_manual(
        values = c("Unique" = "red3", "Non-unique" = "grey70")
      ) +
      labs(
        title = paste0(
          "Significant Hallmarks per Pattern (Unique vs Non-Unique) - ",
          title_prefix, " ",
          sample_name, " (", n_patterns, " patterns)"
        ),
        x = "",
        y = expression(-10*log[10](padj)),
        fill = "Hallmark type"
      ) +
      theme_classic() +
      theme(
        strip.text.y = element_text(angle = 0, face = "bold"),
        axis.text.y = element_text(size = 5),
        legend.position = "right"
      )
  }

  return(list(p_unique, p_nonunique, p_combined))
}

# =========================================================
# READ CSV
# =========================================================
df_enrich <- safe_read(input_csv_enrichment)
df_overr  <- safe_read(input_csv_overr)

# =========================================================
# GENERATE PLOTS
# =========================================================
plots_enrich <- generate_hallmark_plots(df_enrich, "Enrichment")
plots_overr  <- generate_hallmark_plots(df_overr,  "Overrepresentation")

# =========================================================
# SAVE PDFs
# =========================================================
pdf(output_pdf_enrich, width = 12, height = 8)
for(p in plots_enrich){
  if(!is.null(p)) print(p)
}
dev.off()

pdf(output_pdf_overr, width = 12, height = 8)
for(p in plots_overr){
  if(!is.null(p)) print(p)
}
dev.off()
