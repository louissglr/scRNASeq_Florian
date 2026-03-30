library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

markers_file <- snakemake@input[["markers"]]
pdf_file     <- snakemake@output[["pdf"]]

markers <- as.data.frame(read_xlsx(markers_file))

# Number of genes to keep
topN <- 200
markers_top <- markers[1:topN, ]

# Convert each pattern column to gene vector
lst <- lapply(markers_top, unlist)

ego_list <- list()

# GO enrichment per pattern
for (pattern_name in names(lst)) {

    genes <- lst[[pattern_name]]

    ego <- enrichGO(
        gene          = genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.05
    )

    ego_list[[pattern_name]] <- ego
}

# Merge results across patterns
merged <- merge_result(ego_list)

# Generate PDF
pdf(pdf_file, width = 10, height = 8)
print(
    dotplot(
        merged,
        showCategory = 10,
        font.size = 8,
        label_format = 30,
        title = "GO Enrichment merged"
    )
)
dev.off()
