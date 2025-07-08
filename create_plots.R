#!/usr/bin/env Rscript

# Load results
load("deseq2_results.RData")
library(ggplot2)

# PCA plot
vsd <- vst(dds)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("cancer" = "red", "normal" = "blue")) +
    labs(title = "PCA: Cancer vs Normal Samples") +
    theme_minimal()

ggsave("pca_plot.png", p1, width = 8, height = 6)

# Heatmap of top genes
top_genes <- head(order(results$padj), 20)
heatmap_data <- assay(vsd)[top_genes, ]

png("heatmap.png", width = 800, height = 600)
heatmap(heatmap_data, main = "Top 20 Cancer-Associated Genes", scale = "row")
dev.off()
