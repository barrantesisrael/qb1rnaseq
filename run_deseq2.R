#!/usr/bin/env Rscript

# Load previous data
load("data_loaded.RData")
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = samples,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "normal")

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run analysis
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "cancer", "normal"))

# Show results
summary(results)

# Get significant genes
oncogenes <- subset(results, !is.na(padj) & padj < 0.05 & log2FoldChange > 1)
tumor_suppressors <- subset(results, !is.na(padj) & padj < 0.05 & log2FoldChange < -1)

cat("Oncogenes (upregulated in cancer):", nrow(oncogenes), "\n")
cat("Tumor suppressors (downregulated in cancer):", nrow(tumor_suppressors), "\n")

# Save results
save(dds, results, oncogenes, tumor_suppressors, file = "deseq2_results.RData")
