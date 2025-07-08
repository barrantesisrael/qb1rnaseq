#!/usr/bin/env Rscript

# Load libraries
library(DESeq2)
library(ggplot2)

# Create sample information
samples <- data.frame(
    sample = c("cancer_rep1", "cancer_rep2", "cancer_rep3", 
               "normal_rep1", "normal_rep2", "normal_rep3"),
    condition = c("cancer", "cancer", "cancer", 
                  "normal", "normal", "normal")
)

# Read salmon results
files <- file.path("salmon_quant", samples$sample, "quant.sf")
names(files) <- samples$sample

# Create count matrix
first_file <- read.table(files[1], header = TRUE, sep = "\t")
genes <- first_file$Name
count_matrix <- matrix(0, nrow = length(genes), ncol = length(files))
rownames(count_matrix) <- genes
colnames(count_matrix) <- names(files)

for(i in 1:length(files)) {
    data <- read.table(files[i], header = TRUE, sep = "\t")
    count_matrix[, i] <- round(data$NumReads)
}

cat("Analyzing", nrow(count_matrix), "genes across", ncol(count_matrix), "samples\n")

# Save data for next steps
save(count_matrix, samples, file = "data_loaded.RData")
