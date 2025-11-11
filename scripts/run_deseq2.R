#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("--metadata", "-m"), type="character", default=NULL,
              help="Sample metadata file (TSV format)"),
  make_option(c("--expression", "-e"), type="character", default=NULL,
              help="Count matrix file (TSV format)"),
  make_option(c("--output-degs", "-d"), type="character", default="degs.tsv",
              help="Output file for differentially expressed genes [default %default]"),
  make_option(c("--output-plots", "-p"), type="character", 
              default="pca_plot.png,heatmap.png,volcano_plot.png",
              help="Output plot files (comma-separated) [default %default]"),
  make_option(c("--output-image", "-r"), type="character", default="deseq2_results.RData",
              help="Output R workspace file [default %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list,
                          description="Run DESeq2 differential expression analysis")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$metadata)) {
  stop("Error: --metadata argument is required")
}
if (is.null(opt$expression)) {
  stop("Error: --expression argument is required")
}

# Read input files
cat("Reading input files...\n")
samples <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)
count_matrix <- read.table(opt$expression, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Ensure sample order matches between metadata and count matrix
samples <- samples[match(colnames(count_matrix), samples$sample), ]

# Create DESeq2 dataset
suppressWarnings({
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = samples,
                                design = ~ condition)
})

# Set reference level
dds$condition <- relevel(dds$condition, ref = "normal")

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Retained", nrow(dds), "genes for analysis\n")

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
suppressWarnings({
  suppressMessages({
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "cancer", "normal"))
  })
})

# Extract gene symbols from transcript/gene IDs
extract_gene_symbols <- function(transcript_ids) {
  gene_symbols <- rep(NA, length(transcript_ids))
  
  if (file.exists("refs/22.gtf")) {
    cat("Extracting gene symbols from GTF file...\n")
    
    # Read GTF file
    gtf_lines <- readLines("refs/22.gtf")
    gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
    
    # Create a lookup table for transcript_id to gene_name
    transcript_to_gene <- list()
    
    for (line in gtf_lines) {
      # Extract transcript_id and gene_name from attributes
      transcript_match <- regexpr('transcript_id "([^"]+)"', line)
      gene_name_match <- regexpr('gene_name "([^"]+)"', line)
      
      if (transcript_match > 0 && gene_name_match > 0) {
        transcript_id <- regmatches(line, transcript_match)
        transcript_id <- gsub('transcript_id "([^"]+)"', '\\1', transcript_id)
        
        gene_name <- regmatches(line, gene_name_match)
        gene_name <- gsub('gene_name "([^"]+)"', '\\1', gene_name)
        
        # Remove version number from transcript ID for matching
        clean_transcript_id <- sub("\\.[0-9]+$", "", transcript_id)
        transcript_to_gene[[clean_transcript_id]] <- gene_name
      }
    }
    
    # Map transcript IDs to gene symbols
    for (i in 1:length(transcript_ids)) {
      transcript_id <- transcript_ids[i]
      # Remove version number if present
      clean_id <- sub("\\.[0-9]+$", "", transcript_id)
      
      if (clean_id %in% names(transcript_to_gene)) {
        gene_symbols[i] <- transcript_to_gene[[clean_id]]
      }
    }
    
    cat("Mapped", sum(!is.na(gene_symbols)), "out of", length(transcript_ids), "transcripts to gene symbols\n")
  } else {
    cat("GTF file not found, using transcript IDs as gene symbols\n")
  }
  
  # If no symbol found, use the transcript ID without version
  gene_symbols[is.na(gene_symbols)] <- sub("\\.[0-9]+$", "", transcript_ids[is.na(gene_symbols)])
  
  return(gene_symbols)
}

# Create output table with gene symbols
res_df <- data.frame(
  transcript_id = rownames(res),
  gene_symbol = extract_gene_symbols(rownames(res)),
  log2FoldChange = res$log2FoldChange,
  padj = res$padj,
  stringsAsFactors = FALSE
)

# Remove rows with NA padj values and sort by padj
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]

# Write DEGs table
write.table(res_df, file=opt$`output-degs`, sep="\t", quote=FALSE, row.names=FALSE)
cat("Differentially expressed genes written to:", opt$`output-degs`, "\n")

# Parse plot output filenames
plot_files <- strsplit(opt$`output-plots`, ",")[[1]]
if (length(plot_files) != 3) {
  stop("Error: Exactly 3 plot files must be specified (PCA, heatmap, volcano)")
}

pca_file <- plot_files[1]
heatmap_file <- plot_files[2]
volcano_file <- plot_files[3]

# Generate plots
cat("Generating plots...\n")

# Helper function for row variances
rowVars <- function(x) {
  apply(x, 1, var)
}

# Variance stabilizing transformation
suppressWarnings({
  suppressMessages({
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  })
})

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

png(pca_file, width = 800, height = 600)
plot(pca_data$PC1, pca_data$PC2,
     main = "PCA: Cancer vs Normal Tissue",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"),
     col = ifelse(pca_data$condition == "cancer", "red", "blue"),
     pch = 16, cex = 2)
text(pca_data$PC1, pca_data$PC2, labels = pca_data$name, pos = 3, cex = 0.8)
legend("topright", legend = c("Cancer", "Normal"), col = c("red", "blue"), pch = 16)
dev.off()

# Heatmap of top 20 genes
# Get top genes that exist in both results and vsd
cat("Checking genes for heatmap...\n")
cat("Genes in res_df:", nrow(res_df), "\n")
cat("Genes in vsd:", nrow(vsd), "\n")

available_genes <- intersect(rownames(res_df), rownames(vsd))
cat("Overlapping genes:", length(available_genes), "\n")

if (length(available_genes) == 0) {
  # Try using the original gene names from vsd directly
  cat("No overlap found, using top genes from vsd directly...\n")
  vsd_gene_order <- order(rowVars(assay(vsd)), decreasing = TRUE)
  top_genes <- rownames(vsd)[head(vsd_gene_order, min(20, nrow(vsd)))]
} else {
  top_genes <- head(available_genes, min(20, length(available_genes)))
}

cat("Selected", length(top_genes), "genes for heatmap\n")

if (length(top_genes) > 0) {
  heatmap_data <- assay(vsd)[top_genes, ]
  
  png(heatmap_file, width = 800, height = 600)
  heatmap(heatmap_data, 
          main = paste("Top", length(top_genes), "Variable Genes"),
          scale = "row",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          margins = c(8, 8))
  dev.off()
  cat("Heatmap saved to:", heatmap_file, "\n")
} else {
  cat("Error: No genes available for heatmap\n")
}

# Volcano plot
log2FC <- res$log2FoldChange
negLog10P <- -log10(res$padj)
valid_idx <- !is.na(log2FC) & !is.na(negLog10P)

png(volcano_file, width = 800, height = 600)
plot(log2FC[valid_idx], negLog10P[valid_idx],
     main = "Volcano Plot: Cancer vs Normal",
     xlab = "Log2 Fold Change",
     ylab = "-Log10 Adjusted P-value",
     pch = 16, col = "grey", cex = 0.8)

# Highlight significant genes
sig_up <- valid_idx & !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1
sig_down <- valid_idx & !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1

points(log2FC[sig_up], negLog10P[sig_up], col = "red", pch = 16)
points(log2FC[sig_down], negLog10P[sig_down], col = "blue", pch = 16)

# Add significance lines
abline(h = -log10(0.05), lty = 2, col = "black")
abline(v = c(-1, 1), lty = 2, col = "black")

legend("topright", legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "grey"), pch = 16)
dev.off()

# Save R workspace
save(dds, res, vsd, res_df, file = opt$`output-image`)

cat("Analysis complete. Plots saved as:", paste(plot_files, collapse = ", "), "\n")
cat("R workspace saved as:", opt$`output-image`, "\n")