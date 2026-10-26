#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(gprofiler2)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("--input", "-i"), type="character", default="degs.tsv",
              help="Input file with differentially expressed genes [default %default]"),
  make_option(c("--output", "-o"), type="character", default="pathway_results.tsv",
              help="Output file for pathway enrichment results [default %default]"),
  make_option(c("--padj-threshold", "-p"), type="numeric", default=0.05,
              help="Adjusted p-value threshold for significance [default %default]"),
  make_option(c("--fc-threshold", "-f"), type="numeric", default=1,
              help="Log2 fold change threshold [default %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list,
                          description="Perform pathway enrichment analysis on DEGs")
opt <- parse_args(opt_parser)

# Read input file
cat("Reading differential expression results...\n")
degs <- read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Filter significant genes
sig_genes <- degs[!is.na(degs$padj) & degs$padj < opt$`padj-threshold` & 
                  abs(degs$log2FoldChange) > opt$`fc-threshold`, ]

if (nrow(sig_genes) == 0) {
  cat("No significant genes found with current thresholds\n")
  quit(status=1)
}

cat("Found", nrow(sig_genes), "significant genes for pathway analysis\n")

# Separate upregulated and downregulated genes
up_genes <- sig_genes[sig_genes$log2FoldChange > 0, "gene_symbol"]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, "gene_symbol"]

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")

# Function to run pathway enrichment
run_enrichment <- function(gene_list, direction) {
  if (length(gene_list) < 3) {
    cat("Skipping", direction, "genes: too few genes for enrichment\n")
    return(NULL)
  }
  
  cat("Running pathway enrichment for", direction, "genes...\n")
  
  tryCatch({
    result <- gost(
      query = gene_list,
      organism = "hsapiens",
      sources = c("GO:BP", "KEGG", "REAC"),
      correction_method = "fdr",
      significant = TRUE,
      user_threshold = 0.05
    )
    
    if (!is.null(result) && nrow(result$result) > 0) {
      result$result$direction <- direction
      return(result$result)
    } else {
      cat("No significant pathways found for", direction, "genes\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in pathway enrichment for", direction, "genes:", e$message, "\n")
    return(NULL)
  })
}

# Run enrichment for upregulated and downregulated genes
up_pathways <- run_enrichment(up_genes, "upregulated")
down_pathways <- run_enrichment(down_genes, "downregulated")

# Combine results
all_pathways <- NULL
if (!is.null(up_pathways)) {
  all_pathways <- rbind(all_pathways, up_pathways)
}
if (!is.null(down_pathways)) {
  all_pathways <- rbind(all_pathways, down_pathways)
}

if (!is.null(all_pathways)) {
  # Select relevant columns
  output_columns <- c("query", "source", "term_id", "term_name", "term_size", 
                     "query_size", "intersection_size", "p_value", "direction")
  
  # Filter columns that exist
  available_columns <- intersect(output_columns, colnames(all_pathways))
  pathway_output <- all_pathways[, available_columns]
  
  # Sort by p-value
  pathway_output <- pathway_output[order(pathway_output$p_value), ]
  
  # Write results
  write.table(pathway_output, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Pathway enrichment results written to:", opt$output, "\n")
  cat("Total significant pathways found:", nrow(pathway_output), "\n")
  
  # Display top results
  cat("\nTop 10 enriched pathways:\n")
  print(head(pathway_output[, c("term_name", "p_value", "direction")], 10))
  
} else {
  cat("No significant pathways found\n")
  # Create empty output file
  write.table(data.frame(), file=opt$output, sep="\t", quote=FALSE, row.names=FALSE)
}