#!/usr/bin/env Rscript

# Load results
load("deseq2_results.RData")
library(gprofiler2)

# Get gene lists
oncogenes_list <- rownames(oncogenes)
tumor_suppressors_list <- rownames(tumor_suppressors)

# Clean gene names (remove version numbers)
clean_genes <- function(genes) {
    gsub("\\.[0-9]+$", "", genes)
}

oncogenes_clean <- clean_genes(oncogenes_list)
suppressors_clean <- clean_genes(tumor_suppressors_list)

cat("Oncogenes for analysis:", length(oncogenes_clean), "\n")
cat("Tumor suppressors for analysis:", length(suppressors_clean), "\n")

# Save gene lists for drug repositioning
write.table(oncogenes_list, "oncogenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(tumor_suppressors_list, "tumor_suppressors.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Pathway analysis
if(length(oncogenes_clean) > 3) {
    cat("Analyzing oncogene pathways...\n")
    tryCatch({
        oncogene_pathways <- gost(
            query = oncogenes_clean,
            organism = "hsapiens",
            sources = c("GO:BP", "KEGG")
        )
        
        if(!is.null(oncogene_pathways)) {
            cat("Top oncogene pathways:\n")
            print(head(oncogene_pathways$result[,c("term_name", "p_value")], 5))
        }
    }, error = function(e) {
        cat("Error in oncogene pathway analysis:", e$message, "\n")
    })
}

if(length(suppressors_clean) > 3) {
    cat("Analyzing tumor suppressor pathways...\n")
    tryCatch({
        suppressor_pathways <- gost(
            query = suppressors_clean,
            organism = "hsapiens", 
            sources = c("GO:BP", "KEGG")
        )
        
        if(!is.null(suppressor_pathways)) {
            cat("Top tumor suppressor pathways:\n")
            print(head(suppressor_pathways$result[,c("term_name", "p_value")], 5))
        }
    }, error = function(e) {
        cat("Error in tumor suppressor pathway analysis:", e$message, "\n")
    })
}
