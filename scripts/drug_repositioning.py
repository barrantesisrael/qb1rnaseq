#!/usr/bin/env python3

import argparse
import requests
import json
import sys
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description='Drug repositioning analysis using L1000CDS2')
    parser.add_argument('--input', '-i', required=True,
                       help='Input file with differentially expressed genes (TSV format)')
    parser.add_argument('--output', '-o', default='drug_candidates.txt',
                       help='Output file for drug candidates (default: drug_candidates.txt)')
    parser.add_argument('--padj-threshold', '-p', type=float, default=0.05,
                       help='Adjusted p-value threshold (default: 0.05)')
    parser.add_argument('--fc-threshold', '-f', type=float, default=1.0,
                       help='Log2 fold change threshold (default: 1.0)')
    parser.add_argument('--max-genes', '-m', type=int, default=150,
                       help='Maximum number of genes per direction (default: 150)')
    return parser.parse_args()

def read_degs(input_file, padj_threshold, fc_threshold):
    """Read and filter differentially expressed genes"""
    try:
        degs = []
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Convert numeric columns
                try:
                    padj = float(row['padj']) if row['padj'] != 'NA' else None
                    log2fc = float(row['log2FoldChange']) if row['log2FoldChange'] != 'NA' else None
                    
                    if padj is not None and log2fc is not None:
                        row['padj'] = padj
                        row['log2FoldChange'] = log2fc
                        degs.append(row)
                except ValueError:
                    continue
        
        print(f"Read {len(degs)} genes from {input_file}")
        
        # Filter significant genes
        sig_degs = [
            gene for gene in degs
            if gene['padj'] < padj_threshold and abs(gene['log2FoldChange']) > fc_threshold
        ]
        
        print(f"Found {len(sig_degs)} significant genes")
        return sig_degs
        
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

def prepare_l1000cds2_query(sig_degs, max_genes):
    """Prepare gene lists for L1000CDS2 query"""
    # Separate upregulated and downregulated genes
    up_genes = [gene['gene_symbol'].upper() for gene in sig_degs if gene['log2FoldChange'] > 0]
    down_genes = [gene['gene_symbol'].upper() for gene in sig_degs if gene['log2FoldChange'] < 0]
    
    # Limit number of genes
    up_genes = up_genes[:max_genes]
    down_genes = down_genes[:max_genes]
    
    print(f"Using {len(up_genes)} upregulated and {len(down_genes)} downregulated genes")
    
    # Create L1000CDS2 payload
    payload = {
        "data": {
            "upGenes": up_genes,
            "dnGenes": down_genes
        },
        "config": {
            "aggravate": False,  # Reverse mode for therapeutic applications
            "searchMethod": "geneSet",
            "share": True,
            "combination": True,
            "db-version": "latest"
        }
    }
    
    return payload

def query_l1000cds2(payload):
    """Query L1000CDS2 API for drug predictions"""
    url = "https://maayanlab.cloud/L1000CDS2/query"
    
    try:
        print("Querying L1000CDS2 database...")
        response = requests.post(url, json=payload, timeout=60)
        
        if response.status_code == 200:
            result = response.json()
            if 'topMeta' in result and result['topMeta']:
                print(f"Found {len(result['topMeta'])} drug candidates")
                return result['topMeta']
            else:
                print("No drug candidates returned from L1000CDS2")
                return []
        else:
            print(f"L1000CDS2 API error: HTTP {response.status_code}")
            return []
            
    except Exception as e:
        print(f"Error querying L1000CDS2: {e}")
        return []

def write_drug_results(drug_results, output_file):
    """Write drug results to TSV file"""
    if not drug_results:
        # Create empty file
        with open(output_file, 'w') as f:
            f.write("pert_desc\tpert_id\tcell_id\tpert_dose\tpert_dose_unit\tpert_time\tpert_time_unit\tscore\n")
        return
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("pert_desc\tpert_id\tcell_id\tpert_dose\tpert_dose_unit\tpert_time\tpert_time_unit\tscore\n")
        
        # Write data
        for drug in drug_results:
            pert_desc = drug.get('pert_desc', 'Unknown')
            pert_id = drug.get('pert_id', 'Unknown')
            cell_id = drug.get('cell_id', 'Unknown')
            pert_dose = drug.get('pert_dose', 'Unknown')
            pert_dose_unit = drug.get('pert_dose_unit', 'Unknown')
            pert_time = drug.get('pert_time', 'Unknown')
            pert_time_unit = drug.get('pert_time_unit', 'Unknown')
            score = drug.get('score', 0.0)
            
            f.write(f"{pert_desc}\t{pert_id}\t{cell_id}\t{pert_dose}\t{pert_dose_unit}\t{pert_time}\t{pert_time_unit}\t{score}\n")

def main():
    args = parse_arguments()
    
    # Read and filter DEGs
    sig_degs = read_degs(args.input, args.padj_threshold, args.fc_threshold)
    
    if len(sig_degs) == 0:
        print("No significant genes found with current thresholds")
        # Create empty output file
        write_drug_results([], args.output)
        return
    
    # Prepare L1000CDS2 query
    payload = prepare_l1000cds2_query(sig_degs, args.max_genes)
    
    # Query L1000CDS2
    drug_results = query_l1000cds2(payload)
    
    # Write results
    write_drug_results(drug_results, args.output)
    
    if len(drug_results) > 0:
        print(f"Drug candidates written to: {args.output}")
        print(f"Top 5 candidates:")
        for i, drug in enumerate(drug_results[:5]):
            print(f"{i+1}. {drug.get('pert_desc', 'Unknown')} (score: {drug.get('score', 0.0)})")
    else:
        print("No drug candidates found")

if __name__ == "__main__":
    main()