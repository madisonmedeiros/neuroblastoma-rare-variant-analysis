"""
cosmic_check.py

Cross-references rare variant burden testing results against the COSMIC
Cancer Gene Census (CGC) to identify which significant genes have a known
role in cancer.

Inputs:
    - significant_genes.txt   : RV-Excalibur output filtered to P < 0.05 (75 genes)
    - cancer_gene_census.csv  : COSMIC CGC downloaded from https://cancer.sanger.ac.uk/cosmic

Output:
    - significant_cancer_genes.txt : Significant genes overlapping with COSMIC CGC

Usage:
    python3 cosmic_check.py

Author: Madison Medeiros
"""

import pandas as pd

sig_genes = pd.read_csv("significant_genes.txt", sep=r"\s+") 
sig_genes_list = sig_genes['Gene'].tolist() # Add genes to list


cgc = pd.read_csv("cancer_gene_census.csv")
cancer_genes = cgc['Gene Symbol'].dropna().unique().tolist() # Add Cosmic genes to list

# Filter significant genes that are in CGC
cancer_related = sig_genes[sig_genes['Gene'].isin(cancer_genes)]

# Save to new file
cancer_related.to_csv("significant_cancer_genes.txt", sep="\t", index=False)
