"""
get_nb_genes.py

Filters the COSMIC Cancer Gene Census overlap results down to genes where
COSMIC specifically lists "neuroblastoma" in the somatic tumour types column.
This is the final filtering stage in the candidate gene discovery pipeline,
narrowing from cancer-associated genes to neuroblastoma-specific hits.

Pipeline context:
    significant_genes.txt (75 genes, P < 0.05)
        → cosmic_check.py
    significant_cancer_genes.txt (97 genes, COSMIC CGC overlap)
        → get_nb_genes.py  ← this script
    nb_related_genes_cosmic.txt (neuroblastoma-specific hits)
    nb_genes.txt (gene names only, for downstream lookup)

Inputs:
    - significant_cancer_genes.txt : Output of cosmic_check.py
    - cancer_gene_census.csv       : COSMIC CGC downloaded from
                                     https://cancer.sanger.ac.uk/cosmic

Outputs:
    - nb_related_genes_cosmic.txt  : Significant genes where COSMIC lists
                                     neuroblastoma as a somatic tumour type
    - nb_genes.txt                 : Gene names only (one per line),
                                     for use in downstream validation

Usage:
    python3 get_nb_genes.py

Author: Madison Medeiros
"""

import pandas as pd

# Load your filtered list of cancer-related significant genes
sig_cancer_genes = pd.read_csv("significant_cancer_genes.txt", sep=r"\s+")

# Load the Cosmic genes
cgc = pd.read_csv("cancer_gene_census.csv")

# Select columns and clean
cgc_subset = cgc[['Gene Symbol', 'Tumour Types(Somatic)']].dropna()
cgc_subset['Gene Symbol'] = cgc_subset['Gene Symbol'].str.upper()
sig_cancer_genes['Gene'] = sig_cancer_genes['Gene'].str.upper()

# Filter COSMIC rows that mention "neuroblastoma" in the tumor type
cgc_neuro = cgc_subset[cgc_subset['Tumour Types(Somatic)'].str.contains("neuroblastoma", case=False)]

# Get the set of gene symbols associated with neuroblastoma
neuroblastoma_genes = set(cgc_neuro['Gene Symbol'])

# Filter your file for genes matching that list
nb_related_genes = sig_cancer_genes[sig_cancer_genes['Gene'].isin(neuroblastoma_genes)]

# Save the final neuroblastoma-related list
nb_related_genes.to_csv("nb_related_genes_cosmic.txt", sep="\t", index=False)

print(f"Found {len(nb_related_genes)} neuroblastoma-associated genes in your significant cancer gene list.")

# Save only the gene names to nb_genes.txt
nb_related_genes['Gene'].drop_duplicates().to_csv("nb_genes.txt", index=False, header=False)
