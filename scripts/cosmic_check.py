import pandas as pd


sig_genes = pd.read_csv("significant_genes.txt", sep=r"\s+") 
sig_genes_list = sig_genes['Gene'].tolist() # Add genes to list


cgc = pd.read_csv("cancer_gene_census.csv")
cancer_genes = cgc['Gene Symbol'].dropna().unique().tolist() # Add Cosmic genes to list

# Filter significant genes that are in CGC
cancer_related = sig_genes[sig_genes['Gene'].isin(cancer_genes)]

# Save to new file
cancer_related.to_csv("significant_cancer_genes.txt", sep="\t", index=False)