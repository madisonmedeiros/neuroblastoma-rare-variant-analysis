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