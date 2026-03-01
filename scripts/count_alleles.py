import pandas as pd

print("Loading Burden Matrix...", flush=True)

# File Path to the Burden Matrix
burden_matrix_file = "/project/pi_rachel_melamed_uml_edu/neuroblastoma/maddie/double_checking/final_output_RVBurdenMatrix_0.9_0_nfe_rcc.txt"

# Load Burden Matrix
burden_matrix = pd.read_csv(burden_matrix_file, sep=r"\s+", index_col=0, engine="python")

print(f"Loaded burden matrix with {burden_matrix.shape[0]} genes and {burden_matrix.shape[1]} patients", flush=True)

# Compute allele count per gene (sum across patient columns)
allele_counts = burden_matrix.sum(axis=1)  # Sum across columns (patients)

# Convert to DataFrame
allele_counts_df = allele_counts.reset_index()
allele_counts_df.columns = ["Gene", "Allele_Count"]

# Save to output file
output_file = "/project/pi_rachel_melamed_uml_edu/neuroblastoma/maddie/count_alleles/burden_matrix_allele_counts.txt"
allele_counts_df.to_csv(output_file, sep="\t", index=False)

# Print first 10 results
print("Allele counts computed! First 10 genes:")
print(allele_counts_df.head(10))

print(f"Allele counts saved to: {output_file}")
