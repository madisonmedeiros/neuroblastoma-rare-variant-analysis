"""
count_alleles.py

Independently calculates observed allele counts per gene directly from the
RV-Excalibur burden matrix as a QC validation step. Sums rare variant allele
counts across all patients for each gene.

Note: Minor discrepancies between this output and RV-Excalibur summary results
are expected — RV-Excalibur applies additional coverage and M-CAP filters that
this script does not replicate. This script is for QC purposes only.

File paths below are specific to the UMass Lowell UNITY HPC cluster.
Update these paths to match your local environment before running.

Input:
    - final_output_RVBurdenMatrix_0.9_0_nfe_rcc.txt : RV-Excalibur burden matrix
      (rows = genes, columns = patient IDs, values = rare variant allele counts)

Output:
    - burden_matrix_allele_counts.txt : Per-gene allele count totals across all patients

Usage:
    python3 count_alleles.py
    (Designed to be submitted via sbatch on UMass Lowell UNITY HPC)

Author: Madison Medeiros
"""

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
