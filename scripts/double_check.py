import pandas as pd

print("Loading in files", flush=True)

# File Paths
base_path = "/project/pi_rachel_melamed_uml_edu/neuroblastoma/maddie/double_checking/"
vcf_file = base_path + "output_hg19.vcf"
burden_matrix_file = base_path + "final_output_RVBurdenMatrix_0.9_0_nfe_rcc.txt"  
annovar_file = base_path + "hg19_refGene.txt"  

print("Files successfully loaded", flush=True)

# Step 1: Load Gene Regions from ANNOVAR
print("Loading ANNOVAR file", flush=True)
gene_regions = {}

with open(annovar_file, 'r') as annovar:
    for line in annovar:
        parts = line.strip().split("\t")
        if len(parts) < 5:  # Skip malformed lines
            continue
        gene_name = parts[12]  # Gene name
        chrom = parts[2]       # Chromosome
        start = int(parts[4])  # Start position
        stop = int(parts[5])   # Stop position

        if gene_name not in gene_regions:
            gene_regions[gene_name] = []  # Initialize list for regions
        gene_regions[gene_name].append((chrom, start, stop))

print(f"Loaded {len(gene_regions)} genes from ANNOVAR", flush=True)

# Step 2: Load Variants from VCF (ONLY Store Variants)
print("Loading VCF file", flush=True)
patient_variants = {}  # {patient: [(chrom, pos), ...]}
patient_ids = []

with open(vcf_file, 'r') as vcf:
    for line in vcf:
        if line.startswith("#CHROM"):
            header_parts = line.strip().split("\t")
            patient_ids = header_parts[9:]  # Extract patient IDs
            for patient_id in patient_ids:
                patient_variants[patient_id] = []
            continue

        if line.startswith("#"):
            continue  # Skip header lines

        parts = line.strip().split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        samples = parts[9:]

        for idx, genotype in enumerate(samples):
            patient_id = patient_ids[idx]
            genotype_value = genotype.split(":")[0]  # Extract only genotype (GT)

            if "1" in genotype_value or "2" in genotype_value:
                patient_variants[patient_id].append((chrom, pos))  # Store variant location

print(f"Loaded variants for {len(patient_variants)} patients from VCF", flush=True)

# Step 3: Load Burden Matrix
print("Loading Burden Matrix", flush=True)
burden_matrix = pd.read_csv(burden_matrix_file, sep=r"\s+", index_col=0, engine="python")

print(f"Loaded burden matrix with {burden_matrix.shape[0]} genes and {burden_matrix.shape[1]} patients", flush=True)

# Step 4: Normalize Burden Matrix Gene Names
print("Normalizing gene names for comparison...", flush=True)

# Extract first gene name from burden matrix (removing fusion genes)
burden_matrix_genes_cleaned = {gene.split(";")[0]: gene for gene in burden_matrix.index}

# Find missing genes in ANNOVAR
missing_genes = set(burden_matrix_genes_cleaned.keys()) - set(gene_regions.keys())
print(f"Genes in Burden Matrix but NOT in ANNOVAR (first 20): {list(missing_genes)[:20]}", flush=True)

# Find genes in burden matrix that have NO variants in VCF
genes_in_vcf = set()  # Track genes that have at least one variant in the VCF
for patient in patient_variants:
    for gene, regions in gene_regions.items():
        for chrom, start, stop in regions:
            for v_chrom, v_pos in patient_variants[patient]:
                if v_chrom == chrom and start <= v_pos <= stop:
                    genes_in_vcf.add(gene)
                    break
            if gene in genes_in_vcf:
                break  # No need to check further regions for this gene

genes_without_variants = [
    gene for gene in burden_matrix_genes_cleaned if gene in gene_regions and gene not in genes_in_vcf
]
print(f"Genes in Burden Matrix but have NO variants in VCF (first 20): {genes_without_variants[:20]}", flush=True)

print("Filtered down to relevant genes present in both Burden Matrix and VCF.", flush=True)

# Step 5: Validate Burden Matrix Variants Using VCF
errors = []
matches = []

print("Checking for mismatches between Burden Matrix and VCF...", flush=True)

for patient in burden_matrix.columns:
    if patient not in patient_variants:  # If patient is missing in VCF, flag it
        print(f"⚠ Warning: Patient {patient} is in Burden Matrix but missing from VCF!")
        continue

    for gene in burden_matrix.index:
        cleaned_gene = gene.split(";")[0]  # Take the first gene name
        reported_in_matrix = burden_matrix.loc[gene, patient]

        # Only check genes marked as '1' (has a variant in burden matrix)
        if reported_in_matrix == 1:
            has_variant = False

            # Check if the patient has a variant in the gene region
            if cleaned_gene in gene_regions:
                for chrom, start, stop in gene_regions[cleaned_gene]:
                    for v_chrom, v_pos in patient_variants[patient]:
                        if v_chrom == chrom and start <= v_pos <= stop:
                            has_variant = True
                            break
                    if has_variant:
                        break

            # Log results
            if has_variant:
                matches.append(f"Match: {patient} has a variant in {gene} (VCF & Burden Matrix)")
            else:
                errors.append(f"Mismatch: {patient} is marked as 1 in Burden Matrix but has no variant in VCF for {gene}")

print("Variant check completed.", flush=True)

# Step 6: Print Matches & Errors
if matches:
    print(f"✅ Found {len(matches)} correct matches! Showing first 10:")
    for match in matches[:10]:
        print(match)

if errors:
    print(f"❌ Found {len(errors)} mismatches! Showing first 10:")
    for error in errors[:10]:
        print(error)
else:
    print("🎉 No mismatches found! Burden Matrix is consistent with VCF.")

print("Updated variant checking job completed.", flush=True)
