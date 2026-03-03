# Neuroblastoma Rare Variant Burden Analysis

A bioinformatics pipeline for identifying rare genetic variants associated with neuroblastoma using gene-based burden testing. This project was conducted as graduate research at the **University of Massachusetts Lowell** under PI Rachel Melamed.

---

## Background

Neuroblastoma is the most common solid extracranial malignancy of childhood, accounting for ~7% of all cancers in children under 15 and representing the most common cancer in the first year of life. Despite known drivers such as MYCN amplification and ALK mutations, the full genetic landscape remains incompletely characterized.

This pipeline uses **burden testing** — a rare variant aggregation method — to identify genes enriched for rare, potentially pathogenic variants in neuroblastoma cases compared to population-level expectations from gnomAD. Unlike traditional GWAS, which focuses on common variants, burden testing aggregates all rare variants within a gene to detect cumulative effects that individually would not reach significance.

---

## Pipeline Overview

```
Raw VCF (cohort) — genome build hg38
      │
      ▼
[1] VCF Cleaning & QC
      │  - Remove messy headers
      │  - Handle missing genotypes (./.:.:.:.:.:.:.:.:..:.)
      │  - Filter to PASS variants
      ▼
[2] Liftover hg38 → hg19
      │  - Convert coordinates to hg19 (GRCh37)
      │  - Required for ANNOVAR refGene + gnomAD 2.1.1 compatibility
      ▼
[3] ANNOVAR Annotation
      │  - refGene (gene region coordinates)
      │  - gnomAD 2.1.1 (population allele frequencies)
      │  - dbNSFP / M-CAP (pathogenicity scores)
      ▼
[4] PLINK Conversion
      │  - Convert annotated VCF → .bed / .bim / .fam
      ▼
[5] gnomAD Data Download
      │  - Population-level expected allele counts (EAC)
      │  - Coverage harmonization (≥20x in ≥90% of gnomAD samples)
      ▼
[6] RV-Excalibur Burden Testing
      │  - Generate rare variant burden matrix
      │  - Apply MAF × M-CAP filtering thresholds
      │  - Calculate iCF and gCF correction factors
      │  - Generate per-gene Z-scores and P-values
      ▼
[7] Results Filtering & Validation
      │  - Filter significant genes (P < 0.05)
      │  - Cross-reference against known neuroblastoma genes
      │  - Validate via COSMIC, gnomAD, literature
      ▼
Candidate Gene List
```

---

## Repository Structure

```
neuroblastoma-rare-variant-analysis/
├── scripts/
│   ├── count_alleles.py                          # Independent allele count verification script
│   ├── double_check.py                           # Burden matrix validation against VCF data
│   ├── cosmic_check.py                           # Filter significant genes against COSMIC Cancer Gene Census
│   └── get_nb_genes.py                           # Filter cancer genes down to neuroblastoma-specific hits
├── data/
│   ├── cancer_gene_census.csv                    # COSMIC Cancer Gene Census (input for filtering)
│   └── nb_genes.txt                              # Known neuroblastoma-associated genes reference list
├── results/
│   ├── rvexcaliber_summary_results.txt           # Raw RV-Excalibur output (5,953 genes tested)
│   ├── significant_genes.txt                     # P < 0.05 filtered results (75 genes)
│   ├── significant_cancer_genes.txt              # COSMIC CGC overlap (97 genes)
│   ├── new_significant_cancer_genes.txt          # Top hits by enrichment (6 genes)
│   ├── nb_related_genes_cosmic.txt               # COSMIC-confirmed NB genes with stats
│   └── nb_genes.txt                              # Final neuroblastoma gene list (ALK)
└── README.md
```

---

## Requirements

### Software
| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 1.9+ | Genotype file conversion |
| ANNOVAR | latest | Variant annotation |
| RV-Excalibur | — | Gene-based burden testing |
| BCFtools / VCFtools | 1.15+ | VCF merging and filtering |
| BEDTools | 2.31.1+ | Genomic interval operations |
| Python | 3.12+ | Custom scripts |
| pandas | — | Data manipulation |

### Python Packages
```bash
pip install pandas
```

### HPC Environment
This pipeline was developed and run on the **UMass Lowell UNITY HPC cluster** using SLURM for job scheduling. Scripts are designed to be submitted via `sbatch`.

---

## Usage

### Step 0: Download ANNOVAR Databases
```bash
bash scripts/0_get_ANNOVAR_annotations.clean.sh
```
Downloads `refGene`, `gnomAD 2.1.1`, and `dbNSFP` to ANNOVAR's `/humandb` directory. Only needs to be run once.

---

### Step 1: Clean VCF

```bash
# Keep only PASS variants and clean headers
grep -E "^#|PASS" merged_cohort.vcf > merged_pass.vcf
```

---

### Step 2: Liftover hg38 → hg19

The raw neuroblastoma VCF data was originally in **hg38 (GRCh38)** build. Since ANNOVAR's `refGene` database and gnomAD 2.1.1 are both based on **hg19 (GRCh37)**, a liftover was required before annotation.

```bash
# Download the hg38 to hg19 chain file (if not already present)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz

# Run liftover using GATK
gatk LiftoverVcf \
  -I merged_pass.vcf \
  -O output_hg19.vcf \
  --CHAIN hg38ToHg19.over.chain \
  --REJECT rejected_variants.vcf \
  -R hg19_reference.fa

# Alternatively, using CrossMap
CrossMap.py vcf hg38ToHg19.over.chain merged_pass.vcf hg19_reference.fa output_hg19.vcf
```

> **Note:** Some variants may fail liftover and will be written to `rejected_variants.vcf`. Review these to assess how many variants were lost in conversion. The lifted-over file `output_hg19.vcf` is used for all downstream steps.

---

### Step 3: Annotate with ANNOVAR

```bash
# Add ANNOVAR to PATH
export PATH=$PATH:/path/to/annovar

# Annotate with ANNOVAR (hg19)
table_annovar.pl output_hg19.vcf humandb/ \
  -buildver hg19 \
  -out annotated_output \
  -remove \
  -protocol refGene,gnomad211_genome,dbnsfp42a \
  -operation g,f,f \
  -nastring .
```

---

### Step 4: Convert to PLINK Binary Format

```bash
plink --vcf annotated_output.vcf \
      --make-bed \
      --out neuroblastoma_plink

# Copy PLINK files to analysis directory
cp neuroblastoma_plink.bed \
   neuroblastoma_plink.bim \
   neuroblastoma_plink.fam \
   /path/to/analysis_files/plink_files/
```

---

### Step 5: Generate Rare Variant Burden Matrix

```bash
bash scripts/1A_get_RVBurdenMatrix_internal.clean.sh \
  [vcf_file] \
  [plink_prefix] \
  [gnomad_maf_threshold] \
  [mcap_threshold] \
  [gene_list] \
  [coverage_file] \
  [output_prefix] \
  [population] \
  [n_cores]
```

Output: `final_output_RVBurdenMatrix_[MAF]_[MCAP]_[pop]_[suffix].txt`
- Rows = genes
- Columns = patient IDs
- Values = rare variant allele counts per patient per gene

---

### Step 6: Run Association Analysis

```bash
bash scripts/3_get_SummaryAssociations_iCF_gCF.clean.sh \
  [burden_matrix] \
  [gnomad_eac_file] \
  [ranking_dataset] \
  [output_prefix] \
  [maf_threshold] \
  [mcap_threshold] \
  [population] \
  [n_cores]
```

This script:
1. Calculates **observed allele counts (OAC)** from the internal dataset
2. Calculates **expected allele counts (EAC)** from gnomAD population frequencies
3. Computes **iCF** (individual correction factor): ratio of total OAC to total EAC per individual
4. Computes **gCF** (gene correction factor): per-gene-bin ratio of OAC to iCF-adjusted EAC
5. Outputs per-gene Z-scores and P-values
6. Generates a **mountain plot** (cumulative delta allele count plot)

**iCF Formula:**

$$\widehat{iCF}_i = \frac{\sum_{g=1}^{M} \text{OAC}_{i,g}}{\sum_{g=1}^{M} \text{EAC}_{i,g}}$$

---

### Step 7: Verify Allele Counts (Optional QC)

```bash
python3 scripts/count_alleles.py
```

This script independently calculates allele counts per gene directly from the VCF and ANNOVAR output to validate RV-Excalibur results.

> **Note:** Minor discrepancies between this script and RV-Excalibur are expected because RV-Excalibur applies additional quality filters (coverage thresholds, M-CAP filtering). This script is for QC purposes only.

---

### Step 8: Filter Results & Identify Candidate Genes

Results were filtered in three stages:

**Stage 1 — Filter by P-value:** The raw RV-Excalibur summary file (`rvexcaliber_summary_results.txt`) contains results for all 5,953 tested genes. Significant genes (P < 0.05) were extracted with awk:

```bash
awk 'NR==1 || $NF < 0.05' rvexcaliber_summary_results.txt > significant_genes.txt
# Produces 75 significant genes
```

**Stage 2 — `cosmic_check.py`:** Cross-references the 75 significant genes against the COSMIC Cancer Gene Census (CGC), narrowing to genes with a known role in cancer:

```bash
python3 scripts/cosmic_check.py
# Input:  significant_genes.txt + cancer_gene_census.csv
# Output: significant_cancer_genes.txt (97 genes)
```

**Stage 3 — `get_nb_genes.py`:** Filters the cancer gene list further to only genes where COSMIC specifically lists "neuroblastoma" in the somatic tumour types column:

```bash
python3 scripts/get_nb_genes.py
# Input:  significant_cancer_genes.txt + cancer_gene_census.csv
# Output: nb_related_genes_cosmic.txt, nb_genes.txt
```

---

## Results

The pipeline narrowed findings through three progressive filters:

| Stage | Gene Count | Description |
|-------|-----------|-------------|
| RV-Excalibur output | 5,953 genes | Full burden test results with iCF/gCF-adjusted Z-scores |
| `significant_genes.txt` | 75 genes | P < 0.05 after iCF/gCF correction |
| `significant_cancer_genes.txt` | 97 genes | Overlap with COSMIC Cancer Gene Census |
| `new_significant_cancer_genes.txt` | 6 genes | Top hits by significance (P < 0.001 or strong enrichment) |
| `nb_related_genes_cosmic.txt` | 1 gene | COSMIC-confirmed neuroblastoma gene |

---

### COSMIC-Confirmed Neuroblastoma Hit

| Gene | Observed Allele Count | gnomAD Expected (iCF/gCF adjusted) | P-value |
|------|--------------------|--------------------------------------|---------|
| ALK | 4 | 3.85 | 0.469 |

**ALK** (Anaplastic Lymphoma Kinase) is one of the most well-established neuroblastoma driver genes, with somatic mutations in ~8% of sporadic cases and germline mutations in familial neuroblastoma. Its recovery here as the sole COSMIC-confirmed neuroblastoma gene serves as a meaningful validation that the pipeline is correctly detecting biologically relevant signal.

---

### Top Significant Cancer Genes & Literature Findings

The six most significantly enriched genes from `new_significant_cancer_genes.txt` were each investigated through manual PubMed literature review. Two of the six — **PDE4DIP** and **CSMD3** — have direct neuroblastoma associations in the literature despite not being flagged by COSMIC, representing the novel findings of this analysis.

| Gene | Observed | Expected | P-value | Neuroblastoma Relevance |
|------|----------|----------|---------|------------------------|
| MUC16 | 1993 | 352.79 | 0.0 | No direct NB link; promotes metastasis via HuR/cMyc axis in TNBC |
| MUC6 | 2871 | 313.02 | 0.0 | Tumor suppressor via autophagy/β-catenin; found mutated in Wilms tumor alongside PDE4DIP |
| PDE4DIP | 27 | 7.37 | 6.7 × 10⁻⁵ | ✅ **Direct NB link** — identified in a MYCN-related neuroblastoma prognostic model (AKR1C1, CHD5, PDE4DIP, PRKACB) predicting patient survival across 900 NB samples |
| CSMD3 | 26 | 6.35 | 1.5 × 10⁻⁴ | ✅ **Direct NB link** — increased copy number and expression found in metastatic site-derived aggressive neuroblastoma cells; associated with poor overall and relapse-free survival |
| LRP1B | 29 | 13.40 | 0.0024 | Tumor suppressor; suppresses β-catenin/TCF signaling in colon cancer; loss-of-function variants associated with improved immunotherapy response |
| TPR | 9 | 3.40 | 0.030 | Frameshift mutations found in MSI-H colon cancers; role in NB not established |

> **Key finding:** PDE4DIP and CSMD3 are enriched with rare variants in this neuroblastoma cohort AND have independent literature support for involvement in neuroblastoma — making them the strongest novel candidate genes from this analysis. These connections were not captured by COSMIC, highlighting the value of combining computational burden testing with manual literature validation.




---

## Key Parameters

| Parameter | Value Used | Description |
|-----------|-----------|-------------|
| gnomAD MAF threshold | 0.01 (1%) | Rare variant cutoff |
| M-CAP threshold | 0 (all functional) | Pathogenicity score filter |
| gnomAD population | NFE (non-Finnish European) | Reference population for EAC |
| Coverage cutoff | ≥20x in ≥90% of samples | Coverage harmonization |
| Genome build | hg19 (via liftover from hg38) | Reference genome |

---

## Validation Resources

After identifying candidate genes, the following databases were used for follow-up:

- **[COSMIC](https://cancer.sanger.ac.uk/cosmic)** — Filter: Nervous System > Neuroblastoma (not CNS)
- **[gnomAD](https://gnomad.broadinstitute.org)** — Verify population allele frequencies per variant
- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — Check germline variant classifications
- **[PubMed](https://pubmed.ncbi.nlm.nih.gov)** — Literature search: `[GENE] AND neuroblastoma`

---

## File Path Reference (UML UNITY HPC)

```
/project/pi_rachel_melamed_uml_edu/neuroblastoma/maddie/
├── annovar/                    # ANNOVAR installation
├── double_checking/            # Allele count verification scripts
│   ├── output_hg19.vcf
│   ├── hg19_refGene.txt
│   └── final_output_RVBurdenMatrix_0.9_0_nfe_rcc.txt
└── final_steps/
    └── analysis_files/
        └── plink_files/        # .bed / .bim / .fam files
```

---

## Author

**Madison Medeiros**  
M.S. Bioinformatics, University of Massachusetts Lowell  
[medeirosmm@merrimack.edu](mailto:medeirosmm@merrimack.edu)

---

## Acknowledgments

- **PI Rachel Melamed** — University of Massachusetts Lowell
- **RV-Excalibur** — Rare variant burden testing framework
- **gnomAD** — Genome Aggregation Database (Broad Institute)
- **ANNOVAR** — Wang Lab, University of Southern California
