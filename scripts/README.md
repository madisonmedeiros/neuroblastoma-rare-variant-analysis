## Scripts

Run in the following order after completing the main pipeline steps:

count_alleles.py --> QC — independently verify allele counts from burden matrix |
double_check.py -->QC — validate burden matrix against original VCF |
cosmic_check.py --> Filter significant genes against COSMIC Cancer Gene Census |
get_nb_genes.py --> Filter cancer genes to neuroblastoma-specific COSMIC hits |
