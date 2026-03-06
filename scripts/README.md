## Scripts

Run in the following order after completing the main pipeline steps:

1. count_alleles.py --> QC — independently verify allele counts from burden matrix
   
2. double_check.py -->QC — validate burden matrix against original VCF 
   
3. cosmic_check.py --> Filter significant genes against COSMIC Cancer Gene Census 
   
4. get_nb_genes.py --> Filter cancer genes to neuroblastoma-specific COSMIC hits 
