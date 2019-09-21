##############################
# Data prep
##############################

# Add ref and alt variants to each edQTL file
python preprocessing/get_eqtl_ref_alt.py

# Sort and tabix the edQTL files
bash preprocessing/tabix_eqtls.sh

# Get list of SNPs to test for colocalization
python /users/mgloud/projects/gwas-download/overlap/list_snps_to_test.py /users/mgloud/projects/rna_editing/scripts/gwas-overlap-clusters.config 20
cat <(cat /users/mgloud/projects/rna_editing/output/test-snps/rna-editing_all-gwas_gtex-aggro_gwas-pval1e-06_eqtl-pval1e-06_gwas-window500000_eqtl-window0_coloc-tests.txt) <(tail -n +2 /users/mgloud/projects/rna_editing/output/test-snps/rna-editing_all-gwas_gtex-single-snp_gwas-pval1e-06_eqtl-pval1e-06_gwas-window500000_eqtl-window0_coloc-tests.txt) > /users/mgloud/projects/rna_editing/output/test-snps/rna-editing_full-list.txt

# Run colocalization tests
python run_all_tests.py

##############################
# Figure 1
##############################


##############################
# Figure 2
##############################

##############################
# Figure 3
##############################
