##############################
# Obtaining latest immune GWAS
##############################

python preprocessing/upgrade_gwas/get_immune_gwas.sh
python ../bin/gwas-download/munge/custom_munge.py preprocessing/upgrade_gwas/munge_immune_menu.config
python ../bin/gwas-download/overlap/list_snps_to_test.py preprocessing/upgrade_gwas/immune.overlap.config

##############################
# Data prep
##############################

# Add ref and alt variants to each edQTL file
python preprocessing/get_eqtl_ref_alt.py

# Sort and tabix the edQTL files
bash preprocessing/tabix_eqtls.sh

# Get list of SNPs to test for colocalization
python ../bin/gwas-download/overlap/list_snps_to_test.py colocalization/gwas-overlap-clusters.config 20
cat <(cat ../output/complex/test-snps/rna-editing_all-gwas_gtex-aggro_gwas-pval1e-06_eqtl-pval1e-06_gwas-window500000_eqtl-window0_coloc-tests.txt) <(tail -n +2 ../output/complex/test-snps/rna-editing_all-gwas_gtex-single-snp_gwas-pval1e-06_eqtl-pval1e-06_gwas-window500000_eqtl-window0_coloc-tests.txt) > ../output/complex/test-snps/rna-editing_full-list.txt

# Run colocalization tests and assemble them
python colocalization/run_all_tests.py
python colocalization/compilation/concatenate_results.py
# Decide whether using the following filter:
# colocalization/compilation/threshold_by_snp_count.sh

##############################
# General post-analysis
##############################

# Figure 3A
Rscript post-coloc/check_qtl_overlaps/compare_qtl_types.R

# Figure 3B
Rscript post-coloc/interesting_ibd_loci/stack_ibd_loci.R

# Locus-Compare plots for loci in figures 3C and 4B
bash colocalization/custom_coloc_tests/get_ibd_splice_variants.sh > ../output/colocalization/custom_coloc_tests/ibd_splice_variants.txt
python ../bin/coloc-pipeline/dispatch.py colocalization/custom_coloc_tests/figure-3c/interesting-ibd-loci-abridged-3c.config 1
python ../bin/coloc-pipeline/dispatch.py colocalization/custom_coloc_tests/figure-4b/interesting-ibd-loci-abridged-4b.config 1
Rscript colocalization/custom_coloc_tests/plot-3c-4b.R

# Supplementary Tables
Rscript manuscript/supp-tables/create_supp_data.sh	# Compile COLOC results for supplementary table
Rscript post-coloc/interesting_ibd_loci/prioritize_ibd_loci.R  # Used for IBD Supplementary Data Tables


##############################
# Additional analysis not included in paper
##############################
Rscript check_qtl_overlaps/compare_qtl_types/compare_qtl_types.R

# EDA for Huang et al IBD loci
Rscript post-coloc/interesting_ibd_loci/annotate_huang_ibd_loci.R  # Overlapping with existing annotations

# pi1 replication
python post-coloc/check_qtl_p1_overlaps/get_top_sig_edqtls.py
python post-coloc/check_qtl_p1_overlaps/check_qtl_overlaps.py
python post-coloc/check_qtl_p1_overlaps/perform_pi1_tests.R

# Direction-of-effect assessment
python post-coloc/direction_of_effect/edqtl_gwas_direction_comparison.py
