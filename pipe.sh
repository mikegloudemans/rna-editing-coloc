##############################
# Obtaining latest immune GWAS
##############################

# GWAS were downloaded and munged using the default scripts that
# are downloadable with the gwas-download toolkit

# Other GWAS not included here were munged using config files directly in the gwas-download module
python bin/gwas-download/munge/custom_munge.py scripts/preprocessing/munge_config/munge_all_rna_editing.config
python bin/gwas-download/munge/custom_munge.py scripts/preprocessing/munge_config/munge_MS.config

##############################
# Data prep
##############################

# Add ref and alt variants to each edQTL file
python preprocessing/get_eqtl_ref_alt.py

# Sort and tabix the edQTL files
bash preprocessing/tabix_eqtls.sh

##############################
# Determine GWAS / QTL overlap
##############################

# Get list of SNPs to test for colocalization
# NOTE: First, make sure python2 and tabix are loaded
python bin/gwas-download/overlap/list_snps_to_test.py scripts/preprocessing/overlap_config/rna-editing.overlap.2020-09-15.config 20

##############################
# Run colocalization
##############################

# Run colocalization tests
python bin/coloc-pipeline/dispatch.py scripts/colocalization/config/rna-coloc.config 20

# TODO: concatenate all colocalization tests

# Filter down to sites containing 50 or more tested SNPs
python colocalization/compilation/threshold_by_snp_count.sh

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
