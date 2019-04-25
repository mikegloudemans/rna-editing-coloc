# Author: Mike Gloudemans

require(readr)
require(dplyr)

# Load files

# Load file showing all significant GWAS SNPs for a locus
all_gwas_snps = read_delim("/users/mgloud/projects/rna_editing/output/test-snps/rna-editing_all-gwas_gwas-pval1e-06_gwas-window500000_snps-considered.txt", delim = "\t")
all_gwas_snps$id = paste(all_gwas_snps$chr, all_gwas_snps$snp_pos, sep="_")

gwas_hit_counts = all_gwas_snps %>% group_by(trait) %>% summarize(total = length(unique(id)))

# Load file showing all SNPs passing the initial edQTL cutoffs for each trait
all_coloc_tests = read_delim("/users/mgloud/projects/rna_editing/output/test-snps/rna-editing_all-gwas_gtex-single-snp_gwas-pval1e-06_eqtl-pval1e-06_gwas-window500000_eqtl-window0_coloc-tests.txt", delim = "\t")
all_coloc_tests$id = paste(all_coloc_tests$chr, all_coloc_tests$snp_pos, sep="_")

coloc_test_counts = all_coloc_tests %>% group_by(trait) %>% summarize(total = length(unique(id)))

all_coloc_results = read_delim("/users/mgloud/projects/rna_editing/output/single_site_clpp_results.tsv", delim = "\t")
coloc_result_counts = all_coloc_results %>% group_by(gwas_trait) %>% summarize(total = length(unique(ref_snp)))


# Get percentages then...


# See how different or similar this is to the percent overlaps with eQTLs


# TODO: Verify that none of the SNPs in the earlier sets were dropped for invalid reasons (i.e. bugs
# in the pipeline rather than actual lack of colocalization)
