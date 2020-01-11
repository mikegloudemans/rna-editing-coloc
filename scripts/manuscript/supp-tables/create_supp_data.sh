# Probably should change these output file names to the appropriate supplementary table names

coloc_file="../../../output/colocalization/main_coloc_results/aggregated/aggregated_coloc_results_50snps.txt"
coloc_all_out="../../../output/colocalization/main_coloc_results/aggregated/aggregated_coloc_results_50snps_all.txt"
coloc_edqtl_out="../../../output/colocalization/main_coloc_results/aggregated/aggregated_coloc_results_50snps_edqtl_only.txt"

# edQTL, eQTL, and sQTL colocalization results
#
# May be the same as the original file, but will be sorted
# and will have any edQTLs removed that are not the gene-level
# aggregated ones
cat $coloc_file | grep -v ed_clust | grep -v _edQTLs_txt | sort -k11,11 -k10,10gr > $coloc_all_out

# only edQTL colocalization results
cat $coloc_file | grep -v ed_clust | grep -v _edQTLs_txt | grep Fisher | sort -k11,11 -k10,10gr > $coloc_edqtl_out
