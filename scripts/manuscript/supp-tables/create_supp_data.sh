cat  ../output/aggregated_coloc_results_all.txt | grep -v ed_clust | grep -v _edQTLs_txt | sort -k11,11 -k10,10gr > ../output/aggregated_coloc_results_all_filtered.txt
 cat  ../output/aggregated_coloc_results_all.txt | grep -v ed_clust | grep -v _edQTLs_txt | grep Fisher | sort -k11,11 -k10,10gr > ../output/aggregated_coloc_results_all_filtered_edqtl.txt
