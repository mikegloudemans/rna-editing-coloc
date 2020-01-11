coloc_out_dir="../../../output/colocalization/main_coloc_results/aggregated/"

head -n 1 $coloc_out_dir/aggregated_coloc_results.txt > $coloc_out_dir/aggregated_coloc_results_50snps.txt
tail -n +2 $coloc_out_dir/aggregated_coloc_results.txt | sort -k6,6gr | awk '{if ($5 > 50) print $0}' >> $coloc_out_dir/aggregated_coloc_results_50snps.txt
