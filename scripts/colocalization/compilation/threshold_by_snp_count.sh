coloc_out_dir="output/colocalization/rna-editing-revisions/concatenated"

head -n 1 $coloc_out_dir/all_coloc_results.txt > $coloc_out_dir/aggregated_coloc_results_50snps.txt
tail -n +2 $coloc_out_dir/all_coloc_results.txt | sort -k6,6gr | awk '{if ($5 >= 50) print $0}' >> $coloc_out_dir/aggregated_coloc_results_50snps.txt
