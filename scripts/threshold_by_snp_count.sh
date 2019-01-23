head -n 1 ../output/single_site_clpp_results.tsv > ../output/single_site_clpp_results_snp_count_thresholded.tsv
tail -n +2 ../output/single_site_clpp_results.tsv | sort -k6,6gr | awk '{if ($5 > 50) print $0}' >> ../output/single_site_clpp_results_snp_count_thresholded.tsv
