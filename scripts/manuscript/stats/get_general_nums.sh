# How many studies
cat ../../output/aggregated_coloc_results_all.txt | grep Fisher | cut -f11 | sort | uniq | wc -l

# GWAS / locus pairs
cat ../../output/aggregated_coloc_results_all.txt | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f1,11 | lesss | sort | uniq | wc -l

# GWAS / locus / edGene combos
cat ../../output/aggregated_coloc_results_all.txt | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f1,4,11 | lesss | sort | uniq | wc -l

# Number edGenes colocalizing with anything
cat ../../output/aggregated_coloc_results_all.txt | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f4 | lesss | sort | uniq | wc -l

