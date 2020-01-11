# Get a few general counts of features, for inclusion in the main text

coloc_file="../../../output/colocalization/main_coloc_results/aggregated/aggregated_coloc_results_50snps.txt"


# How many studies
cat $coloc_file | grep Fisher | cut -f11 | sort | uniq | wc -l

# GWAS / locus pairs
cat $coloc_file | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f1,11 | lesss | sort | uniq | wc -l

# GWAS / locus / edGene combos
cat $coloc_file | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f1,4,11 | lesss | sort | uniq | wc -l

# Number edGenes colocalizing with anything
cat $coloc_file | grep Fisher | awk '{if ($10 > 0.5) print $0}' |  cut -f4 | lesss | sort | uniq | wc -l

