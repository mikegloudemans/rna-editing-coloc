require(readr)
require(dplyr)

# Compare the colocalizations found using aggregated (by gene) edQTLs
# and single-position edQTLs

# TODO: Right now it's hard to tell whether differences are due to
# difference in method used or difference in actual power/signal.
# Try running the single-SNP analysis with COLOC so that we can
# get a nicer comparison for this part.

# Load in both coloc output files
single_site_coloc = read_delim("/users/mgloud/projects/rna_editing/output/single_site_clpp_results.tsv", delim="\t")
aggro_coloc = read_delim("/users/mgloud/projects/rna_editing/output/aggregated_clpp_results.txt", delim="\t")
single_site_coloc$chr = sapply(single_site_coloc$feature, function(x) {strsplit(x, "_")[[1]][1]})
single_site_coloc$snp_pos = as.numeric(sapply(single_site_coloc$feature, function(x) {strsplit(x, "_")[[1]][2]}))

# Match the single-position edQTLs to the relevant genes
# (just guess for now, and we'll get the true mapping from
# Qin later if needed)
genes = read_delim("/users/mgloud/projects/rna_editing/scripts/compare_aggro_with_single/gencode.v30.annotation.gtf.gz", delim="\t", skip=5, col_names=FALSE)
genes = genes[genes$X3 == "gene",]
genes$gene = sapply(genes$X9, function(x) {strsplit(strsplit(x, "gene_name \"")[[1]][[2]], "\"")[[1]][1]})
single_site_coloc$gene = sapply(1:dim(single_site_coloc)[1], function(i)
{
	gene_chr = genes[genes$X1 == single_site_coloc$chr[i],]
	gene_sub = gene_chr[gene_chr$X4 <= single_site_coloc$snp_pos[i] & gene_chr$X5 >= single_site_coloc$snp_pos[i],]
	return(paste(gene_sub$gene, collapse=";"))
})

single_site_coloc = single_site_coloc[single_site_coloc$n_snps >= 50,]
aggro_coloc = aggro_coloc[aggro_coloc$n_snps >= 50,]

aggro_coloc$coloc_rank = rank(-aggro_coloc$clpp_h4) 
single_site_coloc$coloc_rank = rank(-single_site_coloc$clpp) 

# TODO:
# For both studies...
# Threshold at some point on what colocalizations are considered
# significant (e.g. H4 > CLPP)
dim(single_site_coloc)
single_site_passing = single_site_coloc[!is.na(single_site_coloc$clpp),]
single_site_passing = single_site_passing[single_site_passing$clpp > 0.01,]
dim(single_site_passing)

dim(aggro_coloc)
aggro_passing = aggro_coloc[!is.na(aggro_coloc$clpp_h4),]
aggro_passing = aggro_passing[aggro_passing$clpp_h4 > 0.5,]
dim(aggro_passing)

# If we consider the best colocalizing test for each
# gene from each of the two approaches,
# how do their overall ranks differ?
#
# Only count them if tested in both methods and
# if single site is uniquely assigned to one gene
aggro_sum = aggro_coloc %>% group_by(feature) %>% summarize(coloc_status = min(coloc_rank))
single_site_sum = single_site_coloc %>% group_by(gene) %>% summarize(coloc_status = min(coloc_rank))
tested_in_both = merge(single_site_sum, aggro_sum, by.x="gene", by.y="feature")

dim(aggro_sum)
dim(single_site_sum)
dim(tested_in_both)

# Plot comparison between the two studies and identify drastically different 
# SNPs (i.e. highly significant in one method but not in the other)
plot(tested_in_both$coloc_status.x, tested_in_both$coloc_status.y, xlab="Gene rank in single site analysis", ylab="Gene rank in aggregated analysis", pch=16)
cor.test(tested_in_both$coloc_status.x, tested_in_both$coloc_status.y, method="spearman")

# Draw threshold lines on this plot
abline(a=rank(c(-0.5, -aggro_coloc$clpp_h4))[1], b=0, lty=3, col="red")
abline(v=rank(c(-0.01, -single_site_coloc$clpp))[1], b=0, lty=3, col="red")

# Highlight genes uniquely found in one method
single_site_only = tested_in_both[(tested_in_both$coloc_status.x <= rank(c(-0.01, -single_site_coloc$clpp))[1]) & (tested_in_both$coloc_status.y > rank(c(-0.5, -aggro_coloc$clpp_h4))[1]),]
points(single_site_only$coloc_status.x, single_site_only$coloc_status.y, pch=16, col="red")
aggro_only = tested_in_both[(tested_in_both$coloc_status.x > rank(c(-0.01, -single_site_coloc$clpp))[1]) & (tested_in_both$coloc_status.y <= rank(c(-0.5, -aggro_coloc$clpp_h4))[1]),]
points(aggro_only$coloc_status.x, aggro_only$coloc_status.y, pch=16, col="blue")
both_detected = tested_in_both[(tested_in_both$coloc_status.x <= rank(c(-0.01, -single_site_coloc$clpp))[1]) & (tested_in_both$coloc_status.y <= rank(c(-0.5, -aggro_coloc$clpp_h4))[1]),]
neither_detected = tested_in_both[(tested_in_both$coloc_status.x > rank(c(-0.01, -single_site_coloc$clpp))[1]) & (tested_in_both$coloc_status.y > rank(c(-0.5, -aggro_coloc$clpp_h4))[1]),]

# NOTE: Right now some of the differences may have more to do with the method used than the actual
# underlying biology! We'll run it again using COLOC for both methods, to ensure comparability.

# How many unique to each method? How many agree?
dim(both_detected)
dim(neither_detected)
dim(single_site_only)
dim(aggro_only)

# Examine a few of these sites up close...show the best colocalization
# from either of the two methods (in PPT), one should be much better

# How many genes were tested in only ONE method?
# Of those, how many were significant?
single_test_only = single_site_sum[(!grepl(";", single_site_sum$gene)) & (!grepl(",", single_site_sum$gene)) & (!grepl("Sep", single_site_sum$gene)) & (!grepl("MARCH", single_site_sum$gene)) & (nchar(single_site_sum$gene) > 0) & (sapply(single_site_sum$gene, function(x) { sum(grepl(x, aggro_sum$feature)) == 0 })),]
single_test_only_coloc = single_test_only[single_test_only$coloc_status <= rank(c(-0.01, -single_site_coloc$clpp))[1],]
dim(single_test_only_coloc)

aggro_test_only = aggro_sum[(!grepl(";", aggro_sum$feature)) & (!grepl(",", aggro_sum$feature)) & (!grepl("Sep", aggro_sum$feature)) & (!grepl("Mar", aggro_sum$feature)) & (nchar(aggro_sum$feature) > 0) & (sapply(aggro_sum$feature, function(x) { sum(grepl(x, single_site_sum$gene)) == 0 })),]
aggro_test_only_coloc = aggro_test_only[aggro_test_only$coloc_status <= rank(c(-0.5, -aggro_coloc$clpp_h4))[1],]
dim(aggro_test_only_coloc)

# Takeaway seems to be that while methods agree in many ways, there's
# certainly some value to including both.

