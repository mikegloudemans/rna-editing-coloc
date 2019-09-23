require(gtools)
require(tidyverse)

# Load colocalization matrix
data = read_delim("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", delim="\t")

# Subset to individual QTL types
gene_level_edqtls = data %>% filter(grepl("combined_sorted", eqtl_file))
cluster_level_edqtls = data %>% filter(grepl("clust", eqtl_file))
site_level_edqtls = data %>% filter(grepl("edQTLs", eqtl_file))
eqtls = data %>% filter(grepl("eQTLs", eqtl_file))
sqtls = data %>% filter(grepl("sQTLs", eqtl_file))

# Worth noting that there are many more of the other QTL types than there are
# of edQTLs, even combined across all edQTL types
dim(sqtls)
dim(eqtls)
dim(gene_level_edqtls)
dim(cluster_level_edqtls)
dim(site_level_edqtls)

# Get best coloc score for each GWAS SNP
sqtl_coloc = sqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
eqtl_coloc = eqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
gene_edqtl_coloc = gene_level_edqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
cluster_edqtl_coloc = cluster_level_edqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
site_edqtl_coloc = site_level_edqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))

# Make matrix containing the best COLOC result for each SNP in
# any tissue, any feature of a given QTL type
all_colocs = full_join(sqtl_coloc, eqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_eqtl"))
all_colocs = full_join(all_colocs, gene_edqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_gene_edqtl"))
all_colocs = full_join(all_colocs, cluster_edqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_cluster_edqtl"))
all_colocs = full_join(all_colocs, site_edqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_site_edqtl"))
colnames(all_colocs)[which(colnames(all_colocs) == "best_coloc_h4")] = "best_coloc_h4_sqtl"

# If not tested at all for that type, give it NA, but we'll treat
# that as basically a 0.
all_colocs[is.na(all_colocs)] = 0

# Most values are either very close to 1 or very close to 0...
hist(all_colocs$best_coloc_h4_sqtl)

# Make the pairs plot
plot_matrix = all_colocs[,4:8]
colnames(plot_matrix) = gsub("best_coloc_h4_", "", colnames(plot_matrix))


pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_small.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.2)
dev.off()

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.5)
dev.off()

plot_matrix = apply(plot_matrix, c(1,2), function(x) {max(x, 1e-5)})
plot_matrix = apply(plot_matrix, c(1,2), function(x) {min(x, 1 - 1e-5)})
logit_plot_matrix = logit(plot_matrix)

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_logit_small.pdf", width=8, height=8)
	pairs(logit_plot_matrix, pch = 16, cex=0.2)
dev.off()

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_logit.pdf", width=8, height=8)
	pairs(logit_plot_matrix, pch = 16, cex=0.5)
dev.off()



# Then make individual plots

# If it works, write up what we did along with a brief analysis
# so we can inspect this further later
