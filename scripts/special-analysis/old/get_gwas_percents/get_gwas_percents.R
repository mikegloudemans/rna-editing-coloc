# Author: Mike Gloudemans

require(readr)
require(dplyr)
require(tidyr)
require(reshape2)
require(ggplot2)

# Load files

#
# Get number of GWAS SNPS
#

# Load file showing all significant GWAS SNPs for a locus
all_gwas_snps = read_delim("/users/mgloud/projects/rna_editing/output/test-snps/rna-editing_all-gwas_gwas-pval1e-06_gwas-window500000_snps-considered.txt", delim = "\t", skip=1, col_names=c("chr", "snp_pos", "pvalue", "trait", "gwas_file"))
all_gwas_snps$id = paste(all_gwas_snps$chr, all_gwas_snps$snp_pos, sep="_")
all_gwas_snps$full_trait = paste(all_gwas_snps$gwas_file, all_gwas_snps$trait, sep="_")
all_gwas_snps$full_trait = gsub("/users/mgloud/projects/rna_editing/data/gwas/", "", all_gwas_snps$full_trait)
all_gwas_snps$full_trait = gsub("\\.", "_", all_gwas_snps$full_trait)

gwas_hit_counts = all_gwas_snps %>% group_by(full_trait) %>% summarize(total = length(unique(id)))
dim(gwas_hit_counts)

#
# Get number of colocalization tests in various QTL types...
#

coloc_results = read_delim("/users/mgloud/projects/rna_editing/output/aggregated_clpp_results2.txt", delim="\t")
coloc_results$full_trait = paste(coloc_results$base_gwas_file, coloc_results$gwas_trait, sep="_")
coloc_results$full_trait = gsub("\\.", "_", coloc_results$full_trait)
dim(coloc_results)
coloc_results = coloc_results[coloc_results$clpp_h4 > 0.5,]
dim(coloc_results)

edqtl_results = coloc_results[grepl("edQTLs",coloc_results$eqtl_file),]
edqtl_aggro_results = coloc_results[grepl("Fisher",coloc_results$eqtl_file),]
eqtl_results = coloc_results[grepl("eQTLs",coloc_results$eqtl_file),]
sqtl_results = coloc_results[grepl("sQTLs",coloc_results$eqtl_file),]

edqtl_coloc_counts = edqtl_results %>% group_by(full_trait) %>% summarize(total = length(unique(ref_snp)))
edqtl_aggro_coloc_counts = edqtl_aggro_results %>% group_by(full_trait) %>% summarize(total = length(unique(ref_snp)))

# eQTLs
eqtl_coloc_counts = eqtl_results %>% group_by(full_trait) %>% summarize(total = length(unique(ref_snp)))

# sQTLs
sqtl_coloc_counts = sqtl_results %>% group_by(full_trait) %>% summarize(total = length(unique(ref_snp)))

# Merge all into one table...fill the ones not appearing with 0's
all_counts = full_join(gwas_hit_counts, edqtl_coloc_counts, by="full_trait", suffix = c("", "_edQTL"))
all_counts = full_join(all_counts, edqtl_aggro_coloc_counts, by="full_trait", suffix = c("", "_edQTL_aggro"))
all_counts = full_join(all_counts, eqtl_coloc_counts, by="full_trait", suffix = c("", "_eQTL"))
all_counts = full_join(all_counts, sqtl_coloc_counts, by="full_trait", suffix = c("", "_sQTL"))

#
# Compute percentages explained by each of these different QTL types
#
all_counts = all_counts %>% replace(is.na(.), 0)

all_counts$percent_edQTL = all_counts$total_edQTL / all_counts$total
all_counts$percent_edQTL_aggro = all_counts$total_edQTL_aggro / all_counts$total
all_counts$percent_eQTL = all_counts$total_eQTL / all_counts$total
all_counts$percent_sQTL = all_counts$total_sQTL / all_counts$total

# Let's say there must be at least 10 GWAS hits for us to have
# any confidence in the results
enough_counts = all_counts[all_counts$total >= 10,]

# For quick viewing
rbind(enough_counts$full_trait, all_counts$percent_edQTL_aggro / all_counts$percent_eQTL, all_counts$total)

#
# Plot in a nice bar plot, like the one in GTEx main paper for v6p
#

plot_data = enough_counts[c("full_trait", "percent_edQTL", "percent_edQTL_aggro", "percent_eQTL", "percent_sQTL")]
order_by_max = plot_data[order(-apply(plot_data[,2:5], 1, max)),]
order_by_ratio = plot_data[order(-(plot_data$percent_edQTL + plot_data$percent_edQTL_aggro) / plot_data$percent_eQTL),]
order_by_editing = plot_data[order(-(plot_data$percent_edQTL + plot_data$percent_edQTL_aggro)),]

melted = melt(order_by_max)
melted$full_trait = factor(melted$full_trait, levels = order_by_max$full_trait)

colnames(melted) = c("trait", "qtl_type", "percent_colocalized")	
g = ggplot(data=melted, aes(x=trait, y=percent_colocalized, fill=qtl_type)) +
	geom_bar(stat="identity", position=position_dodge()) +
	theme_bw() +
	theme(plot.margin = unit(c(5,5,5,5), "cm")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("GWAS trait") +
	ylab("Percent colocalized") +
	scale_fill_brewer(palette="Set1", name="QTL type")
g
ggsave(paste0("/users/mgloud/projects/rna_editing/output/plots/get_gwas_percents/sorted_by_max.png"), width=24, height=12)

melted = melt(order_by_ratio)
melted$full_trait = factor(melted$full_trait, levels = order_by_ratio$full_trait)

colnames(melted) = c("trait", "qtl_type", "percent_colocalized")	
g = ggplot(data=melted, aes(x=trait, y=percent_colocalized, fill=qtl_type)) +
	geom_bar(stat="identity", position=position_dodge()) +
	theme_bw() +
	theme(plot.margin = unit(c(5,5,5,5), "cm")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("GWAS trait") +
	ylab("Percent colocalized") +
	scale_fill_brewer(palette="Set1", name="QTL type")
g
ggsave(paste0("/users/mgloud/projects/rna_editing/output/plots/get_gwas_percents/sorted_by_ratio.png"), width=24, height=12)

melted = melt(order_by_editing)
melted$full_trait = factor(melted$full_trait, levels = order_by_editing$full_trait)

colnames(melted) = c("trait", "qtl_type", "percent_colocalized")	
g = ggplot(data=melted, aes(x=trait, y=percent_colocalized, fill=qtl_type)) +
	geom_bar(stat="identity", position=position_dodge()) +
	theme_bw() +
	theme(plot.margin = unit(c(5,5,5,5), "cm")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("GWAS trait") +
	ylab("Percent colocalized") +
	scale_fill_brewer(palette="Set1", name="QTL type")
g
ggsave(paste0("/users/mgloud/projects/rna_editing/output/plots/get_gwas_percents/sorted_by_editing.png"), width=24, height=12)




# (this will be easy to do using melt and ggplot2)

# TODO: Figure out why sQTLs are missing AND figure out why a lot of GWAS aren't actually present here...
# Then do analyses for the other ones too...
