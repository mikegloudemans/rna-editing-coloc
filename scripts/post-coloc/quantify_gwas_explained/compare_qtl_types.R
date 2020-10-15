require(gtools)
require(tidyverse)
require(ggplot2)

# duplicate traits that should not be classified as immune nor non-immune because they're
# already present as immune
blacklist = c("Depressive-Symptoms_Okbay_2016_txt_gz")

# Load colocalization matrix
data = read_delim("output/colocalization/rna-editing-revisions/concatenated/aggregated_coloc_results_50snps.txt", delim="\t")
data = data[!(data$base_gwas_file %in% blacklist),]

all_traits_with_coloc_test = unique(data$base_gwas_file)

# Subset to individual QTL types
gene_level_edqtls = data %>% filter(grepl("combined_sorted", eqtl_file))
eqtls = data %>% filter(grepl("eQTLs", eqtl_file))
sqtls = data %>% filter(grepl("sQTLs", eqtl_file))

# Worth noting that there are many more of the other QTL types than there are
# of edQTLs, even combined across all edQTL types
dim(sqtls)
dim(eqtls)
dim(gene_level_edqtls)

# The same colocalization might be tagged multiple times but just with a slightly different ref_snp.
# This could happen if there are multiple subtraits from the same study that all have slightly different lead SNPs, making the
# GWAS-limiting step not so effective. For an admittedly somewhat crude fix to this problem, we bin the coloc tests
# into 1Mb windows by ref_snp, so we'll only have one hit for each GWAS for each general region

gene_level_edqtls$chr = as.numeric(sapply(gene_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
gene_level_edqtls$snp_pos = as.numeric(sapply(gene_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
gene_level_edqtls$locus = floor(gene_level_edqtls$snp_pos / 1000000)

eqtls$chr = as.numeric(sapply(eqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
eqtls$snp_pos = as.numeric(sapply(eqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
eqtls$locus = floor(eqtls$snp_pos / 1000000)

sqtls$chr = as.numeric(sapply(sqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
sqtls$snp_pos = as.numeric(sapply(sqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
sqtls$locus = floor(sqtls$snp_pos / 1000000)

# Load GWAS SNPs so we know the total number of SNPs at each site
gwas = read_delim("output/snp-lists/rna-gwas-overlap_all-gwas_source-pval-min-0_source-pval-max-5e-08_source-window500000_snps-considered.txt", delim="\t")

# Get locus numbers for GWAS
gwas$source_file = sapply(gwas$source_file, function(x) {s=strsplit(x, "/"); return(s[[1]][length(s[[1]])])})
gwas$source_file = gsub("\\.", "_", gwas$source_file)

# We shouldn't be counting the SNPs from a GWAS if its coloc tests all failed to run in practice
gwas = gwas[gwas$source_file %in% all_traits_with_coloc_test,]
gwas = gwas[!(gwas$source_file %in% blacklist),]
gwas$locus = floor(gwas$snp_pos / 1000000)
gwas$full_id = paste(gwas$chr, gwas$locus, sep="_")
colnames(gwas)[colnames(gwas) == "source_file"] = "base_gwas_file"
num_gwas_hits = gwas %>% group_by(base_gwas_file) %>% summarize(num_hits = length(unique(full_id)))


# Get best coloc score for each GWAS SNP
sqtl_coloc = sqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
eqtl_coloc = eqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))
gene_edqtl_coloc = gene_level_edqtls %>% group_by(ref_snp, base_gwas_file, gwas_trait) %>% summarize(best_coloc_h4 = max(clpp_h4))

# Make matrix containing the best COLOC result for each SNP in
# any tissue, any feature of a given QTL type
all_colocs = full_join(sqtl_coloc, eqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_eqtl"))
all_colocs = full_join(all_colocs, gene_edqtl_coloc, by = c("ref_snp", "base_gwas_file", "gwas_trait"), suffix=c("", "_gene_edqtl"))
colnames(all_colocs)[which(colnames(all_colocs) == "best_coloc_h4")] = "best_coloc_h4_sqtl"

all_colocs$chr = as.numeric(sapply(all_colocs$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
all_colocs$snp_pos = as.numeric(sapply(all_colocs$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
all_colocs$locus = floor(all_colocs$snp_pos / 1000000)

# If not tested at all for that type, give it a 0
# to keep things simple.
all_colocs[is.na(all_colocs)] = 0

write.table(all_colocs, file="output/qtl_comparisons/best_colocs_per_locus.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Repeat the above stuff, but allow at most one colocalization per GWAS to avoid
# duplications like the fact that most BMI GWAS have 10-50 traits measured, all stratifications of just BMI

# Get best coloc score for each GWAS SNP
sqtl_coloc_limited = sqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp_sqtl = ref_snp[which.max(clpp_h4)])
eqtl_coloc_limited = eqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp = ref_snp[which.max(clpp_h4)])
gene_edqtl_coloc_limited = gene_level_edqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4) ,best_snp = ref_snp[which.max(clpp_h4)])

all_colocs_limited = full_join(sqtl_coloc_limited, eqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_eqtl"))
all_colocs_limited = full_join(all_colocs_limited, gene_edqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_gene_edqtl"))
colnames(all_colocs_limited)[which(colnames(all_colocs_limited) == "best_coloc_h4")] = "best_coloc_h4_sqtl"

# If not tested at all for that type, give it a 0
# to keep things simple.
all_colocs_limited[,c(4,6,8)][is.na(all_colocs_limited[,c(4,6,8)])] = 0

print(dim(all_colocs))
print(dim(all_colocs_limited))

write.table(all_colocs_limited, file="output/qtl_comparisons/best_colocs_per_locus_limited.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Make the pairs plot
plot_matrix = all_colocs[,4:8]
colnames(plot_matrix) = gsub("best_coloc_h4_", "", colnames(plot_matrix))


pdf("output/plots/qtl_comparisons/qtl_coloc_comparisons_small.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.3)
dev.off()

pdf("output/plots/qtl_comparisons/qtl_coloc_comparisons.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.5)
dev.off()

plot_matrix = apply(plot_matrix, c(1,2), function(x) {max(x, 1e-5)})
plot_matrix = apply(plot_matrix, c(1,2), function(x) {min(x, 1 - 1e-5)})
logit_plot_matrix = logit(plot_matrix)

pdf("output/plots/qtl_comparisons/qtl_coloc_comparisons_logit_small.pdf", width=8, height=8)
	pairs(logit_plot_matrix, pch = 16, cex=0.3)
dev.off()

pdf("output/plots/qtl_comparisons/qtl_coloc_comparisons_logit.pdf", width=8, height=8)
	pairs(logit_plot_matrix, pch = 16, cex=0.5)
dev.off()

logit_plot_matrix = as.data.frame(logit_plot_matrix)
plot_matrix = as.data.frame(plot_matrix)


plot(logit_plot_matrix$eqtl, logit_plot_matrix$gene_edqtl)

# Then make individual plots

# If it works, write up what we did along with a brief analysis
# so we can inspect this further later

edqtl_only = all_colocs_limited %>% filter(best_coloc_h4_sqtl < 0.3, best_coloc_h4_eqtl < 0.3, best_coloc_h4_gene_edqtl > 0.7)

### What fraction of edQTL colocs are eQTL colocs, and what fraction of
### eQTL colocs are edQTLs?

immune_traits = c("Coronary-Artery-Disease_Nelson_2017_txt_gz",
			"Inflammatory-Bowel-Disease_txt_gz",
			"Crohns-Disease_txt_gz",
			"Ulcerative-Colitis_txt_gz", 
			"Lupus_Bentham_2015_txt_gz", 
			"Multiple-Sclerosis_txt_gz",
			"Amyotrophic-Lateral-Sclerosis_Nicolas_2018_txt_gz",
			"Psoriasis_Tsoi_2012_txt_gz",
			"Atopic-Dermatitis_EAGLE_2015_txt_gz",
			"Celiac-Disease_Trynka_2011_txt_gz",
			"Primary-Sclerosing-Cholangitis_Ji_2017_txt_gz",
			"Primary-Biliary-Cirrhosis_Cordell_2015.txt.gz",
			"Type-1-Diabetes-Meta-Analysis_txt_gz")

immune_disease_code = c("CAD", "IBD", "Crohns", "UC", "Lupus", "MS", "ALS", "Psoriasis", "Eczema", "Celiac", "PSC", "PBC", "T1D") 

other_traits = c("High-Density-Lipoprotein-GWAS-and-Metabochip_txt_gz",
			"Low-Density-Lipoprotein-GWAS-and-Metabochip_txt_gz",
			"Total-Cholesterol-GWAS-and-Metabochip_txt_gz",
			"Triglycerides-GWAS-and-Metabochip_txt_gz",
			"Depression-Meta-Analysis_txt_gz",
			"Neuroticism_Luciano_2017_txt_gz",
			"Parkinsons_Pankratz_2012_txt_gz",
			"Type-2-Diabetes-Adjusted-For-BMI-European_txt_gz",
			"Schizophrenia-European_txt_gz",
			"BMI_txt_gz",
			"Height_txt_gz")

other_disease_code = c("HDL", "LDL", "TC", "TG", "Depression", "Neuroticism", "Parkinsons", "T2D", "Schizophrenia", "BMI", "Height")

### Make a stacked barplot ###

coloc_threshold = 0.9

edqtl_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
edqtl_no_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
no_edqtl_no_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)
no_edqtl_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)

table(edqtl_no_eqtl$base_gwas_file)

edqtl_eqtl_percent = dim(edqtl_eqtl)[1]
edqtl_no_eqtl_percent = dim(edqtl_no_eqtl)[1]
no_edqtl_eqtl_percent = dim(no_edqtl_eqtl)[1]
tested_no_edqtl_no_eqtl_percent = dim(no_edqtl_no_eqtl)[1]
untested_no_edqtl_no_eqtl_percent =  sum(num_gwas_hits$num_hits) - sum(edqtl_eqtl_percent, edqtl_no_eqtl_percent,
										    no_edqtl_eqtl_percent, tested_no_edqtl_no_eqtl_percent)
all_counts = c(edqtl_eqtl_percent, edqtl_no_eqtl_percent, no_edqtl_eqtl_percent, tested_no_edqtl_no_eqtl_percent, untested_no_edqtl_no_eqtl_percent)

immune_colocs_limited = all_colocs_limited[all_colocs_limited$base_gwas_file %in% immune_traits,]
immune_gwas_limited = num_gwas_hits[num_gwas_hits$base_gwas_file %in% immune_traits,]

immune_edqtl_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
immune_edqtl_no_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
immune_no_edqtl_no_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)
immune_no_edqtl_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)

dim(immune_edqtl_eqtl)[1]
dim(immune_edqtl_no_eqtl)[1]
dim(immune_no_edqtl_no_eqtl)[1]
dim(immune_no_edqtl_eqtl)[1]

table(immune_edqtl_no_eqtl$base_gwas_file)

immune_edqtl_eqtl_percent = dim(immune_edqtl_eqtl)[1]
immune_edqtl_no_eqtl_percent = dim(immune_edqtl_no_eqtl)[1]
immune_no_edqtl_eqtl_percent = dim(immune_no_edqtl_eqtl)[1]
immune_tested_no_edqtl_no_eqtl_percent = dim(immune_no_edqtl_no_eqtl)[1]
immune_untested_no_edqtl_no_eqtl_percent =  sum(immune_gwas_limited$num_hits) - sum(immune_edqtl_eqtl_percent, immune_edqtl_no_eqtl_percent,
										    immune_no_edqtl_eqtl_percent, immune_tested_no_edqtl_no_eqtl_percent)
immune_counts = c(immune_edqtl_eqtl_percent, immune_edqtl_no_eqtl_percent, immune_no_edqtl_eqtl_percent, immune_tested_no_edqtl_no_eqtl_percent, immune_untested_no_edqtl_no_eqtl_percent)

disease_specific_colocs = list()
for (disease in immune_traits)
{
	disease_limited = all_colocs_limited[all_colocs_limited$base_gwas_file == disease,]
	disease_gwas_limited = num_gwas_hits[num_gwas_hits$base_gwas_file == disease,]

	disease_edqtl_eqtl = disease_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
	disease_edqtl_no_eqtl = disease_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl > coloc_threshold)
	disease_no_edqtl_no_eqtl = disease_limited %>% filter(best_coloc_h4_eqtl <= coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)
	disease_no_edqtl_eqtl = disease_limited %>% filter(best_coloc_h4_eqtl > coloc_threshold, best_coloc_h4_gene_edqtl <= coloc_threshold)

	disease_specific_colocs[[disease]] = list()
	disease_specific_colocs[[disease]][["edqtl_eqtl_percent"]] = dim(disease_edqtl_eqtl)[1]
	disease_specific_colocs[[disease]][["edqtl_no_eqtl_percent"]] = dim(disease_edqtl_no_eqtl)[1]
	disease_specific_colocs[[disease]][["no_edqtl_eqtl_percent"]] = dim(disease_no_edqtl_eqtl)[1]
	disease_specific_colocs[[disease]][["tested_no_edqtl_no_eqtl_percent"]] = dim(disease_no_edqtl_no_eqtl)[1]
	disease_specific_colocs[[disease]][["untested_no_edqtl_no_eqtl_percent"]] = sum(disease_gwas_limited$num_hits) - sum(disease_specific_colocs[[disease]][["edqtl_eqtl_percent"]],
										   disease_specific_colocs[[disease]][["edqtl_no_eqtl_percent"]],
										   disease_specific_colocs[[disease]][["no_edqtl_eqtl_percent"]],
										   disease_specific_colocs[[disease]][["tested_no_edqtl_no_eqtl_percent"]])

	disease_specific_colocs[[disease]][["counts"]] = c(disease_specific_colocs[[disease]][["edqtl_eqtl_percent"]],
							   disease_specific_colocs[[disease]][["edqtl_no_eqtl_percent"]],
							   disease_specific_colocs[[disease]][["no_edqtl_eqtl_percent"]],
							   disease_specific_colocs[[disease]][["tested_no_edqtl_no_eqtl_percent"]],
							   disease_specific_colocs[[disease]][["untested_no_edqtl_no_eqtl_percent"]])
}


subset = rep(c("all", "autoimmune", immune_disease_code), each=5)
subset = factor(subset, levels=c("all", "autoimmune", immune_disease_code))

condition = rep(c("edQTL and eQTL coloc", "edQTL coloc only", "eQTL coloc only", "QTL, no coloc", "no QTL"), 2+length(immune_traits))
condition = factor(condition, levels=rev(c("eQTL coloc only", "edQTL and eQTL coloc", "edQTL coloc only", "QTL, no coloc", "no QTL")))
counts = c(all_counts, immune_counts, as.vector(sapply(immune_traits, function(x)
					     {
						     disease_specific_colocs[[x]][["counts"]]
					     })))
stack_data = data.frame(subset, condition, counts)

g = ggplot(stack_data, aes(fill=condition, y=counts, x=subset)) +
	geom_bar(position="fill", stat="identity") + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1))# +
	#geom_text(aes(label=round(counts,2)), vjust=1.6, color="white", size=3.5)
g

###############################################
# A different approach 
# Forget about eQTLs, just focus on edQTLs
###############################################

edqtl_coloc = gene_edqtl_coloc_limited %>% filter(best_coloc_h4 > coloc_threshold)
no_edqtl_coloc = gene_edqtl_coloc_limited %>% filter(best_coloc_h4 <= coloc_threshold)

edqtl_count = dim(edqtl_coloc)[1]
tested_no_edqtl_count = dim(no_edqtl_coloc)[1]
untested_no_edqtl_count =  sum(num_gwas_hits$num_hits) - edqtl_count - tested_no_edqtl_count
all_count = c(edqtl_count, tested_no_edqtl_count, untested_no_edqtl_count)
all_percents = all_count / sum(all_count)

nonimmune_edqtl_colocs_limited = gene_edqtl_coloc_limited[!(gene_edqtl_coloc_limited$base_gwas_file %in% immune_traits),]
nonimmune_gwas_limited = num_gwas_hits[!(num_gwas_hits$base_gwas_file %in% immune_traits),]

nonimmune_edqtl = nonimmune_edqtl_colocs_limited %>% filter(best_coloc_h4 > coloc_threshold)
nonimmune_no_edqtl = nonimmune_edqtl_colocs_limited %>% filter(best_coloc_h4 <= coloc_threshold)

nonimmune_edqtl_count = dim(nonimmune_edqtl)[1]
nonimmune_tested_no_edqtl_count = dim(nonimmune_no_edqtl)[1]
nonimmune_untested_no_edqtl_count =  sum(nonimmune_gwas_limited$num_hits) - nonimmune_edqtl_count - nonimmune_tested_no_edqtl_count

nonimmune_count = c(nonimmune_edqtl_count, nonimmune_tested_no_edqtl_count, nonimmune_untested_no_edqtl_count)
nonimmune_percents = nonimmune_count / sum(nonimmune_count)

immune_edqtl_colocs_limited = gene_edqtl_coloc_limited[gene_edqtl_coloc_limited$base_gwas_file %in% immune_traits,]

immune_edqtl = immune_edqtl_colocs_limited %>% filter(best_coloc_h4 > coloc_threshold)
immune_no_edqtl = immune_edqtl_colocs_limited %>% filter(best_coloc_h4 <= coloc_threshold)

immune_edqtl_count = dim(immune_edqtl)[1]
immune_tested_no_edqtl_count = dim(immune_no_edqtl)[1]
immune_untested_no_edqtl_count =  sum(immune_gwas_limited$num_hits) - immune_edqtl_count - immune_tested_no_edqtl_count

immune_count = c(immune_edqtl_count, immune_tested_no_edqtl_count, immune_untested_no_edqtl_count)
immune_percents = immune_count / sum(immune_count)

disease_specific_edqtl_colocs = list()
for (disease in c(immune_traits, other_traits))
{
	disease_limited = gene_edqtl_coloc_limited[gene_edqtl_coloc_limited$base_gwas_file == disease,]
	disease_gwas_limited = num_gwas_hits[num_gwas_hits$base_gwas_file == disease,]

	disease_edqtl = disease_limited %>% filter(best_coloc_h4 > coloc_threshold)
	disease_no_edqtl = disease_limited %>% filter(best_coloc_h4 <= coloc_threshold)

	disease_specific_edqtl_colocs[[disease]] = list()
	disease_specific_edqtl_colocs[[disease]][["edqtl_count"]] = dim(disease_edqtl)[1]
	disease_specific_edqtl_colocs[[disease]][["tested_no_edqtl_count"]] = dim(disease_no_edqtl)[1]
	disease_specific_edqtl_colocs[[disease]][["untested_no_edqtl_count"]] = sum(disease_gwas_limited$num_hits) - disease_specific_edqtl_colocs[[disease]][["edqtl_count"]] - disease_specific_edqtl_colocs[[disease]][["tested_no_edqtl_count"]]

	disease_specific_edqtl_colocs[[disease]][["count"]] = c(disease_specific_edqtl_colocs[[disease]][["edqtl_count"]],
							   disease_specific_edqtl_colocs[[disease]][["tested_no_edqtl_count"]],
							   disease_specific_edqtl_colocs[[disease]][["untested_no_edqtl_count"]])
	disease_specific_edqtl_colocs[[disease]][["percents"]] = disease_specific_edqtl_colocs[[disease]][["count"]] / sum(disease_specific_edqtl_colocs[[disease]][["count"]])
}

subset = rep(c("non-immune", other_disease_code, "immune-related", immune_disease_code), each=3)
subset = factor(subset, levels=c("non-immune", other_disease_code, "immune-related", immune_disease_code))

condition = rep(c("edQTL coloc", "edQTL, no coloc", "no edQTL"), 2+length(immune_traits) + length(other_traits))
condition = factor(condition, levels=rev(c("edQTL coloc", "edQTL, no coloc", "no edQTL")))
percent_hits_explained = c(nonimmune_percents,
	   as.vector(sapply(other_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["percents"]]
					     })),
	   immune_percents, 
	   as.vector(sapply(immune_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["percents"]]
					     })))
stack_data = data.frame(subset, condition, percent_hits_explained)

g = ggplot(stack_data, aes(fill=condition, y=percent_hits_explained, x=subset)) +
	#geom_bar(position="fill", stat="identity") + 
	geom_bar(stat="identity") + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1))# +
	#geom_text(aes(label=round(counts,2)), vjust=1.6, color="white", size=3.5)
g

##################################################
# Replot it without the non-GWAS hits
##################################################

immune_order = rev(order(sapply(immune_traits, function(x) {
				     sum(disease_specific_edqtl_colocs[[x]][["percents"]][c(1,2)])
			     	})))
other_order = rev(order(sapply(other_traits, function(x) {
				     sum(disease_specific_edqtl_colocs[[x]][["percents"]][c(1,2)])
			     	})))


disease = rep(c("non-immune", other_disease_code, "immune-related", immune_disease_code), each=2)
disease = factor(disease, levels=c("non-immune", other_disease_code[other_order], "immune-related", immune_disease_code[immune_order]))

condition = rep(c("edQTL, coloc", "edQTL, no coloc"), 2+length(immune_traits) + length(other_traits))
condition = factor(condition, levels=rev(c("edQTL, coloc", "edQTL, no coloc")))
percent_hits_explained = c(nonimmune_percents[c(1,2)],
	   as.vector(sapply(other_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["percents"]][c(1,2)]
					     })),
	   immune_percents[c(1,2)], 
	   as.vector(sapply(immune_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["percents"]][c(1,2)]
					     })))
total_percent_hits_explained = percent_hits_explained
total_percent_hits_explained[(1:(length(total_percent_hits_explained)/2))*2] = total_percent_hits_explained[(1:(length(total_percent_hits_explained)/2))*2] +
								total_percent_hits_explained[(1:(length(total_percent_hits_explained)/2))*2-1]
total_percent_hits_explained[(1:(length(total_percent_hits_explained)/2))*2-1] = 0

count_hits_explained = c(nonimmune_count[c(1,2)],
	   as.vector(sapply(other_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["count"]][c(1,2)]
					     })),
	   immune_count[c(1,2)], 
	   as.vector(sapply(immune_traits, function(x)
					     {
						     disease_specific_edqtl_colocs[[x]][["count"]][c(1,2)]
					     })))
count_hits_explained[(1:(length(count_hits_explained)/2))*2] = count_hits_explained[(1:(length(count_hits_explained)/2))*2] +
								count_hits_explained[(1:(length(count_hits_explained)/2))*2-1]
count_hits_explained[(1:(length(count_hits_explained)/2))*2-1] = 0
count_hits_explained = as.character(count_hits_explained)
count_hits_explained[count_hits_explained == "0"] = ""

stack_data = data.frame(disease, condition, percent_hits_explained, count_hits_explained)

g = ggplot(stack_data, aes(fill=condition, y=percent_hits_explained, x=disease), colour="black") +
	theme_classic() +
	scale_fill_manual(values=c("gray80", "gray60")) +
	#geom_bar(position="fill", stat="identity") + 
	geom_bar(stat="identity") + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1, face=ifelse(levels(disease)=="non-immune" | levels(disease)=="immune-related","bold","plain"))) +
	theme(legend.title=element_blank()) + 
	xlab("GWAS trait") +
	ylab("% GWAS hits explained") +
	geom_text(aes(label=count_hits_explained, y=total_percent_hits_explained+0.02), vjust=1.6, color="black", size=3.5)
g

ggsave("output/plots/quantify_gwas_explained/panel3a.pdf", height=6, width=6)

stopifnot(FALSE)

####
all_disease_edqtl_colocs = list()
for (disease in unique(gene_edqtl_coloc_limited$base_gwas_file))
{
	disease_limited = gene_edqtl_coloc_limited[gene_edqtl_coloc_limited$base_gwas_file == disease,]
	disease_gwas_limited = num_gwas_hits[num_gwas_hits$base_gwas_file == disease,]

	disease_edqtl = disease_limited %>% filter(best_coloc_h4 > coloc_threshold)
	disease_no_edqtl = disease_limited %>% filter(best_coloc_h4 <= coloc_threshold)

	all_disease_edqtl_colocs[[disease]] = list()
	all_disease_edqtl_colocs[[disease]][["edqtl_count"]] = dim(disease_edqtl)[1]
	all_disease_edqtl_colocs[[disease]][["tested_no_edqtl_count"]] = dim(disease_no_edqtl)[1]
	all_disease_edqtl_colocs[[disease]][["untested_no_edqtl_count"]] = sum(disease_gwas_limited$num_hits) - all_disease_edqtl_colocs[[disease]][["edqtl_count"]] - all_disease_edqtl_colocs[[disease]][["tested_no_edqtl_count"]]

	all_disease_edqtl_colocs[[disease]][["count"]] = c(all_disease_edqtl_colocs[[disease]][["edqtl_count"]],
							   all_disease_edqtl_colocs[[disease]][["tested_no_edqtl_count"]],
							   all_disease_edqtl_colocs[[disease]][["untested_no_edqtl_count"]])
	all_disease_edqtl_colocs[[disease]][["percents"]] = all_disease_edqtl_colocs[[disease]][["count"]] / sum(all_disease_edqtl_colocs[[disease]][["count"]])
}

subset = rep(unique(gene_edqtl_coloc_limited$base_gwas_file), each=3)
subset = factor(subset, levels= unique(gene_edqtl_coloc_limited$base_gwas_file))

condition = rep(c("edQTL, coloc", "edQTL, no coloc", "no edQTL"), length(unique(gene_edqtl_coloc_limited$base_gwas_file)))
condition = factor(condition, levels=rev(c("edQTL, coloc", "edQTL, no coloc", "no edQTL")))

percent_hits_explained = as.vector(sapply(unique(gene_edqtl_coloc_limited$base_gwas_file), function(x)
					     {
						     all_disease_edqtl_colocs[[x]][["percents"]]
					     }))
count_hits_explained = as.vector(sapply(unique(gene_edqtl_coloc_limited$base_gwas_file), function(x)
					     {
						     all_disease_edqtl_colocs[[x]][["count"]]
					     }))



stack_data = data.frame(subset, condition, percent_hits_explained, count_hits_explained)

write.table(stack_data, file="output/plot-data/edqtl_coloc_info.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

