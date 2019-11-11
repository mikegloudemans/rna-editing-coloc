require(gtools)
require(tidyverse)
require(ggplot2)

# TODO: Figure out a way to deal with the "multiple traits in a GWAS" issue
# Maybe just take the best trait from each base file?

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

# The same colocalization might be tagged multiple times but just with a slightly different ref_snp.
# This could happen if there are multiple subtraits from the same study that all have slightly different lead SNPs, making the
# GWAS-limiting step not so effective. For an admittedly somewhat crude fix to this problem, we bin the coloc tests
# into 1Mb windows by ref_snp, so we'll only have one hit for each GWAS for each general region

gene_level_edqtls$chr = as.numeric(sapply(gene_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
gene_level_edqtls$snp_pos = as.numeric(sapply(gene_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
gene_level_edqtls$locus = floor(gene_level_edqtls$snp_pos / 1000000)

cluster_level_edqtls$chr = as.numeric(sapply(cluster_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
cluster_level_edqtls$snp_pos = as.numeric(sapply(cluster_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
cluster_level_edqtls$locus = floor(cluster_level_edqtls$snp_pos / 1000000)

site_level_edqtls$chr = as.numeric(sapply(site_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
site_level_edqtls$snp_pos = as.numeric(sapply(site_level_edqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
site_level_edqtls$locus = floor(site_level_edqtls$snp_pos / 1000000)

eqtls$chr = as.numeric(sapply(eqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
eqtls$snp_pos = as.numeric(sapply(eqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
eqtls$locus = floor(eqtls$snp_pos / 1000000)

sqtls$chr = as.numeric(sapply(sqtls$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
sqtls$snp_pos = as.numeric(sapply(sqtls$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
sqtls$locus = floor(sqtls$snp_pos / 1000000)

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

all_colocs$chr = as.numeric(sapply(all_colocs$ref_snp, function(x){strsplit(x, "_")[[1]][1]}))
all_colocs$snp_pos = as.numeric(sapply(all_colocs$ref_snp, function(x){strsplit(x, "_")[[1]][2]}))
all_colocs$locus = floor(all_colocs$snp_pos / 1000000)

# If not tested at all for that type, give it a 0
# to keep things simple.
all_colocs[is.na(all_colocs)] = 0

# Most values are either very close to 1 or very close to 0...
hist(all_colocs$best_coloc_h4_sqtl)

write.table(all_colocs, file="/users/mgloud/projects/rna_editing/output/qtl_comparisons/best_colocs_per_locus.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Repeat the above stuff, but allow at most one colocalization per GWAS to avoid
# duplications like the fact that most BMI GWAS have 10-50 traits measured, all stratifications of just BMI

# Get best coloc score for each GWAS SNP
sqtl_coloc_limited = sqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp_sqtl = ref_snp[which.max(clpp_h4)])
eqtl_coloc_limited = eqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp = ref_snp[which.max(clpp_h4)])
gene_edqtl_coloc_limited = gene_level_edqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4) ,best_snp = ref_snp[which.max(clpp_h4)])
cluster_edqtl_coloc_limited = cluster_level_edqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp = ref_snp[which.max(clpp_h4)])
site_edqtl_coloc_limited = site_level_edqtls %>% group_by(chr, locus, base_gwas_file) %>% summarize(best_coloc_h4 = max(clpp_h4), best_snp = ref_snp[which.max(clpp_h4)])

all_colocs_limited = full_join(sqtl_coloc_limited, eqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_eqtl"))
all_colocs_limited = full_join(all_colocs_limited, gene_edqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_gene_edqtl"))
all_colocs_limited = full_join(all_colocs_limited, cluster_edqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_cluster_edqtl"))
all_colocs_limited = full_join(all_colocs_limited, site_edqtl_coloc_limited, by = c("chr", "locus", "base_gwas_file"), suffix=c("", "_site_edqtl"))
colnames(all_colocs_limited)[which(colnames(all_colocs_limited) == "best_coloc_h4")] = "best_coloc_h4_sqtl"

# If not tested at all for that type, give it a 0
# to keep things simple.
all_colocs_limited[is.na(all_colocs_limited)] = 0

print(dim(all_colocs))
print(dim(all_colocs_limited))

# It turns out I guess in most (about 66%) cases we're not dealing with all that much of this doubling where we have multiple
# traits all colocalizing for the same trait

# Most values are either very close to 1 or very close to 0...
hist(all_colocs_limited$best_coloc_h4_sqtl)

write.table(all_colocs_limited, file="/users/mgloud/projects/rna_editing/output/qtl_comparisons/best_colocs_per_locus_limited.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Make the pairs plot
plot_matrix = all_colocs[,4:8]
colnames(plot_matrix) = gsub("best_coloc_h4_", "", colnames(plot_matrix))


pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_small.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.3)
dev.off()

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons.pdf", width=8, height=8)
	pairs(plot_matrix, pch = 16, cex=0.5)
dev.off()

plot_matrix = apply(plot_matrix, c(1,2), function(x) {max(x, 1e-5)})
plot_matrix = apply(plot_matrix, c(1,2), function(x) {min(x, 1 - 1e-5)})
logit_plot_matrix = logit(plot_matrix)

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_logit_small.pdf", width=8, height=8)
	pairs(logit_plot_matrix, pch = 16, cex=0.3)
dev.off()

pdf("/users/mgloud/projects/rna_editing/output/plots/qtl_comparisons/qtl_coloc_comparisons_logit.pdf", width=8, height=8)
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

immune_traits = c("Inflammatory-Bowel-Disease_Liu_2015_txt_gz")

### Make a stacked barplot ###

edqtl_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl > 0.5, best_coloc_h4_gene_edqtl > 0.5)
edqtl_no_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl <= 0.5, best_coloc_h4_gene_edqtl > 0.5)
no_edqtl_no_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl <= 0.5, best_coloc_h4_gene_edqtl <= 0.5)
no_edqtl_eqtl = all_colocs_limited %>% filter(best_coloc_h4_eqtl > 0.5, best_coloc_h4_gene_edqtl <= 0.5)

dim(edqtl_eqtl)[1]
dim(edqtl_no_eqtl)[1]
dim(no_edqtl_no_eqtl)[1]
dim(no_edqtl_eqtl)[1]

table(edqtl_no_eqtl$base_gwas_file)

edqtl_eqtl_percent = dim(edqtl_eqtl)[1] / dim(all_colocs_limited)[1]
edqtl_no_eqtl_percent = dim(edqtl_no_eqtl)[1] / dim(all_colocs_limited)[1]
no_edqtl_no_eqtl_percent = dim(no_edqtl_no_eqtl)[1] / dim(all_colocs_limited)[1]
no_edqtl_eqtl_percent = dim(no_edqtl_eqtl)[1] / dim(all_colocs_limited)[1]
all_counts = c(edqtl_eqtl_percent, edqtl_no_eqtl_percent, no_edqtl_no_eqtl_percent, no_edqtl_eqtl_percent)

immune_colocs_limited = all_colocs_limited[all_colocs_limited$base_gwas_file %in% immune_traits,]

immune_edqtl_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl > 0.5, best_coloc_h4_gene_edqtl > 0.5)
immune_edqtl_no_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl <= 0.5, best_coloc_h4_gene_edqtl > 0.5)
immune_no_edqtl_no_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl <= 0.5, best_coloc_h4_gene_edqtl <= 0.5)
immune_no_edqtl_eqtl = immune_colocs_limited %>% filter(best_coloc_h4_eqtl > 0.5, best_coloc_h4_gene_edqtl <= 0.5)

dim(immune_edqtl_eqtl)[1]
dim(immune_edqtl_no_eqtl)[1]
dim(immune_no_edqtl_no_eqtl)[1]
dim(immune_no_edqtl_eqtl)[1]

table(immune_edqtl_no_eqtl$base_gwas_file)

immune_edqtl_eqtl_percent = dim(immune_edqtl_eqtl)[1] / dim(immune_colocs_limited)[1]
immune_edqtl_no_eqtl_percent = dim(immune_edqtl_no_eqtl)[1] / dim(immune_colocs_limited)[1]
immune_no_edqtl_no_eqtl_percent = dim(immune_no_edqtl_no_eqtl)[1] / dim(immune_colocs_limited)[1]
immune_no_edqtl_eqtl_percent = dim(immune_no_edqtl_eqtl)[1] / dim(immune_colocs_limited)[1]
immune_counts = c(immune_edqtl_eqtl_percent, immune_edqtl_no_eqtl_percent, immune_no_edqtl_no_eqtl_percent, immune_no_edqtl_eqtl_percent)


subset = c(rep("all", 4), rep("autoimmune", 4))
condition = rep(c("both", "edQTL only", "eQTL only", "neither"), 2)
counts = c(all_counts, immune_counts)
stack_data = data.frame(subset, condition, counts)

ggplot(stack_data, aes(fill=condition, y=counts, x=subset)) +
	geom_bar(position="fill", stat="identity")
# TODO: I should really invert this thing to make it go from bottom up for colocs
