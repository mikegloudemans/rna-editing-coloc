require(tidyverse)

# Load and clean data
gwas_hits = read.csv("data/ibd-huang-annotations/Huang2017GWAS-IBD.csv", header=TRUE)
gwas_hits$credible_start = as.numeric(gsub(",", "", as.character(gwas_hits$credible_start)))
gwas_hits$credible_end = as.numeric(gsub(",", "", as.character(gwas_hits$credible_end)))+1 # Add one to end to make it an UCSC-style interval
gwas_hits$chr = paste0("chr", gwas_hits$chr)

bed_only = gwas_hits[c("chr", "credible_start", "credible_end")]

# First we need to liftOver all the loci to hg38
write.table(bed_only, file="tmp/bed_for_liftover.tmp", sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)
system("liftOver tmp/bed_for_liftover.tmp data/liftOver/chains/hg19ToHg38.over.chain.gz output/ibd-huang-annotations/liftover/Huang2017GWAS-IBD-lifted.bed output/ibd-huang-annotations/liftover/Huang2017GWAS-IBD-lifted-liftover-failed.bed")

lifted_gwas = read.table("output/ibd-huang-annotations/liftover/Huang2017GWAS-IBD-lifted.bed", header=FALSE)
# Check below, but the following command hopefully won't even run because none fail
#lifted_failed = read.table("/users/mgloud/projects/rna_editing/output/liftover/Huang2017GWAS-IBD-lifted-liftover-failed.bed", header=FALSE)
colnames(lifted_gwas) = c("chr_hg38", "credible_start_hg38", "credible_end_hg38")
lifted_hits = cbind(gwas_hits, lifted_gwas)

# Select the first point from each GWAS hit as the representative

# Then we need to load all coloc results for IBD only
coloc_results = read.table("output/colocalization/rna-editing-revisions/concatenated/all_coloc_results.txt", header=TRUE, stringsAsFactors=FALSE)
#ibd_coloc_only = coloc_results[grep("Bowel", coloc_results$gwas_trait),]
ibd_coloc_only = coloc_results[grep("Bowel", coloc_results$base_gwas_file),]

# Determine which IBD locus, if any, each coloc result belongs to
ibd_coloc_only$HD = sapply(ibd_coloc_only$ref_snp, function(x)
			   {
				   chrom = paste0("chr", strsplit(x, "_")[[1]][1])
				   pos = as.numeric(strsplit(x, "_")[[1]][2])
				   sub = lifted_hits[lifted_hits$chr_hg38 == chrom,]
				   sub$dist = abs(((sub$credible_start + sub$credible_end) / 2) - pos)
				   sub = sub[sub$dist < 500000,]
				   if (dim(sub)[1] == 0)
				   {
					   return(NA)
				   }
				   sub = sub[order(sub$dist),]
				   return(sub$HD[1])
			   })
ibd_coloc_only = ibd_coloc_only[!is.na(ibd_coloc_only$HD),]

# Percent of loci where we have at least tested for colocalization
print(length(unique(ibd_coloc_only$HD)) / length(unique(lifted_hits$HD)))

coloc_threshold = 0.5
coloc_threshold = 0.9

# Report top tissue, top colocalization, binary colocalization status for each coloc type.
# Also report top QTL feature within each category.
ibd_coloc_summary = ibd_coloc_only %>% group_by(HD) %>% summarize(best_eqtl_coloc = max(clpp_h4[grepl("eQTL", eqtl_file)]), 
					      best_eqtl_tissue = gsub("_allpairs_txt_gz_eQTLs_txt_gz", "", c("",eqtl_file)[c(which((clpp_h4 == best_eqtl_coloc) & (grepl("eQTL", eqtl_file)))+1,1)[1]]),
					      best_eqtl_feature = c("",feature)[c(which((clpp_h4 == best_eqtl_coloc) & (grepl("eQTL", eqtl_file)))+1,1)[1]],
					      best_eqtl_ref_snp = c("",ref_snp)[c(which((clpp_h4 == best_eqtl_coloc) & (grepl("eQTL", eqtl_file)))+1,1)[1]],
					      best_eqtl_gwas_trait = c("",gwas_trait)[c(which((clpp_h4 == best_eqtl_coloc) & (grepl("eQTL", eqtl_file)))+1,1)[1]],
					      eqtl_coloc_status = best_eqtl_coloc > coloc_threshold,

					      best_sqtl_coloc = max(clpp_h4[grepl("sQTL", eqtl_file)]), 
					      best_sqtl_tissue = gsub("_sQTLs_txt_gz", "", c("",eqtl_file)[c(which((clpp_h4 == best_sqtl_coloc) & (grepl("sQTL", eqtl_file)))+1,1)[1]]),
					      best_sqtl_feature = c("",feature)[c(which((clpp_h4 == best_sqtl_coloc) & (grepl("sQTL", eqtl_file)))+1,1)[1]],
					      best_sqtl_ref_snp = c("",ref_snp)[c(which((clpp_h4 == best_sqtl_coloc) & (grepl("sQTL", eqtl_file)))+1,1)[1]],
					      best_sqtl_gwas_trait = c("",gwas_trait)[c(which((clpp_h4 == best_sqtl_coloc) & (grepl("sQTL", eqtl_file)))+1,1)[1]],
					      sqtl_coloc_status = best_sqtl_coloc > coloc_threshold,
					      best_edqtl_coloc = max(clpp_h4[grepl("Fisher_combined", eqtl_file)]), 
					      best_edqtl_tissue = gsub("_Fisher_combined_sorted_txt_gz", "", c("",eqtl_file)[c(which((clpp_h4 == best_edqtl_coloc) & (grepl("Fisher_combined", eqtl_file)))+1,1)[1]]),
					      best_edqtl_feature = c("",feature)[c(which((clpp_h4 == best_edqtl_coloc) & (grepl("Fisher_combined", eqtl_file)))+1,1)[1]],
					      best_edqtl_ref_snp = c("",ref_snp)[c(which((clpp_h4 == best_edqtl_coloc) & (grepl("Fisher_combined", eqtl_file)))+1,1)[1]],
					      best_edqtl_gwas_trait = c("",gwas_trait)[c(which((clpp_h4 == best_edqtl_coloc) & (grepl("Fisher_combined", eqtl_file)))+1,1)[1]],
					      edqtl_coloc_status = best_edqtl_coloc > coloc_threshold				      
					      )
ibd_coloc_summary$best_eqtl_feature = sapply(ibd_coloc_summary$best_eqtl_feature, function(x) {strsplit(x, "\\.")[[1]][1]})
hgnc = read.table("data/hgnc/hgnc_biomart.txt", header=TRUE, fill=TRUE)
hgnc = hgnc[c("Gene","ID")]
colnames(hgnc) = c("best_eqtl_feature", "best_eqtl_feature_hgnc")
ibd_coloc_summary = left_join(ibd_coloc_summary, hgnc)

# Get hg38 coordinates for each locus
lifted_map = lifted_hits[!duplicated(lifted_hits$HD),]
lifted_map = lifted_map[c("HD", "chr_hg38", "credible_start_hg38", "credible_end_hg38")]

ibd_coloc_summary = left_join(ibd_coloc_summary, lifted_map)


ibd_coloc_summary$best_eqtl_feature_hgnc[is.na(ibd_coloc_summary$best_eqtl_feature_hgnc)] = ""
ibd_coloc_summary = ibd_coloc_summary[order(as.character(ibd_coloc_summary$best_eqtl_feature_hgnc), decreasing=TRUE),]
ibd_coloc_summary = ibd_coloc_summary[!duplicated(ibd_coloc_summary$HD),]
ibd_coloc_summary = ibd_coloc_summary[order(ibd_coloc_summary$HD),]

print(dim(ibd_coloc_summary))
print(colSums(ibd_coloc_summary[c("eqtl_coloc_status", "sqtl_coloc_status", "edqtl_coloc_status")]))
print(colMeans(ibd_coloc_summary[c("eqtl_coloc_status", "sqtl_coloc_status", "edqtl_coloc_status")]))

# Get total % of loci colocalized in each QTL category 
sum(!(ibd_coloc_summary$eqtl_coloc_status | ibd_coloc_summary$sqtl_coloc_status) & ibd_coloc_summary$edqtl_coloc_status)
sum(!(ibd_coloc_summary$edqtl_coloc_status | ibd_coloc_summary$sqtl_coloc_status) & ibd_coloc_summary$eqtl_coloc_status)
sum(!(ibd_coloc_summary$edqtl_coloc_status | ibd_coloc_summary$eqtl_coloc_status) & ibd_coloc_summary$sqtl_coloc_status)

# Or what percent colocalized with any category at all
sum(ibd_coloc_summary$edqtl_coloc_status | ibd_coloc_summary$eqtl_coloc_status | ibd_coloc_summary$sqtl_coloc_status)

# Rank QTL types from top to bottom for each locus
ibd_coloc_summary[c("eqtl_coloc_rank","sqtl_coloc_rank","edqtl_coloc_rank")] = t(apply(ibd_coloc_summary[c("best_eqtl_coloc", "best_sqtl_coloc", "best_edqtl_coloc")], 1, function(x) {rank(-x)}))

ibd_coloc_summary = ibd_coloc_summary[c("HD", "chr_hg38", "credible_start_hg38", "credible_end_hg38", "best_edqtl_coloc", "best_edqtl_tissue", "best_edqtl_feature", "best_edqtl_ref_snp", "best_edqtl_gwas_trait", "edqtl_coloc_status", "best_sqtl_coloc", "best_sqtl_tissue", "best_sqtl_feature", "best_sqtl_ref_snp", "best_sqtl_gwas_trait", "sqtl_coloc_status", "best_eqtl_coloc", "best_eqtl_tissue", "best_eqtl_feature", "best_eqtl_ref_snp", "best_eqtl_gwas_trait", "eqtl_coloc_status", "best_eqtl_feature_hgnc", "edqtl_coloc_rank","sqtl_coloc_rank","eqtl_coloc_rank")]


write.table(ibd_coloc_summary, "output/ibd-huang-annotations/ibd_coloc_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Check on some of the GWAS loci that we didn't see any tests for...
lifted_hits[!(lifted_hits$HD %in% ibd_coloc_only$HD),][c("HD", "chr_hg38", "credible_start_hg38", "credible_end_hg38")]
