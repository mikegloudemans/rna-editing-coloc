require(tidyverse)
require(ggplot2)

coloc_threshold = 0.5

# Load all loci, coloc status
loci = read.csv("/users/mgloud/projects/rna_editing/data/HuangIBDAnnotationSummary.csv", header=TRUE)
coloc_summary = read.table("/users/mgloud/projects/rna_editing/output/ibd_analysis/ibd_coloc_summary.tsv", header=TRUE, fill=TRUE, sep="\t")
colocs = coloc_summary$HD[coloc_summary$best_edqtl_coloc > coloc_threshold]
coloc_genes = coloc_summary$best_edqtl_feature[coloc_summary$best_edqtl_coloc > coloc_threshold]

# Sort by p-value
loci = loci[rev(order(loci$p_multi)),]

# Get our annotations for each locus
coding_loci = unique(loci[loci$Coding != "",]$HD)

loci$has_coding = loci$HD %in% coding_loci
loci$has_coloc = loci$HD %in% colocs

loci$coding_genes = ""
loci$coding_genes[loci$has_coding] = sapply(loci[loci$has_coding,]$HD, function(x) 
		      {
			      idx = which((loci$HD==x) & (loci$Coding!=""))[1]
			      return(strsplit(as.character(loci$Coding[idx]), "\\(")[[1]][1])
		      })

loci$edqtl_genes = ""
loci$edqtl_genes[loci$has_coloc] = sapply(loci[loci$has_coloc,]$HD, function(x) 
		      {
			      return(as.character(coloc_summary$best_edqtl_feature[which(coloc_summary$HD==x)]))
		      })

coloc_results = read.table("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", header=TRUE, stringsAsFactors=FALSE)
coloc_results = coloc_results[coloc_results$base_gwas_file == "Inflammatory-Bowel-Disease_Liu_2015_txt_gz",]

relevant_tissues = c("Adipose-Subcutaneous",
		     "Pancreas",
		     "Esophagus-Mucosa",
		     "Adipose-Visceral-Omentum",
		     "Colon-Transverse",
		     "Esophagus-Muscularis",
		     "Thyroid",
		     "Brain-Hippocampus",
		     "Whole-Blood",
		     "Cells-EBV-transformed-lymphocytes",
		     "Colon-Sigmoid"
		     )

loci$coloc_in_relevant_tissue = sapply(loci$edqtl_gene, function(x)
	{
		sub = coloc_results[coloc_results$feature == x,]
		sub = sub[grepl("Fisher_combined", sub$eqtl_file),]
		sub = sub[sub$clpp_h4 > 0.5,]
		sub$eqtl_file = gsub("_Fisher_combined_sorted_txt_gz", "", sub$eqtl_file)
		return(sum(sub$eqtl_file %in% relevant_tissues) > 0)

	})

# Remove duplicates now (since they're already sorted by p-value)
loci = loci[!duplicated(loci$HD),]

loci$annotation = ""
loci[loci$has_coding & loci$has_coloc,]$annotation = "both"
loci[loci$has_coding & !loci$has_coloc,]$annotation = "coding"
loci[!loci$has_coding & loci$has_coloc,]$annotation = "edQTL coloc"
loci[!loci$has_coding & !loci$has_coloc,]$annotation = "neither"
loci$annotation = factor(loci$annotation, levels=c("edQTL coloc", "coding", "both", "neither"))

loci$significance = "yes"

our_loci = loci[loci$annotation != "neither",]

# Plot as a colored heatmap stack
g = ggplot(loci,aes(x=significance,y=factor(p_multi),fill=annotation)) +
	theme_classic() +
	theme(legend.position="top", legend.title = element_blank()) +
	scale_fill_manual(values = c("red", "orange2", "yellow3", "gray20")) +
	theme(axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank()) +
	theme(axis.title.y=element_blank(),
	      axis.text.y=element_blank(),
	      axis.ticks.y=element_blank()) +
	theme(plot.margin = unit(c(0, 7, 0, 7), "cm")) + 
	geom_tile()
g


# Add annotation labels on top


barplot(pmin(loci$p_multi, 50))


g = ggplot(loci,aes(x=factor(1:length(p_multi)),y=pmin(p_multi, 50),fill=annotation)) +
	theme_classic() +
	theme(legend.position="top", legend.title = element_blank()) +
	scale_fill_manual(values = c("red", "orange2", "yellow3", "gray20")) +
	theme(axis.text.x=element_blank()) +
	theme(axis.ticks.y=element_blank()) +
	ylab("-log10(P-value)") +
	xlab("IBD GWAS Loci") +
	geom_bar(stat="identity")

g

g = ggplot(our_loci,aes(x=factor(p_multi, levels=rev(p_multi)),y=pmin(p_multi, 50),fill=annotation)) +
	theme_classic() +
	theme(legend.position="top", legend.title = element_blank()) +
	scale_fill_manual(values = c("forestgreen", "orange2", "red", "gray20")) +
	theme(axis.text.y=element_blank()) +
	ylab("-log10(GWAS P-value)") +
	xlab("IBD GWAS Loci") +
	geom_bar(stat="identity") +
       	geom_text(aes(label=edqtl_genes, y=55, hjust=0, fontface=ifelse(coloc_in_relevant_tissue,"bold","plain"))) + 
       	#geom_text(aes(label=coding_genes, y=-5, hjust=1)) + 
	ylim(0,80) +
	coord_flip()

g

ggsave("panel3b.pdf", height=6, width=6)

