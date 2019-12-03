require(tidyverse)
require(ggplot2)

gwas_hits = read.csv("/users/mgloud/projects/rna_editing/data/HuangIBDAnnotationSummary.csv", header=TRUE, fill=TRUE)
coloc_summary = read.table("/users/mgloud/projects/rna_editing/output/ibd_analysis/ibd_coloc_summary.tsv", header=TRUE, fill=TRUE, sep="\t")

# Test whether a locus is fully annotated, ambiguous, or unannotated
is_annotated = function(x)
{
	if (sum(x) == 0)
	{
		return("un-annotated")
	}
	else if (sum(!x) == 0)
	{
		return("annotated")
	}
	else
	{
		return("ambiguous")
	}
}

# Make custom annotations to illustrate what level of annotation we see

# Has any annotation other than eQTL
gwas_hits$has_noneqtl_annotation = (gwas_hits$Coding != "") | (gwas_hits$Immune_Cell != "") | (gwas_hits$Metabolite != "") | (gwas_hits$epigenetic != "") | (gwas_hits$epigenetic_all != "") | (gwas_hits$TFBS != "")

# Has any annotation that's not simply an epigenetic mark
gwas_hits$has_nonepigenetic_annotation = (gwas_hits$Coding != "") | (gwas_hits$Immune_Cell != "") | (gwas_hits$Metabolite != "") | (gwas_hits$eQTL != "")

# Has any coding annotation
gwas_hits$has_coding = (gwas_hits$Coding != "") 

# Has any annotation that's not eQTL or epigenetic annotation
gwas_hits$has_nonepigenetic_noneqtl_annotation = (gwas_hits$Coding != "") | (gwas_hits$Immune_Cell != "") | (gwas_hits$Metabolite != "")


# Now aggregate across independent signals at each locus
annotation_summary = gwas_hits %>% group_by(HD) %>% summarize(has_annotation = is_annotated(Annotated), has_coding = is_annotated(has_coding), has_noneqtl_annotation = is_annotated(has_noneqtl_annotation), has_nonepigenetic_annotation = is_annotated(has_nonepigenetic_annotation), has_nonepigenetic_noneqtl_annotation = is_annotated(has_nonepigenetic_noneqtl_annotation))

has_new_edqtl_coloc = function(x)
{
	if (!(x %in% coloc_summary$HD))
	{
		return(FALSE)
	}
	else
	{
		return(coloc_summary$edqtl_coloc_status[which(coloc_summary$HD == x)])
	}
}

has_new_eqtl_coloc = function(x)
{
	if (!(x %in% coloc_summary$HD))
	{
		return(FALSE)
	}
	else
	{
		return(coloc_summary$eqtl_coloc_status[which(coloc_summary$HD == x)])
	}
}

has_new_sqtl_coloc = function(x)
{
	if (!(x %in% coloc_summary$HD))
	{
		return(FALSE)
	}
	else
	{
		return(coloc_summary$sqtl_coloc_status[which(coloc_summary$HD == x)])
	}
}



# Add colocalization information
annotation_summary$has_new_edqtl_coloc = sapply(annotation_summary$HD, has_new_edqtl_coloc)
annotation_summary$has_new_eqtl_coloc = sapply(annotation_summary$HD, has_new_eqtl_coloc)
annotation_summary$has_new_sqtl_coloc = sapply(annotation_summary$HD, has_new_sqtl_coloc)

annotation_summary$has_no_coloc = !annotation_summary$has_new_edqtl_coloc & (!annotation_summary$has_new_eqtl_coloc) & (!annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_edqtl_coloc_only = annotation_summary$has_new_edqtl_coloc & (!annotation_summary$has_new_eqtl_coloc) & (!annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_sqtl_coloc_only = !annotation_summary$has_new_edqtl_coloc & (!annotation_summary$has_new_eqtl_coloc) & (annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_eqtl_coloc_only = !annotation_summary$has_new_edqtl_coloc & (annotation_summary$has_new_eqtl_coloc) & (!annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_edqtl_sqtl_coloc_only = annotation_summary$has_new_edqtl_coloc & (!annotation_summary$has_new_eqtl_coloc) & (annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_edqtl_eqtl_coloc_only = annotation_summary$has_new_edqtl_coloc & (annotation_summary$has_new_eqtl_coloc) & (!annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_sqtl_eqtl_coloc_only = !annotation_summary$has_new_edqtl_coloc & (annotation_summary$has_new_eqtl_coloc) & (annotation_summary$has_new_sqtl_coloc)
annotation_summary$has_new_sqtl_eqtl_edqtl_coloc = annotation_summary$has_new_edqtl_coloc & (annotation_summary$has_new_eqtl_coloc) & (annotation_summary$has_new_sqtl_coloc)


print(table(annotation_summary[c("has_annotation", "has_new_edqtl_coloc")]))
print(table(annotation_summary[c("has_coding", "has_new_edqtl_coloc")]))
print(table(annotation_summary[c("has_noneqtl_annotation", "has_new_edqtl_coloc")]))
print(table(annotation_summary[c("has_nonepigenetic_annotation", "has_new_edqtl_coloc")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_edqtl_coloc")]))

print(table(annotation_summary[c("has_annotation", "has_new_eqtl_coloc")]))
print(table(annotation_summary[c("has_coding", "has_new_eqtl_coloc")]))
print(table(annotation_summary[c("has_noneqtl_annotation", "has_new_eqtl_coloc")]))
print(table(annotation_summary[c("has_nonepigenetic_annotation", "has_new_eqtl_coloc")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_eqtl_coloc")]))

print(table(annotation_summary[c("has_annotation", "has_new_edqtl_coloc_only")]))
print(table(annotation_summary[c("has_coding", "has_new_edqtl_coloc_only")]))
print(table(annotation_summary[c("has_noneqtl_annotation", "has_new_edqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_annotation", "has_new_edqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_edqtl_coloc_only")]))

print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_edqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_sqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_eqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_edqtl_sqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_edqtl_eqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_sqtl_eqtl_coloc_only")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_new_sqtl_eqtl_edqtl_coloc")]))
print(table(annotation_summary[c("has_nonepigenetic_noneqtl_annotation", "has_no_coloc")]))

# Turn it into a plot

# Totally a template example here, unrelated
#
#subset = c(rep("all", 4), rep("autoimmune", 4))
#condition = rep(c("both", "edQTL only", "eQTL only", "neither"), 2)
#counts = c(all_counts, immune_counts)
#stack_data = data.frame(subset, condition, counts)
#
#ggplot(stack_data, aes(fill=condition, y=counts, x=subset)) +
#	geom_bar(position="fill", stat="identity")


g = ggplot(annotation_summary, aes(has_annotation))
g + geom_bar()
g + geom_bar(aes(fill = has_new_edqtl_coloc))


g = ggplot(annotation_summary, aes(has_nonepigenetic_noneqtl_annotation))
g + geom_bar()
g + geom_bar(aes(fill = has_new_edqtl_coloc))

