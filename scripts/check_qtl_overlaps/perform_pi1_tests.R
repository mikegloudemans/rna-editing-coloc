require(dplyr)
require(readr)
require(qvalue)
require(metap)
require(ggplot2)
require(reshape2)
require(RColorBrewer)

# TODO: To be sure that this pipeline is SOMEWHAT helping, pipe through one set of p-values that's uniformly random
# and maybe one that's highly correlated, to see what we get there instead...show this in slides

fishers = function(x)
{
	if (length(x) < 2)
	{
		return(min(x))
	}
	else
	{
		return(sumlog(x)$p)
	}
}

# TODO: remove comment from here
full_data = read_delim("/users/mgloud/projects/rna_editing/output/qtl_comparisons.txt", delim="\t")

possible_sources = c("egenes", "sgenes", "edFeats")
possible_targets = c("eQTL", "sQTL", "edQTL", "aggro")

for (ps in possible_sources)
{
	pt_trunc = possible_targets[-which(possible_sources == ps)]

	data = full_data[grepl(ps, full_data$source_file),]

	overlaps = array(0, dim=c(length(unique(data$source_file)), 5))
	rownames(overlaps) = sapply(unique(data$source_file), function(x) 
				    {
					    s = strsplit(x, "/")[[1]]
					    return(s[length(s)])
				    })
	rownames(overlaps) = gsub("\\.v8\\.egenes\\.txt\\.gz", "", rownames(overlaps))

	colnames(overlaps) = c(pt_trunc, "all", "null_control")

	for (sf in 1:length(unique(data$source_file)))
	{
		src_file = unique(data$source_file)[sf]

		#hist(data$target_pval)
		#pi1 = 1 - qvalue(data$target_pval)$pi0 ## pi1 if we don't consider sharing across tissues
		
		tests = data[data$source_file == src_file,]
		tests = tests[tests$target_pval != 0,]

		#hist(tests$target_pval)
		#pi1 = 1 - qvalue(tests$target_pval)$pi0 ## pi1 if we don't consider sharing across tissues

		# These two tests suffer from the same problem -- inflation of the large values...
		#meta = tests %>% group_by(target_file, source_feature) %>% summarize(fishers_p = fishers(target_pval))
		#meta = tests %>% group_by(target_file, source_feature) %>% summarize(meta_pval = min(p.adjust(target_pval, method="BH")))
		#meta = tests %>% group_by(target_file, source_feature) %>% summarize(bh_pval = min(p.adjust(target_pval, method="bonferroni")))

		# TODO: Figure out a better way to perform this analysis!
		# Cut the meta pvals into bins
		# For now, we'll throw away the top bin and then do a very rough 
		# geometric estimate of the amount of significant p-values
		#bins = cut(meta$meta_pval, 20)
		#non_sig = 20 * mean(table_bins[1:18])

		#hist(meta$meta_pval)
		#pi1 = 1 - qvalue(meta$meta_pval)$pi0 ## pi1 if we don't consider sharing across tissues

		meta = tests %>% group_by(target_file, source_feature) %>% summarize(significant = min(p.adjust(target_pval, method="BH")) < 0.05)

		# Get overlap with each type of target p-val
		for (target in unique(meta$target_file))
		{
			# Get the percentage of targets in this source / target pair that led to a BH corrected p-value of < 0.05
			print(src_file)
			print(target)
			mini = meta[meta$target_file == target,]
			print(sum(mini$significant) / length(mini$significant))
			print("")

			target_type = -1
			if (grepl(pt_trunc[1], target))
			{
				target_type = 1	
			}
			else if (grepl(pt_trunc[2], target))
			{
				target_type = 2
			}
			else if (grepl(pt_trunc[3], target))
			{
				target_type = 3
			}

			overlaps[sf, target_type] = sum(mini$significant) / length(mini$significant)
		}

		# Now combine all to get one metascore saying whether the eQTL's duplicated in any other eQTL type
		#sub = overlaps[sf, 1:3]
		#sub = sub[sub != 0]
		#overlaps[sf, 4] = min(p.adjust(sub, method="BH"))

		# Just collapse across all target studies. This might not be exactly "fair" because some studies have a lot more potential targets than others
		# TODO: Figure out what exactly to do about this to make it more fair
		any_test = tests %>% group_by(source_feature) %>% summarize(significant = min(p.adjust(target_pval, method="BH")) < 0.05)
		overlaps[sf, 4] = sum(any_test$significant) / length(any_test$significant)

		null_tests = tests
		null_tests$target_pval = runif(length(null_tests$target_pval))		
		null_summary = null_tests %>% group_by(source_feature) %>% summarize(significant = min(p.adjust(target_pval, method="BH")) < 0.05)
		overlaps[sf, 5] = sum(null_summary$significant) / length(null_summary$significant) 
	}

	print(overlaps)


	# Plot the replication results in a bar graph
	melted = melt(overlaps)
	colnames(melted) = c("qtl_tissue", "target_type", "percent_replicated")	
	g = ggplot(data=melted, aes(x=qtl_tissue, y=percent_replicated, fill=target_type)) +
		geom_bar(stat="identity", position=position_dodge()) +
		theme_bw() +
		theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_y_continuous(limits = c(0, 1)) +
		xlab(paste(ps, "QTL source tissue")) +
		ylab("Percent replication") +
		scale_fill_brewer(palette="Set1", name="Replication\nQTL type")
	g

	ggsave(paste0("/users/mgloud/projects/rna_editing/output/plots/check_qtl_overlaps/", ps, "_replication.png"), width=18, height=6)

}



# A little simulation; not relevant to the above code...
#fish_test = rep(0, 10000)
#for (i in 1:10000)
#{
#	r = c(rep(runif(1), 20)	, rep(runif(1), 30))
	#r = runif(50
#	fish_test[i] = (sumlog(r)$p)
#}

