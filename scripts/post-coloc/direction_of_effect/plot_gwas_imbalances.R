require(ggplot2)

results = read.table("output/direction_of_effect/main_results_table.txt", header=TRUE)
results = results[(results$test_mode == "binomial-best-cluster") & (results$snp_set == "gwas_hits"),]

#results = results[(results$test_mode == "binomial-tissue-site-combined") & (results$snp_set == "gwas_hits"),]
#results = results[!grepl("Parkinsons", results$gwas_trait),]

code = c("Coronary-Artery-Disease_Nelson_2017",
	  "Inflammatory-Bowel-Disease",
          "Crohns-Disease",
          "Ulcerative-Colitis",
	"Lupus_Bentham_2015",
	"High-Density-Lipoprotein-GWAS-and-Metabochip",
	"Low-Density-Lipoprotein-GWAS-and-Metabochip",
	"Total-Cholesterol-GWAS-and-Metabochip",
	"Triglycerides-GWAS-and-Metabochip",
	"Multiple-Sclerosis",
	"Rheumatoid-Arthritis-European",
	"Amyotrophic-Lateral-Sclerosis-Meta-Analysis",
	"Psoriasis_Tsoi_2012",
	"Atopic-Dermatitis_EAGLE_2015",
	"Depression-Meta-Analysis",
	"Neuroticism_Luciano_2017",
	"Depressive-Symptoms_Okbay_2016",
	"Asthma-European-Fixed-Effect",
	"Celiac-Disease_Trynka_2011",
	"Parkinsons_Pankratz_2012",
	"Primary-Sclerosing-Cholangitis_Ji_2017",
	"Type-1-Diabetes-Meta-Analysis",
	"Alzheimers_Jansen_2018",
	"Schizophrenia-European",
	"BMI",
	"Height",
	"Primary-Biliary-Cirrhosis_Cordell_2015",
	"all"
)


#results$short_names = c("Asthma", "Psoriasis", "Lupus", "Celiac", "Eczema", "CAD", "MS", "IBD", "Crohns", "UC", "Allergies")
results$short_names = code
results = results[results$num_negative + results$num_positive > 5,]

results = results[c(order(results$short_names[1:length(results$short_names)-1]), length(results$short_names)),]

#results = rbind(results,data.frame(list(gwas_trait="ALL", test_mode="binomial", snp_set="gwas_hits", num_positive=sum(results$num_positive), num_negative = sum(results$num_negative), num_na = sum(results$num_na), short_names="All combined")))
results$short_names = factor(results$short_names, levels=rev(results$short_names))

results$frac_positive = results$num_positive / (results$num_negative + results$num_positive)

# Temporary hack
results[results$short_names == "Atopic-Dermatitis_EAGLE_2015",]$frac_positive = 1-results[results$short_names == "Atopic-Dermatitis_EAGLE_2015",]$frac_positive
results = results[!results$short_names == "all",]

# Get width of confidence intervals for the binomial test for imbalance,
# for each trait
results$ci_width = qnorm(0.975)*sqrt((results$frac_positive) * (1-results$frac_positive) / (results$num_negative + results$num_positive))
results$point_shape = 19
#results$point_shape[length(results$point_shape)] = 25
results$point_shape = factor(results$point_shape)
results$bar_color = "grey60"
#results$bar_color[length(results$bar_color)] = "black"
results$size=3
#results$size[length(results$size)] = 4
results$size=factor(results$size)

# Thanks for the code, Nicole Ferraro
fig = ggplot(results, aes(x=short_names,y=frac_positive)) +
	coord_flip() +
	geom_hline(yintercept=0.5, linetype="dotted") +
	geom_point(aes(shape=point_shape, size=size, color=bar_color), position=position_dodge(width=0)) +
	scale_shape_manual(values=c(19,18)) +
	scale_color_manual(values=c('black', 'grey50')) +
	scale_size_manual(values=c(4,6)) +
	geom_errorbar(aes(ymin = pmax(0, frac_positive - ci_width), ymax = pmin(1, frac_positive + ci_width), col=bar_color), width=0, lwd=1, position=position_dodge(width=0)) +
	theme_bw() + ylab('Editing imbalance of GWAS variant') + xlab('') +
	theme(axis.text.y = element_text(size=11), legend.position = "none")
fig

ggsave("gwas_plot.pdf", width=4, height=4)
