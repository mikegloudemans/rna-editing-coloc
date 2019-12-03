require(locuscomparer)

file1 = args[1]
file2 = args[2]
out_file = args[3]
title1 = args[4]
title2 = args[5]
vcf = args[6]
snp = args[7]

m = locuscompare(in_fn1 = file1, in_fn2 = file2, title=title1, title2=title2, snp=snp, genome="hg19") # vcf_fn=vcf, snp=snp)
m

# Plot with low quality to reduce file size; re-plot later if needing a high-quality figure

ggsave("/users/mgloud/projects/rna_editing/output/custom-ibd-plots/", width=12, height=6, dpi=50)
#ggsave(out_file, width=12, height=6, dpi=100)

