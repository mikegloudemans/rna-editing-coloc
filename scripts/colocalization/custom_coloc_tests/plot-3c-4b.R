require(ggplot2)
require(cowplot)
require(locuscomparer)

# Some added stuff that allows me to modify axes in locuscomparer
source("lc_mod.R")

# TODO: Figure out a way to make this visible for the public version without
# having to modify the code to keep tmp files.

# Probably the way to go is to just provide these files straight out

# Should also be separate functions

# And should not create a "trunc.eqtl" file in the scripts directory

edqtl = "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_Fisher_combined_sorted_txt_gz/TNFRSF14_eqtl_lc_data.txt"
edqtl_gwas =  "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_Fisher_combined_sorted_txt_gz/TNFRSF14_gwas_lc_data.txt"

edqtl_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Spleen edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38")
m

ggsave("../output/colocalization/custom_coloc_tests/figure-4b/4b_edQTL.pdf", width=12, height=6)

edqtl = "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000157873.17_eqtl_lc_data.txt"
edqtl_gwas = "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000157873.17_gwas_lc_data.txt"

edqtl_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, 10, 10, title="Transverse Colon edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Spleen eQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m

ggsave("../output/colocalization/custom_coloc_tests/figure-4b/4b_eQTL.pdf", width=12, height=6)


edqtl = "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_sQTLs_txt_gz/1.2558468.2559823.clu_eqtl_lc_data.txt"
edqtl_gwas = "/users/mgloud/projects/brain_gwas/tmp/2019-11-15_22-41-02.058618/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Inflammatory-Bowel-Disease/1_2540642/Spleen_sQTLs_txt_gz/1.2558468.2559823.clu_gwas_lc_data.txt"

edqtl_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, 10, 10, title="Transverse Colon edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Spleen sQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m

ggsave("../output/colocalization/custom_coloc_tests/figure-4b/4b_sQTL.pdf", width=12, height=6)

system("echo \"rsid\\tpval\" > trunc.eqtl")
system("tail -n +4000 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon_Transverse_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000249669.9_eqtl_lc_data.txt >> trunc.eqtl")
system("echo \"rsid\\tpval\" > trunc.eqtl.gwas")
system("tail -n +4000 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon_Transverse_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000249669.9_gwas_lc_data.txt >> trunc.eqtl.gwas")

edqtl = "trunc.eqtl"
edqtl_gwas = "trunc.eqtl.gwas"
edqtl_table = read.table(edqtl_gwas)
g_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, 10, 10, title="Transverse Colon edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Transverse Colon eQTL, CARMN", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m

ggsave("../output/colocalization/custom_coloc_tests/figure-3c/3c_eQTL.pdf", width=12, height=6)


system("echo \"rsid\\tpval\" > trunc.eqtl")
system("tail -n +2500 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon-Transverse_Fisher_combined_sorted_txt_gz/CARMN_eqtl_lc_data.txt >> trunc.eqtl")
system("echo \"rsid\\tpval\" > trunc.eqtl.gwas")
system("tail -n +2500 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon-Transverse_Fisher_combined_sorted_txt_gz/CARMN_gwas_lc_data.txt >> trunc.eqtl.gwas")

edqtl = "trunc.eqtl"
edqtl_gwas = "trunc.eqtl.gwas"
edqtl_table = read.table(edqtl_gwas)
g_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, 10, 10, title="Transverse Colon edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Transverse Colon edQTL, CARMN", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m

ggsave("../output/colocalization/custom_coloc_tests/figure-3c/3c_edQTL.pdf", width=12, height=6)


system("echo \"rsid\\tpval\" > trunc.eqtl")
system("tail -n +3500 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon_Transverse_sQTLs_txt_gz/5.149429039.149430915.clu_eqtl_lc_data.txt >> trunc.eqtl")
system("echo \"rsid\\tpval\" > trunc.eqtl.gwas")
system("tail -n +3500 /users/mgloud/projects/brain_gwas/tmp/2019-11-16_00-12-30.543891/locuscompare/Inflammatory-Bowel-Disease_Liu_2015_txt_gz/Ulcerative-Colitis/5_150236753/Colon_Transverse_sQTLs_txt_gz/5.149429039.149430915.clu_gwas_lc_data.txt >> trunc.eqtl.gwas")

edqtl = "trunc.eqtl"
edqtl_gwas = "trunc.eqtl.gwas"
edqtl_table = read.table(edqtl_gwas)
g_table = read.table(edqtl_gwas)
snp = edqtl_table$V1[which.min(as.numeric(as.character(edqtl_table$V2), na.rm=TRUE))]

#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, 10, 10, title="Transverse Colon edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Transverse Colon sQTL, CARMN", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38") # vcf_fn=vcf, snp=snp)
m

ggsave("../output/colocalization/custom_coloc_tests/figure-3c/3c_sQTL.pdf", width=12, height=6)
