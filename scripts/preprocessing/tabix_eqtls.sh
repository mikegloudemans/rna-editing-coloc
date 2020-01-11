# Sort, zip, tabix edQTL files

# WARNING:
# WARNING:
# WARNING:
#
# If all of these are run at the same time, it will take a very long time
# to complete (a few days up to even a week or two). I recommend starting
# this ASAP, or first trying to obtain the already-tabixed files.
# 

# Tabix single-site editing QTLs
for f in `ls /users/mgloud/projects/rna_editing/data/eqtls/`; do
	echo -e "feature\tfeature_chr\tfeature_start\tfeature_end\tstrand\tN\tdistance\trsid\tchr\tsnp_pos\tsnp_pos_end\tpvalue\tbeta\tis_top\tref\talt" > /users/mgloud/projects/rna_editing/data/tabix_eqtls/$f.edQTLs.txt 
	zcat /users/mgloud/projects/rna_editing/data/eqtls/$f/$f.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.with.alleles.txt.gz | sed s/\ /\\t/g | sort -k9,9 -k10,10n >> /users/mgloud/projects/rna_editing/data/tabix_eqtls/$f.edQTLs.txt
	bgzip -f /users/mgloud/projects/rna_editing/data/tabix_eqtls/$f.edQTLs.txt
	tabix -f -S 1 -s 9 -b 10 -e 10 /users/mgloud/projects/rna_editing/data/tabix_eqtls/$f.edQTLs.txt.gz
done

# Do it for aggregated editing values too
for f in `ls /users/mgloud/projects/rna_editing/data/eqtls_aggro/*with.alleles*`; do
	filebase=`echo "$f" | sed s/.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Fisher_combined.with.alleles.txt.gz//g`
	filebase=`basename $filebase`
	echo -e "SNP\tgene\tchr\tsnp_pos\tmean_dist\tpvalue\tbeta\tref\talt" > /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	zcat /users/mgloud/projects/rna_editing/data/eqtls_aggro/$filebase.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Fisher_combined.with.alleles.txt.gz | tail -n +2 | sort -k3,3 -k4,4n >> /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	bgzip -f /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	tabix -f -S 1 -s 3 -b 4 -e 4 /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt.gz
done

# Do it for aggregated (cluster) editing values too
for f in `ls ../../data/edqtl_clusters/editingQLTs_agg_lancaster/*with.alleles*`; do
	filebase=`echo "$f" | sed s/.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Lancaster_agg.cluster.with.alleles.txt.gz//g`
	filebase=`basename $filebase`
	echo -e "SNP\tgene\tchr\tsnp_pos\tmean_dist\tpvalue\tbeta\tref\talt" > ../../data/edqtl_clusters/tabix/$filebase.ed.clust.txt
	zcat ../../data/edqtl_clusters/editingQLTs_agg_lancaster/$filebase.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Lancaster_agg.cluster.with.alleles.txt.gz | tail -n +2 | sort -k3,3 -k4,4n >> ../../data/edqtl_clusters/tabix/$filebase.ed.clust.txt
	bgzip -f ../../data/edqtl_clusters/tabix/$filebase.ed.clust.txt
	tabix -f -S 1 -s 3 -b 4 -e 4 ../../data/edqtl_clusters/tabix/$filebase.ed.clust.txt.gz
done

# While we're at it, let's do the v8 GTEx eQTLs too
for f in `ls ../../data/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations`; do
	echo -e "feature\tchr\tsnp_pos\tref\talt\tbuild\ttss_distance\tma_samples\tma_count\tmaf\tpvalue\tbeta\tse" > ../../data/shared_data/eqtls/gtex_v8/$f.eQTLs.txt
	zcat ../../data/gtex_v8/eqtl/GTEx_Analysis_v8_eQTL_all_associations/$f | sed s/_/\\t/g | sed s/chr//g | tail -n +2 | sort -k2,2 -k3,3n >> ../../data/shared_data/eqtls/gtex_v8/$f.eQTLs.txt
	bgzip -f ../../data/shared_data/eqtls/gtex_v8/$f.eQTLs.txt
	tabix -f -S 1 -s 2 -b 3 -e 3 ../../data/shared_data/eqtls/gtex_v8/$f.eQTLs.txt.gz
done

# And why not also do the v8 sQTLs, mostly the same thing
for f in `ls ../../data/gtex_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations`; do
	echo -e "feature\tfeature2\tchr\tsnp_pos\tref\talt\tbuild\ttss_distance\tma_samples\tma_count\tmaf\tpvalue\tbeta\tse" > ../../data/shared_data/sqtls/gtex_v8/$f.sQTLs.txt
	zcat ../../data/gtex_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations/$f | sed s/_/\\t/g | sed s/chr//g | tail -n +2 | sort -k3,3 -k4,4n >> ../../data/shared_data/sqtls/gtex_v8/$f.sQTLs.txt
	bgzip -f ../../data/shared_data/sqtls/gtex_v8/$f.sQTLs.txt
	tabix -f -S 1 -s 3 -b 4 -e 4 ../../data/shared_data/sqtls/gtex_v8/$f.sQTLs.txt.gz
done

