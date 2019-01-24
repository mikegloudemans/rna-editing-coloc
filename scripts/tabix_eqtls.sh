# Sort, zip, tabix edQTL files

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
	zcat $f | head -n 1 > /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	zcat /users/mgloud/projects/rna_editing/data/eqtls_aggro/$filebase.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Fisher_combined.with.alleles.txt.gz | tail -n +2 | sort -k3,3 -k4,4n >> /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	bgzip -f /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt
	tabix -f -S 1 -s 3 -b 4 -e 4 /users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/$filebase.Fisher.combined.sorted.txt.gz
done

