#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 6/7/2018

# Print key attributes about GWAS, that we can use to select best GWAS
# (May also be useful for others in lab to make this table and post)

import glob
import gzip
import subprocess
import sys
import operator
import pandas as pd
sys.path.insert(0, '/users/mgloud/projects/brain_gwas/scripts')
import SNP 

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


def main():
    files = glob.glob("/users/mgloud/projects/rna_editing/data/gwas/*.gz")

    with open("/users/mgloud/projects/rna_editing/output/snps_to_test_aggro_gwas1e6_eqtl1e6.txt", "w") as w:

        w.write("chr\tsnp_pos\tgwas_file\teqtl_file\ttrait\tgwas_pvalue\teqtl_pvalue\tgene\n")
        w.flush()

        for file in sorted(files):

            with gzip.open(file) as f:
                if not "effect_allele" in f.readline().strip().split():
                    continue
            
            info = snps_by_threshold(file, 1e-6, file)

            for snp in info:
                print snp

                for tissue in glob.glob("/users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/*.gz"):
                    # Look up this SNP in this tissue to see if it's any good
                    matches = subprocess.check_output("tabix {0} chr{1}:{2}-{2}".format(tissue, snp[0], snp[1]), shell=True)
                    if matches == "":
                        continue
                    matches = matches.strip().split("\n")
                    for match in matches:
                        data = match.strip().split("\t")
                        if float(data[5]) < 1e-6:
                            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(snp[0], snp[1], file, tissue, snp[3], snp[2], data[5], data[1]))
                            w.flush()


def snps_by_threshold(gwas_file, gwas_threshold, trait, window=1000000):

    snp_counts = {}
    hit_counts = {}

    with gzip.open(gwas_file) as f:
        header = f.readline().strip().split()

        trait_index = -1
        if "trait" in header:

            trait_index = header.index("trait")

        pval_index = header.index("pvalue")
        chr_index = header.index("chr")
        snp_pos_index = header.index("snp_pos")

        all_snps = []

        #i = 0
        for line in f:
            #i += 1
            #if i > 1000000:
            #    break
            data = line.strip().split("\t")
            if trait_index != -1:
                trait = data[trait_index]
            try:
                pvalue = float(data[pval_index])
            except:
                continue
            chr = data[chr_index]
            pos = int(data[snp_pos_index])
            snp_counts[trait] = snp_counts.get(trait, 0) + 1
            if pvalue > gwas_threshold:
                continue

            all_snps.append((chr, pos, pvalue, trait))
    
    # For now, include only autosomal SNPs.
    filtered = []
    for s in all_snps:
        if "chr" in str(s[0]):
            filtered.append((s[0][3:], s[1], s[2], s[3]))
        else:
            filtered.append((s[0], s[1], s[2], s[3]))

    all_snps = sorted(filtered, key=operator.itemgetter(2)) 

    # Go through the list of SNPs in order, adding the ones
    # passing our criteria.
    snps_to_test = []
    for snp in all_snps:

        # For now, ignore a SNP if it's in the MHC region -- this
        # would require alternative methods.
        if (snp[0] == "6") and snp[1] > 25000000 and snp[1] < 35000000:
                continue

        # Before adding a SNP, make sure it's not right next
        # to another SNP that we've already selected.
        skip = False
        for kept_snp in snps_to_test:
                if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < window and kept_snp[3] == snp[3]:
                        skip = True
                        break
        if not skip:
            snps_to_test.append(snp)
            
    return(snps_to_test)
        
if __name__ == "__main__":
    main()
