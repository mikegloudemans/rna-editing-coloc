import gzip
import glob

snp_map = {}
snp_set = set([])

# Load all variants into memory with ref/alt definitions
for i in range(1, 23):
    print i
    with gzip.open("../../data/1KG/hg38/ALL.chr{0}_GRCh38.genotypes.20170504.vcf.gz".format(i)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            data = line.strip().split()
            if (data[0], data[1]) in snp_set:
                try:
                    del snp_map[(data[0], data[1])]
                except:
                    pass
            else:
                snp_map[(data[0], data[1])] = (data[3], data[4])
                snp_set.add((data[0], data[1]))

# Now cycle through eQTL files
for tissue_dir in glob.glob("../data/eqtls/*"):
    tissue = tissue_dir.split("/")[-1]

    with gzip.open("../data/eqtls/{0}/{0}.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.with.alleles.txt.gz".format(tissue), "w") as w:
        with gzip.open("../data/eqtls/{0}/{0}.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.txt.gz".format(tissue)) as f:
            for line in f:
                data = line.strip().split()
                # Only write lines if we can obtain ref/alt info for them
                if (data[8].replace("chr", ""), data[9]) in snp_map:
                    variants = snp_map[(data[8].replace("chr", ""), data[9])]
                    w.write(line.strip() + " {0} {1}\n".format(variants[0], variants[1]))

# Now cycle through eQTL aggro files
for tissue_dir in glob.glob("../data/edqtl_clusters/editingQLTs_agg_lancaster/*"):
    tissue = tissue_dir.split("/")[-1].split(".")[0]
    print tissue

    if tissue == "All":
        continue

    with gzip.open("../data/edqtl_clusters/editingQLTs_agg_lancaster/{0}.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Lancaster_agg.cluster.with.alleles.txt.gz".format(tissue), "w") as w:
        with gzip.open("../data/edqtl_clusters/editingQLTs_agg_lancaster/{0}.edMat.20cov.60samps.noXYM.qqnorm.bed.QTLtools.nominal.Lancaster_agg.cluster.txt.gz".format(tissue)) as f:
            w.write(f.readline().strip() + "\tref\talt\n")
            for line in f:
                data = line.strip().split()
                # Only write lines if we can obtain ref/alt info for them
                if (data[2].replace("chr", ""), data[3]) in snp_map:
                    variants = snp_map[(data[2].replace("chr", ""), data[3])]
                    w.write(line.strip() + "\t{0}\t{1}\n".format(variants[0], variants[1]))




