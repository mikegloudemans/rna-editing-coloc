import glob
import gzip

# Simple script to keep just the top SNP for each edQTL feature
edqtls = glob.glob("/users/mgloud/projects/rna_editing/data/edqtl_signif/*.gz")

for e in edqtls:
    top = {}    # Store top signifiicant edQTLs
    tissue = e.strip().split(".")[0].split("/")[-1]
    with gzip.open("/users/mgloud/projects/rna_editing/data/edqtl_signif/top_features/{0}.edFeats.txt.gz".format(tissue), "w") as w:
        with gzip.open(e) as f:
            w.write("chr\tvariant_pos\tref\talt\tbuild\t" + "\t".join(f.readline().strip().split("\t")[1:]) + "\n")
            for line in f:
                data = line.strip().split()
                feature = data[1]
                pval = float(data[-1])

                if feature in top:
                    #print feature, top[feature], pval
                    assert top[feature] >= pval
                else:
                    top[feature] = pval
                    line = "\t".join(data[0].split("_") + data[1:]) + "\n"
                    w.write(line)
