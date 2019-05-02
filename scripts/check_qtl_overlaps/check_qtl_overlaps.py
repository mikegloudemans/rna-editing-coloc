import glob
import gzip
import subprocess
import os
import sys
import time
import traceback

from multiprocessing import Pool

# Test QTLs from a source file to see whether they're also
# QTLs of a different (target) type, in the same tissue

output_file = "/users/mgloud/projects/rna_editing/output/qtl_comparisons.txt"
max_cores = 15

## Dictionary mapping source file to target file (target is the tabixed file
## where we want to try to replicate our QTLs)
significant_sites = {
                        "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/*egenes.txt.gz":
                        [   
                            "/users/mgloud/projects/rna_editing/data/tabix_eqtls/{0}.edQTLs.txt.gz",
                            "/users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/{0}.Fisher.combined.sorted.txt.gz",
                            "/users/mgloud/projects/brain_gwas/data/sqtls/gtex_v8/{0}.allpairs.txt.gz.sQTLs.txt.gz"
                        ],
                        "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sqtl/GTEx_Analysis_v8_sQTL/*sgenes.txt.gz":
                        [
                            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v8/{0}.allpairs.txt.gz.eQTLs.txt.gz",
                            "/users/mgloud/projects/rna_editing/data/tabix_eqtls/{0}.edQTLs.txt.gz",
                            "/users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/{0}.Fisher.combined.sorted.txt.gz"
                        ],
                        "/users/mgloud/projects/rna_editing/data/edqtl_signif/top_features/{0}.edFeats.txt.gz":
                        [
                            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v8/{0}.allpairs.txt.gz.eQTLs.txt.gz",
                            "/users/mgloud/projects/rna_editing/data/tabix_eqtls_aggro/{0}.Fisher.combined.sorted.txt.gz",
                            "/users/mgloud/projects/brain_gwas/data/sqtls/gtex_v8/{0}.allpairs.txt.gz.sQTLs.txt.gz"
                        ]
                    }


def main():
    # Write header of output file
    with open(output_file, "w") as w:
        w.write("source_file\ttarget_file\tsource_id\tsource_feature\ttarget_feature\ttarget_pval\n")

    # For each source QTL study
    for source in significant_sites:

        all_source = glob.glob(source)

        # For each file in source QTL study...
        for source_qtl in all_source:
            print source_qtl

            tissue = source_qtl.split("/")[-1].split(".")[0]

            # Get all target files corresponding to the same tissue
            all_target = []
            for s in significant_sites[source]:
                # Check what separation characters are used for tissue spaces
                base_dir = "/".join(s.strip().split("/")[:-1])
                examples = glob.glob(base_dir + "/*.gz")
                if sum(["-" in e for e in examples]) > 0:
                    all_target.append(s.format(tissue.replace("_", "-")))
                else:
                    all_target.append(s.format(tissue))

            # Go through file and store significant SNPs for each study
            with gzip.open(source_qtl) as f:
                header = f.readline().strip().split()

                chrom_column = header.index("chr")
                pos_column = header.index("variant_pos")
                id_column = header.index("gene_id")
                if "gene_name" in header:
                    feature_column = header.index("gene_name")
                else:
                    feature_column = id_column

                # Do this part in parallel
                pool = Pool(max_cores)
                for line in f:
                    pool.apply_async(parallel_wrapper, args=(line, source_qtl, chrom_column, pos_column, id_column, feature_column, all_target))
                pool.close()
                pool.join()

def parallel_wrapper(line, source_qtl, chrom_column, pos_column, id_column, feature_column, all_target):
    try:
        process_lines_parallel(line, source_qtl, chrom_column, pos_column, id_column, feature_column, all_target)
    except Exception:
        traceback.print_exc(file=sys.stdout)


def process_lines_parallel(line, source_qtl, chrom_column, pos_column, id_column, feature_column, all_target):
    chunks = line.strip().split()
    print source_qtl
    print chunks[:5]

    chrom = chunks[chrom_column].replace("chr","")
    pos = chunks[pos_column]
    id = chunks[id_column]
    feature = chunks[feature_column]

    for t in range(len(all_target))[::-1]:
        if not os.path.isfile(all_target[t]):
            del all_target[t]

    # For each target file, get the corresponding test(s)
    for target_qtl in all_target:

        if not os.path.isfile(target_qtl):
            continue

        # TODO: Move this...
        with gzip.open(target_qtl) as tq:
            header = tq.readline().strip().split()
            second = tq.readline()

        header = [h.lower() for h in header]
        if "gene" in header:
            target_feature_index = header.index("gene")
        else:
            target_feature_index = header.index("feature")
        target_pval_index = header.index("pvalue")

        if "chr" in second:
            data = subprocess.check_output("tabix {0} chr{1}:{2}-{2}".format(target_qtl, chrom, pos), shell=True)

        else:
            data = subprocess.check_output("tabix {0} {1}:{2}-{2}".format(target_qtl, chrom, pos), shell=True)

        if data == "":
            continue

        # Now write the relevant pvalue and metadata to the output file
        for d in data.strip().split("\n"):
            s = d.split("\t")
            target_feature = s[target_feature_index]
            target_pval = s[target_pval_index]
            with open(output_file, "a") as a:
                a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(source_qtl, target_qtl, id, feature, target_feature, target_pval))


if __name__ == "__main__":
    main()
