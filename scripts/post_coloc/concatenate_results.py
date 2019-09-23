import glob
import subprocess

# Collect results from all analysis (they're all just done with COLOC at this point for comparability)

clpp_file = glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/*/*coloc*")[0]
with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "w") as w:
    with open(clpp_file) as f:
        w.write(f.readline())


with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "a") as a:
    for clpp_file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/*/*coloc*"):
        with open(clpp_file) as f:
            f.readline()
            for line in f:
                a.write(line)

for file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/*/plots/*"):
    subprocess.call("cp -r {0} /users/mgloud/projects/rna_editing/output/plots/aggregated_clpp".format(file), shell=True)