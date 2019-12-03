import glob
import subprocess

# Collect results from all analysis (they're all just done with COLOC at this point for comparability)

results = subprocess.check_output("ls /users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc", shell=True).strip().split()

clpp_file = glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/{0}/*coloc*".format(results[0]))[0]
with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "w") as w:
    with open(clpp_file) as f:
        w.write(f.readline())

with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "a") as a:
    for res in results:
        for clpp_file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/{0}/*coloc*".format(res)):
            with open(clpp_file) as f:
                f.readline()
                for line in f:
                    a.write(line)

other_results = subprocess.check_output("ls /users/mgloud/projects/brain_gwas/output/completed/rna_editing/moved-2019-11-13/rna-editing-tests-all-coloc", shell=True).strip().split()
with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "a") as a:
    for res in other_results:
        for clpp_file in glob.glob("/users/mgloud/projects/brain_gwas/output/completed/rna_editing/moved-2019-11-13/rna-editing-tests-all-coloc/{0}/*coloc*".format(res)):
            with open(clpp_file) as f:
                f.readline()
                for line in f:
                    a.write(line)

other_results = subprocess.check_output("ls /users/mgloud/projects/brain_gwas/output/rna-editing-tests-immune-bonus", shell=True).strip().split()
with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt", "a") as a:
    for res in other_results:
        for clpp_file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-immune-bonus/{0}/*coloc*".format(res)):
            with open(clpp_file) as f:
                f.readline()
                for line in f:
                    a.write(line)






#for file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests-all-coloc/*/plots/*"):
#    subprocess.call("cp -r {0} /users/mgloud/projects/rna_editing/output/plots/aggregated_clpp".format(file), shell=True)
