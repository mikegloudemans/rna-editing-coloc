import glob
import subprocess

# Collect results from all analysis (they're all just done with COLOC at this point for comparability)

# Raw input files
coloc_raw_output_directory = "../../../output/colocalization/main_coloc_results/raw_coloc_output"
# Where to place processed output files
coloc_processed_output_directory = "../../../output/colocalization/main_coloc_results/aggregated"

results = subprocess.check_output("ls {0}".format(coloc_raw_output_directory), shell=True).strip().split()

clpp_file = glob.glob("{1}/{0}/*coloc*".format(results[0], coloc_raw_output_directory))[0]
with open("{0}/aggregated_coloc_results.txt".format(coloc_processed_output_directory), "w") as w:
    with open(clpp_file) as f:
        w.write(f.readline())

with open("{0}/aggregated_coloc_results.txt".format(coloc_processed_output_directory), "a") as a:
    for res in results:
        for clpp_file in glob.glob("{1}/{0}/*coloc*".format(res, coloc_raw_output_directory)):
            with open(clpp_file) as f:
                f.readline()
                for line in f:
                    a.write(line)

for file in glob.glob("{0}/*/plots/*".format(coloc_raw_output_directory)):
    subprocess.call("cp -r {0} {1}/plots/".format(file, coloc_processed_output_directory), shell=True)

# The following section is only needed if the colocalization took multiple runs to complete
'''
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
'''





