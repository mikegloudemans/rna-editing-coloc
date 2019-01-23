import glob
import subprocess

clpp_file = glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests/*/*clpp*")[0]
print clpp_file
with open("/users/mgloud/projects/rna_editing/output/single_site_clpp_results.txt", "w") as w:
    with open(clpp_file) as f:
        w.write(f.readline())


with open("/users/mgloud/projects/rna_editing/output/single_site_clpp_results.txt", "a") as a:
    for clpp_file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests/*/*clpp*"):
        with open(clpp_file) as f:
            f.readline()
            for line in f:
                a.write(line)

for file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-tests/*/plots/*"):
    subprocess.call("cp -r {0} /users/mgloud/projects/rna_editing/output/plots/single_site_clpp".format(file), shell=True)
