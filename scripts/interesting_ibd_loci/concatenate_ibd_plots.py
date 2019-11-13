import glob
import subprocess

for file in glob.glob("/users/mgloud/projects/brain_gwas/output/rna-editing-ibd-plots/*/plots/*"):
    subprocess.call("cp -r {0} /users/mgloud/projects/rna_editing/output/plots/ibd-plots".format(file), shell=True)
