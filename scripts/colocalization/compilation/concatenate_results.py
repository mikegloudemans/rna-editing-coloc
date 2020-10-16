import glob
import subprocess

# Collect results from all analysis (they're all just done with COLOC at this point for comparability)

header_file = "output/colocalization/rna-editing-revisions/2020-10-04_14-09-43.014967_rna-coloc/rna-gwas-overlap_all-gwas_gtex_source-pval-min-0_source-pval-max-5e-08_lookup-pval1e-05_source-window500000_lookup-window10000_coloc-tests_phase8_coloc_status.txt"
coloc_output_files = glob.glob("output/colocalization/rna-editing-revisions/*/*coloc_status.txt")
coloc_concat_file = "output/colocalization/rna-editing-revisions/concatenated/all_coloc_results.txt"

# No line should be written twice
lines = set([])

with open(coloc_concat_file, "w") as w:
	with open(header_file) as f:
		header = f.readline()
		w.write(header)
		lines.add(header)

	for coloc_file in coloc_output_files:
		with open(coloc_file) as f:
			for line in f:
				# Write the line only if it hasn't been written before
				if line not in lines:
					lines.add(line)
					w.write(line)

# If needing to copy plots too, then fix this up
'''
for file in glob.glob("{0}/*/plots/*".format(coloc_raw_output_directory)):
    subprocess.call("cp -r {0} {1}/plots/".format(file, coloc_processed_output_directory), shell=True)
'''
