import subprocess
import gzip

# Write file header

# Run "overlap" script to pull top edQTLs from each tissue
#subprocess.check_call(["python", "/users/mgloud/projects/gwas-download/overlap/list_snps_to_test.py", "/users/mgloud/projects/rna_editing/scripts/edqtl_prs/top-edqtls.config"])

# Look up these edQTLs in the corresponding tissue using tabix
# Write results out to a file

gtex_genos = ""

with open("/users/mgloud/projects/rna_editing/output/top-edqtls/prs_list.txt", "w") as w:

    # Write header
    with gzip.open("/users/mgloud/projects/rna_editing/data/tabix_eqtls/Muscle-Skeletal.edQTLs.txt.gz") as head:
        w.write("file\t" + head.readline())
    with open("/users/mgloud/projects/rna_editing/output/top-edqtls/top-edqtls_all-gwas_gwas-pval1e-20_gwas-window500000_snps-considered.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            
            # Tabix lookup and write to file for later reference
            tbx_result = subprocess.check_output("tabix {0} chr{1}:{2}-{2} | grep {3}".format(data[4], data[0], data[1], data[3]), shell=True)
            w.write("{0}\t".format(data[4]) + tbx_result)
            tbx_split = tbx_result.split("\t")

            # TODO: Somehow make sure we're actually looking at the right variant and not
            # just some other one at the same site

            # Look up all these files in GTEx, concatenate them into a single table
            vcf_line = subprocess.check_output("tabix /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz chr{0}:{1}-{1}".format(data[0], data[1]), shell=True)

            vcf_line = [g.split(":")[0] for g in vcf_line.split("\t")]
            # Genotypes only
            genos = vcf_line[9:]
            num_alts = [sum(["1" == s for s in g.split("/")]) for g in genos]
            num_refs = [sum(["0" == s for s in g.split("/")]) for g in genos]

            # Figure out which ones have the risk allele
            # For missing individuals, just assign them the average I guess...

            # (For now, don't worry about the fact that many individuals have contributed more than one tissue)


# View the overall distribution of # of top SNPs in different individuals
