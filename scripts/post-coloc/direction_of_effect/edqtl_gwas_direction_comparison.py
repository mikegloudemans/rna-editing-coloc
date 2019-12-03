# Author: Mike Gloudemans
# Date created: 2019-09-04

import subprocess
import operator
import gzip
import sys
import glob
from scipy import stats
import numpy as np
from multiprocessing import Pool
import traceback

# Declare constants

top_only_mode = True

max_cores = 24

gwas_pval_threshold = 1e-5
#gwas_pval_threshold = 5e-8

max_edqtl_distance = 100000

# List GWAS to test (or even just do all?)
all_gwas = glob.glob("/users/mgloud/projects/rna_editing/data/gwas/*.txt.gz")

all_gwas = [ag for ag in all_gwas if "Okada_2014" in ag or "Multiple-Sclerosis" in ag or "Ulcerative-Colitis" in ag or "Coronary-Artery-Disease_Nelson" in ag or "Crohns" in ag or "Inflammatory-Bowel-Disease" in ag or "Type-2-Diabetes" in ag]

# List directory with eQTLs
all_eqtls = glob.glob("/mnt/lab_data/montgomery/mgloud/rna_editing/data/tabix_eqtls/*.txt.gz")

# temporary

def main():

    if not top_only_mode:
        valid_qtls = get_eqtl_set(all_eqtls)
    else:
        valid_qtls = set([])

    coloc_results = get_coloc_snps()
    print set([cr[3] for cr in coloc_results])

    # For each GWAS file...
    all_test_agreements = {}
    for gwas in all_gwas:
        print
        print 
        print gwas
        # Select GWAS SNPs to compare
        # 
        # Two options for this: could go through the genome in a tiled fashion, breaking into 1MB chunks,
        # pulling best from each. OR could get top SNPs for the trait up to some certain level,
        # and see what this gives us
        #
        # Let's try the second option to start, since we have code to do it already
        #

        all_traits = set([])
        # Get list traits for GWAS file
        with gzip.open(gwas) as f:
            head = f.readline().strip().split("\t")
            if "effect_direction" not in head or "effect_allele" not in head or "non_effect_allele" not in head:
                continue
            if "trait" in head:
                trait_index = head.index("trait")
                i = 0
                for line in f:
                    if i > 100000:
                        break
                    data = line.strip().split("\t")
                    all_traits.add(data[trait_index])
                    i += 1
            else:
                all_traits.add(gwas)

        all_traits = list(all_traits)

        # For all traits measured in this single GWAS file...
        for trait in all_traits:
            print 
            print trait
            gwas_hits = snps_by_threshold(valid_qtls, gwas, gwas_pval_threshold, trait)
            # Test edQTL directions!
            test_results = test_edqtl_directions(gwas_hits)

            for tr in test_results:
                if tr not in all_test_agreements:
                    all_test_agreements[tr] = {"+": 0, "na": 0, "-": 0}
                for item in test_results[tr]:
                    all_test_agreements[tr][item] += test_results[tr][item]

            # Give only the colocs for the current trait
            print "Coloc results only:"
            coloc_sub = set([(cr[0], cr[1]) for cr in list(coloc_results) if cr[3] == trait])
            test_results = test_edqtl_directions(gwas_hits, coloc_sub)

            for tr in test_results:
                if tr + "coloc_only" not in all_test_agreements:
                    all_test_agreements[tr + "coloc_only"] = {"+": 0, "na": 0, "-": 0}
                for item in test_results[tr]:
                    all_test_agreements[tr + "coloc_only"][item] += test_results[tr][item]

            print "Tiled mode:",

            gwas_hits = snps_by_tiling(valid_qtls, gwas, trait)
            # Test edQTL directions!
            test_results = test_edqtl_directions(gwas_hits)

            for tr in test_results:
                if tr + "tiled" not in all_test_agreements:
                    all_test_agreements[tr + "tiled"] = {"+": 0, "na": 0, "-": 0}
                for item in test_results[tr]:
                    all_test_agreements[tr + "tiled"][item] += test_results[tr][item]

    # See if binomial test passes in either direction (two-tailed)
    print "All tests together:"
    for mode in all_test_agreements:
        print mode
        binom_k = min(all_test_agreements[mode]["+"], all_test_agreements[mode]["-"])
        # Get prob of seeing a result more extreme than this, in either direction
        binom_p = 2 * stats.binom.cdf(binom_k, all_test_agreements[mode]["+"] + all_test_agreements[mode]["-"], 0.5)   # Fails in case where they're exactly even, but that's OK because that's insignificant anyway
            
        print all_test_agreements[mode], binom_p

def snps_by_threshold(valid_qtls, gwas_file, gwas_threshold, default_trait, window=1000000):

    trait = gwas_file

    with gzip.open(gwas_file) as f:
        header = f.readline().strip().split()

        trait_index = -1
        if "trait" in header:
            trait_index = header.index("trait")
        pval_index = header.index("pvalue")
        chr_index = header.index("chr")
        snp_pos_index = header.index("snp_pos")
        direction_index = header.index("effect_direction")
        alt_index = header.index("effect_allele")
        ref_index = header.index("non_effect_allele")

        all_snps = []

        for line in f:
            data = line.strip().split("\t")
            if trait_index != -1:
                trait = data[trait_index]
                if trait != default_trait:
                    continue
            try:
                pvalue = float(data[pval_index])
            except:
                continue

            chr = data[chr_index]
            pos = data[snp_pos_index]
            
            if not top_only_mode:
                if (chr.replace("chr", ""), pos) not in valid_qtls:
                    continue

            pos = int(pos)

            if pvalue > gwas_threshold:
                continue



            direction = data[direction_index]
            ref = data[ref_index]
            alt = data[alt_index]

            all_snps.append((chr, pos, pvalue, trait, direction, ref, alt))
    
    # For now, include only autosomal SNPs.
    filtered = []
    for s in all_snps:
        if "chr" in str(s[0]):
            filtered.append((s[0][3:], s[1], s[2], s[3], s[4], s[5], s[6]))
        else:
            filtered.append((s[0], s[1], s[2], s[3], s[4], s[5], s[6]))

    all_snps = sorted(filtered, key=operator.itemgetter(2)) 

    # Go through the list of SNPs in order, adding the ones
    # passing our criteria.
    snps_to_test = []
    for snp in all_snps:

        # For now, ignore a SNP if it's in the MHC region -- this
        # would require alternative methods.
        if (snp[0] == "6") and snp[1] > 25000000 and snp[1] < 35000000:
                continue

        # Before adding a SNP, make sure it's not right next
        # to another SNP that we've already selected.
        
        skip = False
        for kept_snp in snps_to_test:
                if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < window and kept_snp[3] == snp[3]:
                        skip = True
                        break
        if not skip:
            snps_to_test.append(snp)
    
    return(snps_to_test)


def snps_by_tiling(valid_qtls, gwas_file, default_trait, window=1000000):

    trait = gwas_file

    with gzip.open(gwas_file) as f:
        header = f.readline().strip().split()

        trait_index = -1
        if "trait" in header:
            trait_index = header.index("trait")
        pval_index = header.index("pvalue")
        chr_index = header.index("chr")
        snp_pos_index = header.index("snp_pos")
        direction_index = header.index("effect_direction")
        alt_index = header.index("effect_allele")
        ref_index = header.index("non_effect_allele")
        
        snps_to_test = []

        current_chrom = "1"
        current_position = 0
        best_pval = 1
        best_snp = None
        for line in f:
            data = line.strip().split("\t")
            if trait_index != -1:
                trait = data[trait_index]
                if trait != default_trait:
                    continue
            try:
                pvalue = float(data[pval_index])
            except:
                continue
            chr = data[chr_index].replace("chr", "")
            pos = data[snp_pos_index]
            if (chr, pos) not in valid_qtls:
                continue

            pos = int(pos)

            # If we're on the next chromosome, reset
            if current_chrom != chr:
                if not (best_snp is None):
                    snps_to_test.append(best_snp)
                current_chrom = chr
                current_position = 0
                best_pval = 1
                best_snp = None

            if pos - current_position > window:
                if not (best_snp is None):
                    snps_to_test.append(best_snp)
                current_position = pos - (pos % window)
                best_pval = 1
                best_snp = None

            if pvalue < best_pval:
                direction = data[direction_index]
                ref = data[ref_index]
                alt = data[alt_index]
                best_pval = pvalue
                best_snp = (chr, pos, pvalue, trait, direction, ref, alt)

    return(snps_to_test)

def test_edqtl_directions(gwas_hits, coloc_results = None):

    agreement = {}

    modes = ["binomial", \
            "best", \
            "majority",
            "combined-tissue", \
            "combined-tissue-and-sites"]
    for mode in modes:
        agreement[mode] = {"+": 0, "na": 0, "-": 0}

    # Run key SNPs in parallel
    pool = Pool(max_cores)
    results = pool.map(test_hit_wrapper, [(gh, coloc_results) for gh in gwas_hits])

    all_beta = []
    for result_pair in results:
        result = result_pair[0]
        for mode in result:
            for value in result[mode]:
                agreement[mode][value] += result[mode][value]
        if not (result_pair[1] is None):
            all_beta.append(result_pair[1])

    print "1-sample T-test with weights:"
    print stats.ttest_1samp(all_beta, 0)

    ranks = stats.rankdata([abs(a) for a in all_beta])
    neg_ranks = [ranks[i] for i in range(len(ranks)) if all_beta[i] < 0]
    pos_ranks = [ranks[i] for i in range(len(ranks)) if all_beta[i] > 0]
    print "Wilcoxon rank-sum test, neg. vs. positive ranks"
    print stats.ranksums(neg_ranks, pos_ranks)

    # See if binomial test passes in either direction (two-tailed)
    for mode in agreement:
        print mode
        binom_k = min(agreement[mode]["+"], agreement[mode]["-"])
        # Get prob of seeing a result more extreme than this, in either direction
        binom_p = 2 * stats.binom.cdf(binom_k, agreement[mode]["+"] + agreement[mode]["-"], 0.5)   # Fails in case where they're exactly even, but that's OK because that's insignificant anyway
        print agreement[mode], binom_p

    

    return agreement


def test_hit_wrapper(params):
    hit = params[0]
    coloc_results = params[1]
    try:
        return test_hit(hit, coloc_results)
    except:
        traceback.print_exc(file=sys.stdout)


def test_hit(hit, coloc_results):

    #print hit
    # Get the top editing site across all tissues
    # But also save all editing sites in all tissues
    all_edits = {}
    my_agreement = {}
    modes = ["binomial", \
            "best", \
            "majority",
            "combined-tissue", \
            "combined-tissue-and-sites"]
    for mode in modes:
        my_agreement[mode] = {"+": 0, "na": 0, "-": 0}

    if not (coloc_results is None):
        print list(coloc_results)[:3]
        print (str(hit[0]), str(hit[1]).replace("chr", ""))
        if not (str(hit[0]), str(hit[1]).replace("chr", "")) in coloc_results:
            for mode in my_agreement:
                my_agreement[mode]["na"] += 1
            return my_agreement, None


    best_pval = 1
    best_feature = None
    best_dir = "na"
    for aq in all_eqtls:
        out = subprocess.check_output(["tabix", aq, "chr{0}:{1}-{2}".format(hit[0], hit[1], hit[1])]).strip()
        if out == "":
            continue
        info = out.strip().split("\n")
        for line in info:

            data = line.strip().split()

            # Make sure it's close enough to be worth testing
            if abs(int(data[6])) > max_edqtl_distance:
                continue

            edit_site = data[0]
            pval = float(data[11])
            beta = float(data[12])

            zscore = stats.norm.ppf(pval / 2)
            if beta > 0:
                zscore *= -1

            if edit_site not in all_edits:
                all_edits[edit_site] = []
            all_edits[edit_site].append(round(zscore,2))

            if pval < best_pval:
                best_pval = pval
                best_feature = edit_site
                best_dir = beta

    if best_pval == 1:
        #print "No edQTL tested"
        for mode in my_agreement:
            my_agreement[mode]["na"] += 1
        return my_agreement, None
    
    if best_dir > 0:
        my_agreement["best"]["+"] += 1
    elif best_dir < 0:
        my_agreement["best"]["-"] += 1
    else:
        my_agreement["best"]["na"] += 1

    # Stouffer's combined
    # Start by just taking the MEAN of all Z-scores, since
    # we don't really

    # Combining the best site (from one individual tissue) across all tissues
    tissue_combined = sum(all_edits[best_feature]) / (len(all_edits[best_feature])**(0.5))

    #print all_edits[best_feature]
    #print tissue_combined

    if tissue_combined > 0:
        my_agreement["combined-tissue"]["+"] += 1
    elif tissue_combined < 0:
        my_agreement["combined-tissue"]["-"] += 1
    else:
        my_agreement["combined-tissue"]["na"] += 1

    # Combined each site across all tissues, and then combining all sites
    z_vec = []
    for site in all_edits:
        tissue_combined = sum(all_edits[site]) / (len(all_edits[site]) ** (0.5))
        z_vec.append(round(tissue_combined,2))
        #print tissue_combined, all_edits[site]
    tissue_combined_site_combined = sum(z_vec) / (len(z_vec) ** (0.5))
    #print "-------"
    #print z_vec
    #print tissue_combined_site_combined
    #print

    if tissue_combined_site_combined > 0:
        my_agreement["combined-tissue-and-sites"]["+"] += 1
    elif tissue_combined_site_combined < 0:
        my_agreement["combined-tissue-and-sites"]["-"] += 1
    else:
        my_agreement["combined-tissue-and-sites"]["na"] += 1


    #print best_pval

    # Then at that top editing site, count support in favor of either
    # direction, across tissues.
    # For now, we'll only look at that one single editing site.
    this_snp = {"+": 0, "-": 0}
    for aq in all_eqtls:
        out = subprocess.check_output(["tabix", aq, "chr{0}:{1}-{2}".format(hit[0], hit[1], hit[1])]).strip()
        if out == "":
            # This feature wasn't tested for this tissue at this SNP
            continue
        info = out.strip().split("\n")
        for line in info:
            data = line.strip().split()
            if data[0] != best_feature:
                continue
            # Once we find the effect size of the GWAS SNP on the
            # feature of interest in this tissue...
            beta = float(data[12])
            #print beta,
            ref = data[14].upper()
            alt = data[15].upper()

            if hit[4] == "+":
                gwas_alt_dir = True
            elif hit[4] == "-":
                gwas_alt_dir = False
            else:
                #print "Invalid GWAS effect direction", hit[4]
                continue

            qtl_alt_dir = beta > 0

            if hit[5].upper() == ref and hit[6].upper() == alt:
                pass
            elif hit[5].upper() == alt and hit[6].upper() == ref:
                # Flip direction to match pval with GWAS
                #print "flip", ref, alt,
                qtl_alt_dir = not qtl_alt_dir
            else:
                #print "Ref/alt don't match:", hit[5], hit[6], alt, ref
                continue
            same_direction = not (qtl_alt_dir ^ gwas_alt_dir)
            if same_direction:
                this_snp["+"] += 1
            else:
                this_snp["-"] += 1
    
    # Binomial test mode
    # Take this test if binomial test passes in either direction (two-tailed)
    binom_k = min(this_snp["+"], this_snp["-"])
    # Get prob of seeing a result more extreme than this, in either direction
    binom_p = 2 * stats.binom.cdf(binom_k, this_snp["+"]+this_snp["-"], 0.5)   # Fails in case where they're exactly even, but that's OK because that's insignificant anyway
    
    if this_snp["+"] > this_snp["-"] and binom_p < 0.1:
        my_agreement["binomial"]["+"] += 1
    elif this_snp["-"] > this_snp["+"] and binom_p < 0.1:
        my_agreement["binomial"]["-"] += 1
    else:
        my_agreement["binomial"]["na"] += 1

    # No binomial mode
    if this_snp["+"] > this_snp["-"]:
       my_agreement["majority"]["+"] += 1
    elif this_snp["-"] > this_snp["+"]:
        my_agreement["majority"]["-"] += 1
    else:
        my_agreement["majority"]["na"] += 1

    return my_agreement, tissue_combined 


def get_eqtl_set(qtls):
    snps = set([])
    for q in qtls:
        print q
        with gzip.open(q) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                # We're going to limit how far away the SNPs can actually be...
                if abs(int(data[6])) > max_edqtl_distance:
                    continue
                snps.add((data[8].replace("chr", ""), data[9]))

    return snps

def get_coloc_snps():
    coloc_snps = set([])
    with open("/users/mgloud/projects/rna_editing/output/aggregated_coloc_results_all.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split()

            # Has to be an editing QTL locus
            if "eQTL" in data[1] or "sQTL" in data[1]:
                continue

            # Get only colocalizations
            if float(data[9]) > 0.9:
                chrom, snp = data[0].split("_")
                coloc_snps.add((chrom, snp, data[1], data[2], data[3]))
    return coloc_snps

if __name__ == "__main__":
    
    main()
