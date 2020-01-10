# Author: Mike Gloudemans
# Date created: 2019-09-04

import os
import subprocess
import operator
import gzip
import sys
import glob
from scipy import stats
import numpy as np
from multiprocessing import Pool
import traceback

import pickle

#############################
# Declare settings
#############################

# Output z-scores for every site?
output_full_zscores = True

# Should we look only at GWS hits, rather than tile across
# the entire genome and test effects in every 1MB window?
skip_tiled_mode = False

# Should we skip the step of filtering down to only colocalized
# SNPs?
skip_coloc_filtering = False

# If the lead GWAS SNP wasn't even tested for edQTLs, should
# we call it an NA? Or should be continue selecting nearby
# SNPs until one is found that was tested for edQTLs?
top_only_mode = False

max_cores = 36

gwas_pval_threshold = 1e-5
#gwas_pval_threshold = 5e-8

max_edqtl_distance = 100000

# How big of chunks should we break the genome into for the tiling approach?
tiling_window = 5000000

# Minimum H4 score at which we consider a locus
# to have colocalization
#min_coloc_threshold = 0.9
min_coloc_threshold = 0.5

#############################

# List GWAS to test (or even just do all?)
all_gwas = glob.glob("/users/mgloud/projects/rna_editing/data/gwas/*.txt.gz")

all_gwas = [ag for ag in all_gwas if "Okada_2014" in ag or "Multiple-Sclerosis" in ag or "Ulcerative-Colitis" in ag or "Coronary-Artery-Disease_Nelson" in ag or "Crohns" in ag or "Inflammatory-Bowel-Disease" in ag or "Type-2-Diabetes" in ag]
#all_gwas = [ag for ag in all_gwas if "Bowel" in ag]

# Add the immune ones that we downloaded just for this paper
all_gwas += glob.glob("/users/mgloud/projects/rna_editing/scripts/preprocessing/upgrade_gwas/munged/*/*.gz")

# Vitilogo GWAS has problems because not filtered properly for allele frequency; deal with this one later
all_gwas = [ag for ag in all_gwas if "Viti" not in ag]

# TODO: Add some negative controls (non-immune traits)

# List directory with eQTLs
all_eqtls = glob.glob("/mnt/lab_data/montgomery/mgloud/rna_editing/data/tabix_eqtls/*.txt.gz")

# temporary
#all_eqtls = all_eqtls[:3]

def main():

    if not top_only_mode:
        valid_qtls = get_eqtl_set(all_eqtls)
    else:
        valid_qtls = set([])


    if output_full_zscores:
        with open("zscore_table.txt", "w") as w:
            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format("gwas_chr", "gwas_pos", "gwas_pval", "gwas_trait", "edit_tissue", "edit_site", "gwas_risk_edqtl_zscore"))

    coloc_results = get_coloc_snps()
    # ^ list of snps in tuple form (chrom, snp, eqtl_file, gwas_trait, feature)
    # It's already filtered down to the list of SNPs that actually colocalized

    # Dictionary to record concordance agreement status for several different types of concordance tests
    # across ALL GWAS (not just a single study for this one)
    all_test_agreements = {}
    # For each GWAS file...
    for gwas in all_gwas:
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

        # Get list traits for GWAS file
        all_traits = get_traits_in_file(gwas)

        if all_traits is None:
            continue
        
        # For all traits measured in this single GWAS file...
        for trait in all_traits:

            ################################################
            # Test concordance, GWAS lead SNPs only,
            # no coloc filtering
            ################################################
            print 
            print gwas
            print trait
            # Get all SNPs with p-value below our chosen threshold
            gwas_hits = snps_by_threshold(valid_qtls, gwas, gwas_pval_threshold, trait)
            # Test edQTL directions!
            # Also, output all the z-scores for these tests if it's been requested
            test_results = test_edqtl_directions(gwas_hits, output_full_zscores = output_full_zscores)

            for tr in test_results:
                if tr not in all_test_agreements:
                    all_test_agreements[tr] = {"+": 0, "na": 0, "-": 0}
                for item in test_results[tr]:
                    all_test_agreements[tr][item] += test_results[tr][item]

            if not skip_coloc_filtering:

                ########################################################
                # Test concordance, GWAS lead SNPs only,
                # subsetting to SNPs colocalized for current GWAS trait
                ########################################################

                print "Coloc results only:"
                coloc_sub = set([(cr[0], cr[1]) for cr in list(coloc_results) if cr[3] == trait])
                test_results = test_edqtl_directions(gwas_hits, coloc_sub)

                for tr in test_results:
                    if tr + "coloc_only" not in all_test_agreements:
                        all_test_agreements[tr + "coloc_only"] = {"+": 0, "na": 0, "-": 0}
                    for item in test_results[tr]:
                        all_test_agreements[tr + "coloc_only"][item] += test_results[tr][item]

            # We can opt to skip the next part to save time
            if not skip_tiled_mode:
            
                ################################################
                # Test concordance, one SNP from every chunk,
                # no coloc filtering
                ################################################
                print "Tiled mode:",

                gwas_hits = snps_by_tiling(valid_qtls, gwas, trait, tiling_window)
                # Test edQTL directions!
                test_results = test_edqtl_directions(gwas_hits)

                for tr in test_results:
                    if tr + "tiled" not in all_test_agreements:
                        all_test_agreements[tr + "tiled"] = {"+": 0, "na": 0, "-": 0}
                    for item in test_results[tr]:
                        all_test_agreements[tr + "tiled"][item] += test_results[tr][item]

    # For ALL traits combined...
    # See if binomial test passes in either direction (two-tailed)
    print "All tests together:"
    for mode in all_test_agreements:
        print mode
        binom_k = min(all_test_agreements[mode]["+"], all_test_agreements[mode]["-"])
        # Get prob of seeing a result more extreme than this, in either direction
        binom_p = 2 * stats.binom.cdf(binom_k, all_test_agreements[mode]["+"] + all_test_agreements[mode]["-"], 0.5)   
        # Fails in case where they're exactly even, but that's OK because that's insignificant anyway
            
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
           
            # If this speed-up parameter isn't set, then we should keep
            # successively picking additional GWAS SNPs until we run out
            # or until we find one that was actually tested for eQTLs in
            # the first place
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
            
            # If this speed-up parameter isn't set, then we should keep
            # successively picking additional GWAS SNPs until we run out
            # or until we find one that was actually tested for eQTLs in
            # the first place
            if not top_only_mode:
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

def test_edqtl_directions(gwas_hits, coloc_results = None, output_full_zscores=False):

    # Names for all the ways we want to test the directional agreement
    modes = ["binomial", \
            "best", \
            "majority",
            "combined-tissue", \
            "combined-tissue-and-sites"]

    # Dictionary to show the number of agreements vs. disagreements for each
    # test type
    #
    # "agreement" or "+" means the allele that increases RNA editing levels
    # also increases the risk of disease.
    #
    # "disagreement" or "-" means the allele that increases RNA editing levels
    # decreases the risk of disease.
    #
    agreement = {}
    for mode in modes:
        agreement[mode] = {"+": 0, "na": 0, "-": 0}

    # Test SNPs for agreement in parallel, to save time
    pool = Pool(max_cores)
    results = pool.map(test_hit_wrapper, [(gh, coloc_results, output_full_zscores) for gh in gwas_hits])

    # Total up the test results for all different tested SNPs
    # and also keep a record of all the combined z-scores for the
    # tissue-combined tests
    all_beta = []
    for result_pair in results:
        result = result_pair[0]
        for mode in result:
            for value in result[mode]:
                agreement[mode][value] += result[mode][value]
        if not (result_pair[1] is None):
            all_beta.append(result_pair[1])

    # Is the average edQTL effect size of a GWAS SNP non-zero?
    print "1-sample T-test with weights:"
    print stats.ttest_1samp(all_beta, 0)

    # Not sure if this is that great of a test for this but...
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

# Nothing special, just a wrapper that prints the
# traceback message if our spawned threads fail.
def test_hit_wrapper(params):
    hit = params[0]
    coloc_results = params[1]
    output_full_zscores = params[2]
    try:
        return test_hit(hit, coloc_results, output_full_zscores)
    except:
        traceback.print_exc(file=sys.stdout)

# If coloc_results is None, include the hit no matter what 
# Otherwise, skip the hit if it is not colocalized
def test_hit(hit, coloc_results=None, output_full_zscores=False):

    # Z-scores for all editing sites in all tissues tested for association with the GWAS hit
    all_edits = {}
    # Results of concordance testing under various approaches
    my_agreement = {}
    # The different ways of testing for concordance
    modes = ["binomial", \
            "best", \
            "majority",
            "combined-tissue", \
            "combined-tissue-and-sites"]
    for mode in modes:
        my_agreement[mode] = {"+": 0, "na": 0, "-": 0}

    # If we are filtering by coloc results (i.e. we have passed a list
    # of coloc SNPs to filter down to), then we run this block to skip SNPs that
    # didn't colocalize. Otherwise, coloc_results defaults to None and we skip
    # this block.
    if not (coloc_results is None):
        if not (str(hit[0]), str(hit[1]).replace("chr", "")) in coloc_results:
            for mode in my_agreement:
                my_agreement[mode]["na"] += 1
            return my_agreement, None
        
    # In the GWAS, is the alt allele the one that increases risk?
    if hit[4] == "+":
        gwas_alt_increases_risk = True
    elif hit[4] == "-":
        gwas_alt_increases_risk = False
    else:
        #print "Invalid GWAS effect direction", hit[4]
        for mode in my_agreement:
            my_agreement[mode]["na"] += 1
        return my_agreement, None

    gwas_ref = hit[5].upper()
    gwas_alt = hit[6].upper()

    best_edqtl_pval = 1
    best_edit_site = None
    best_edqtl_direction = "na"

    # Loop through all edQTL tissues
    for edqtl_tissue in all_eqtls:
        tissue_short_name = edqtl_tissue.strip().split("/")[-1].split("\.")[0]
        # Get the editing sites that were tested for association with the GWAS SNP
        out = subprocess.check_output(["tabix", edqtl_tissue, "chr{0}:{1}-{2}".format(hit[0], hit[1], hit[1])]).strip()
        if out == "":
            continue
        info = out.strip().split("\n")
        for line in info:

            data = line.strip().split()

            # Make sure it's close enough to be worth testing
            if abs(int(data[6])) > max_edqtl_distance:
                continue

            edit_site = data[0]
            edqtl_pval = float(data[11])
            edqtl_beta = float(data[12])

            edqtl_ref = data[14].upper()
            edqtl_alt = data[15].upper()

            # In the QTL study, is the alt allele the one that increases editing?
            edqtl_alt_increases_editing = edqtl_beta > 0

            # If the GWAS and the QTL study have opposite designations for ref/alt,
            # then flip the direction of effect so it corresponds to matching 
            # ref/alt designation
            if gwas_ref == edqtl_ref and gwas_alt == edqtl_alt:
                gwas_alt_increases_editing = edqtl_alt_increases_editing
            elif gwas_ref == edqtl_alt and gwas_alt == edqtl_ref:
                # Flip direction to match pval with GWAS
                #print "flip", ref, alt,
                gwas_alt_increases_editing = not edqtl_alt_increases_editing
            else:
                # Ref/alt don't match
                continue

            # Does the editing-increasing allele also increase GWAS risk?
            risk_increases_editing = (gwas_alt_increases_risk == gwas_alt_increases_editing)

            # Convert the p-value to a signed z-score, which will be necessary
            # when we combine across tissues or across editing sites
            gwas_risk_edqtl_zscore = stats.norm.ppf(edqtl_pval / 2) 
            if not risk_increases_editing:
                # Get the direction of the effect
                gwas_risk_edqtl_zscore *= -1

            # Keep track of all the signed zscores for this editing site,
            # across tissues
            if edit_site not in all_edits:
                all_edits[edit_site] = []
            all_edits[edit_site].append(round(gwas_risk_edqtl_zscore,2))

            if output_full_zscores:
                with open("zscore_table.txt", "a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(hit[0], hit[1], hit[2], hit[3].split("/")[-1], tissue_short_name, edit_site, gwas_risk_edqtl_zscore))

            # Track which editing site has the strongest association with this GWAS hit,
            # across all sites and tissues
            if edqtl_pval < best_edqtl_pval:
                best_edqtl_pval = edqtl_pval
                best_edit_site = edit_site
                best_edqtl_direction = gwas_risk_edqtl_zscore

    if best_edqtl_pval == 1:
        # No edQTLs were tested for this site, so we can't say anything about it
        # for any of the possible concordance testing methods
        for mode in my_agreement:
            my_agreement[mode]["na"] += 1
        return my_agreement, None

    ###################################
    # Test: Best snp only
    ###################################
    
    # Test: If we look at only the MOST impactful edQTL association, was it an agreement
    # or a disagreement?
    if best_edqtl_direction > 0:
        my_agreement["best"]["+"] += 1
    elif best_edqtl_direction < 0:
        my_agreement["best"]["-"] += 1
    else:
        my_agreement["best"]["na"] += 1

    ###################################
    # Test: Best site, across tissues
    ##################################
    
    # Stouffer's combined z-score

    # Combining the best site (from one individual tissue) across all tissues
    tissue_combined = sum(all_edits[best_edit_site]) / (len(all_edits[best_edit_site])**(0.5))

    if tissue_combined > 0:
        my_agreement["combined-tissue"]["+"] += 1
    elif tissue_combined < 0:
        my_agreement["combined-tissue"]["-"] += 1
    else:
        my_agreement["combined-tissue"]["na"] += 1
    
    ###################################
    # Test: Combined across tissues, then across sites
    ###################################

    # Stouffer's combined z-score
    z_vec = []
    for site in all_edits:
        tissue_combined = sum(all_edits[site]) / (len(all_edits[site]) ** (0.5))
        z_vec.append(round(tissue_combined,2))
    tissue_combined_site_combined = sum(z_vec) / (len(z_vec) ** (0.5))

    if tissue_combined_site_combined > 0:
        my_agreement["combined-tissue-and-sites"]["+"] += 1
    elif tissue_combined_site_combined < 0:
        my_agreement["combined-tissue-and-sites"]["-"] += 1
    else:
        my_agreement["combined-tissue-and-sites"]["na"] += 1
    
    ###################################
    # Test: Best site, vote across tissues
    ###################################
    
    # Then at that top editing site, count support in favor of either
    # direction, across tissues.
    #
    # For now, we'll only look at that one single best editing site that 
    # we already determined in a previous step
    #

    # In how many tissues does the edQTL and GWAS effect direction match vs. not?
    this_snp = {"+": sum([z > 0 for z in all_edits[best_edit_site]]), "-": sum([z < 0 for z in all_edits[best_edit_site]])}
    
    if this_snp["+"] > this_snp["-"]:
       my_agreement["majority"]["+"] += 1
    elif this_snp["-"] > this_snp["+"]:
        my_agreement["majority"]["-"] += 1
    else:
        my_agreement["majority"]["na"] += 1

    #########################################
    # Test: binomial
    #########################################

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

    # Final test results include all the different test modes
    return my_agreement, tissue_combined 

# Create a set containing every single SNP that was tested for
# edQTL association in at least one tissue
# {(chrom, pos), ....}
# 
# It would be nice to make this faster, but probably not possible because
# you can only read the disk so fast...
def get_eqtl_set(qtls):
    pickle_file = "all_edqtl_{0}.pkl".format(max_edqtl_distance)
    if os.path.exists(pickle_file):
        with open(pickle_file) as f:
            return pickle.load(f)
    snps = set([])
    # For each edQTL tissue...
    for q in qtls:
        print q
        with gzip.open(q) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                # We're going to limit how far away the SNPs can actually be...
                # this improves the reliability of our set of edQTL associations
                if abs(int(data[6])) > max_edqtl_distance:
                    continue
                snps.add((data[8].replace("chr", ""), data[9]))

    with open(pickle_file, "w") as pkl:
        pickle.dump(snps, pkl)

    return snps

# List of tuples representing SNPs that achieved colocalization
# [ (chrom, snp, eqtl_file, gwas_trait, feature) , ...]
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
            if float(data[9]) >= min_coloc_threshold:
                chrom, snp = data[0].split("_")
                # (chrom, snp, eqtl_file, gwas_trait, feature)
                coloc_snps.add((chrom, snp, data[1], data[2], data[3]))
    return coloc_snps

# List all the traits in a gwas file
# by reading from the "trait" column if it exists
def get_traits_in_file(gwas):
        all_traits = set([])
        with gzip.open(gwas) as f:
            head = f.readline().strip().split("\t")
            if "effect_direction" not in head or "effect_allele" not in head or "non_effect_allele" not in head:
                # Signal to skip this trait if the effect direction isn't even obtainable
                return None
            if "trait" in head:
                trait_index = head.index("trait")
                i = 0
                for line in f:
                    # We assume the trait appears at least once in the first 100k lines of the file, since
                    # the file should be sorted by position and not by trait
                    if i > 100000:
                        break
                    data = line.strip().split("\t")
                    all_traits.add(data[trait_index])
                    i += 1
            else:
                all_traits.add(gwas)

        all_traits = list(all_traits)

        return all_traits

if __name__ == "__main__":
    
    main()
