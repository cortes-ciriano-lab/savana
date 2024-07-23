"""
Script to fit absolute copy number
Created: 27/06/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import numpy as np
from scipy.stats import gaussian_kde
import statistics
import math
import copy
import argparse
import timeit
import os
import sys


import savana.cn_functions as cnfitter


def process_allele_counts(allele_counts_bed_path):
    allele_counts = []
    with open(allele_counts_bed_path, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            A,C,G,T,N = int(fields[5]),int(fields[6]),int(fields[7]),int(fields[8]),int(fields[9])
            ps_id = f"{fields[0]}_{fields[-1]}" # assign phase sets IDs as same phase set might accur in different chromosomes
            fields.append(ps_id)
            DP = A+C+G+T+N
            if DP > 20: # only include het SNPs allele counts with a total depth of > 20.
                fields.append(DP)
                # print(fields.append(str(DP)))
                allele_counts.append(fields)
    # columns:"chr", "start", "end","REF", "ALT", "A","C","G","T","N","AF_0","AF_1","GT","PS","ps_id","DP"
    # get unique phase sets to iterate (or multiprocess through)
    phasesets = list(dict.fromkeys([x[-2] for x in allele_counts]))
    file.close()
    return allele_counts, phasesets

# def estimate_cellularity_phased_hetSNPs(phasesets,allele_counts,dp_cutoff,min_ps_size=10,min_ps_length=500000):
def process_phased_hetSNPs(phasesets, allele_counts, dp_cutoff, min_ps_size=10, min_ps_length=500000):
    '''
    Estimate sample cellularity/purity using allele counts of phased hetSNPs.
    '''
    phasesets_dict = {}
    ps_summary = []
    for ps in phasesets:
        if "None" not in ps:
            # print(ps)
            # subset allele count data to ps and remove SNPs with AF0/AF1 of 0 or 1
            ac_ps = [x for x in allele_counts if x[-2] == ps and float(x[10]) != 0 and float(x[10]) != 1 and float(x[11]) != 0 and float(x[11]) != 1]
            ps_snps = len(ac_ps) # number of SNPs in ps
            if ps_snps < min_ps_size: # skip ps with less than min_ps_size SNPs and remove empty ps
                continue
            else:
                ps_length = max([int(x[2]) for x in ac_ps]) - min([int(x[1]) for x in ac_ps]) # genomic lenght of ps_length
                ps_depth = statistics.mean([int(x[-1]) for x in ac_ps]) # mean depth of ps
                # ps_weight = ps_snps/no_hetSNPs # estimate weight of PS based on number of hetSNPs present in PS compared to all hetSNPs across sample
                # all AFs at all het SNP position within ps
                af = [float(x[10]) for x in ac_ps] + [float(x[11]) for x in ac_ps]
                if ps_snps > min_ps_size and ps_length > min_ps_length and ps_depth < dp_cutoff and cnfitter.is_unimodal(af) == False:
                    ps_bc = cnfitter.bimodality_coefficient(af)
                # if ps_snps > min_ps_size and ps_length > min_ps_length and ps_depth < dp_cutoff and is_unimodal(af) == False:
                #     ps_bc = bimodality_coefficient(af)
                    phasesets_dict[ps] = ac_ps
                    ps_summary.append([ps, ps_depth, ps_length, ps_snps, ps_bc])
    return phasesets_dict, ps_summary


def estimate_cellularity(phasesets_dict, ps_summary):
    ps_summary = sorted(ps_summary,key=lambda x: x[-1]) # sort ps by bimodality coefficient
    # ps_summary = sorted(sorted(ps_summary,key=lambda x: x[-2]), key=lambda x: round(x[-1],3)) # sort by ps_bc and weight
    # select top 10% PSs
    # # ps_top = [x[0] for x in ps_summary[-(round(len(ps_summary)*0.01)):]]
    #
    ps_top = [x[0] for x in ps_summary[-10:]]
    # pull out data for ps_top
    ps_top_acs = []
    for pst in ps_top:
        ps_top_acs += phasesets_dict[pst]
    # Estimate cellularity
    af_cutoff = 0.5
    pur0 = statistics.median([float(x[10]) for x in ps_top_acs if float(x[10]) > af_cutoff] + [float(x[11]) for x in ps_top_acs if float(x[11]) > af_cutoff])
    pur1 = statistics.median([float(x[10]) for x in ps_top_acs if float(x[10]) < (1-af_cutoff)] + [float(x[11]) for x in ps_top_acs if float(x[11]) < (1-af_cutoff)])
    cellularity = statistics.mean([1-(1-max(pur0,pur1))*2]+[1-(min(pur0,pur1))*2])
    return cellularity

def process_log2r_input(log2r_cn_path):
    rel_copy_number = []
    with open(log2r_cn_path, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            chrom,start,end,seg_id = fields[1],int(fields[2]),int(fields[3]),fields[-2]
            copy_number=2**float(fields[-1]) # revert log2 of relatitve copy number
            bin_length = end-start+1 
            # print(chrom,start,end,seg_id,copy_number)
            rel_copy_number.append([chrom,start,end,bin_length,seg_id,copy_number])
    # Collapse copy number to segments
    segs = list(dict.fromkeys([x[-2] for x in rel_copy_number])) # unique segment ids
    med_length = statistics.median([x[3] for x in rel_copy_number]) # Estimate median length of all bins to determine weight of segment
    rel_copy_number_segments = []
    for s in segs:
        cn_seg = [x for x in rel_copy_number if x[-2] == s]
        bin_count = len(cn_seg)
        sum_of_bin_lengths = sum([x[3] for x in cn_seg])
        weight = sum_of_bin_lengths/med_length
        min_start = min([x[1] for x in cn_seg])
        max_end = max([x[2] for x in cn_seg])
        # sanity check copy number values for each segment are all the same!
        if all(i[-1] == cn_seg[0][-1] for i in cn_seg) == True:
            rel_copy_number_segments.append([cn_seg[0][0], min_start, max_end, s, bin_count, sum_of_bin_lengths, weight, cn_seg[0][-1]])
        else:
            print(f"    ERROR: segment {s} contains multiple copy number values.")
            break
    return rel_copy_number_segments

def fit_absolute_cn(outdir, log2r_cn_path, allele_counts_bed_path, sample,
    min_ploidy, max_ploidy, ploidy_step, min_cellularity, max_cellularity, cellularity_step, cellularity_buffer,
    distance_function, distance_filter_scale_factor, distance_precision,
    max_proportion_zero, min_proportion_close_to_whole_number, max_distance_from_whole_number,
    min_ps_size, min_ps_length, threads):
    '''
    # 3. Estimate purity using phased heterozygous SNPs
    ### Method/principles based on: https://doi.org/10.1371/journal.pone.0045835
    '''
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... CN fitting will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads

    # Open and prepare data to estimate cellularity if phased hetSNPs allele counts were provided.
    if allele_counts_bed_path == None:
        print("     ... No phased hetSNPs allele counts were provided. Cellularity is estimated using grid search. Note that this is less accurate. Please provide hetSNPs allele counts if possible. See documentation for instructions on how to generate these.")
    elif allele_counts_bed_path != None:
        print("     ... Allele counts for phased hetSNPs provided and being processed to estimate sample purity ...")

        allele_counts, phasesets = process_allele_counts(allele_counts_bed_path)
        # Estimate depth cutoff to exclude potentially amplified/gained regions as those will impact BAF distribution
        ac_mean_depth = statistics.mean([x[-1] for x in allele_counts])
        dp_cutoff = 2 * ac_mean_depth
        # no_hetSNPs = len(allele_counts)

        phasesets_dict, ps_summary = process_phased_hetSNPs(phasesets, allele_counts, dp_cutoff, min_ps_size=min_ps_size, min_ps_length=min_ps_length)

        cellularity = estimate_cellularity(phasesets_dict, ps_summary)
        digs = len(str(cellularity_step))-2 if isinstance(cellularity_step,int) != True else 1
        print(f"        estimated cellularity using hetSNPs = {round(cellularity,digs)}.")
        min_cellularity = round(max(0,cellularity - cellularity_buffer),digs)
        max_cellularity = round(min(1,cellularity + cellularity_buffer),digs)

    #----
    # 4. Estimate ploidy and fit ACN using estimated sample purity
    ### Method/principles based on rascal R package: https://www.biorxiv.org/content/10.1101/2021.07.19.452658v1
    #----

    rel_copy_number_segments = process_log2r_input(log2r_cn_path)

    # Prepare input for copy number fitting
    relative_CN = [x[-1] for x in rel_copy_number_segments]
    weights = [x[-2] for x in rel_copy_number_segments]

    # Copy number fitting
    fits = cnfitter.estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=weights, distance_function=distance_function)
    fits_r = cnfitter.reduce_grid(fits,distance_filter_scale_factor = distance_filter_scale_factor)
    solutions = cnfitter.viable_solutions(fits_r, relative_CN, weights,
                        max_proportion_zero = max_proportion_zero,
                        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
                        max_distance_from_whole_number = max_distance_from_whole_number)

    # check if viable solutions were found. If not, terminate script and write out error message and arguments to file for inspection and adjustment
    if len(solutions) == 0:
        print("No fits found. See No_fit_found_PARAMS_out.tsv in output")
        with open(f"{outdir}/No_fit_found_PARAMS.tsv", 'w') as params_out:
            params_out.write(f'No viable solution fullfilling set paramaters was found. \nPerform QC, review parameters and rerun if required/appropriate with adjusted parameters.\n')
            # TODO: add this back in (either manually or one step up)
            """
            params_out.write(f'\nPARAMETERS:\n')
            for key, value in vars(args).items():
                    params_out.write(f'{key}: {value}\n')
            """
            params_out.write(f'\nCandidate fits found prior to checking viability/acceptability using set parameters:\n')
            header=['purity','ploidy','x','y','distance']
            params_out.write('\t'.join(header)+'\n')
            for r in fits_r:
                Line = '\t'.join(str(e) for e in r) + '\n'
                params_out.write(Line)
        sys.exit(1) # Exit the script with a status code of 1 (indicating an error)

    solutions_ranked = cnfitter.rank_solutions(solutions,distance_precision)
    final_fit = solutions_ranked[0]
    print(f"        Data fitted to purity = {final_fit[0]} and ploidy = {final_fit[1]}.")


    # Convert relative to absolute copy number and obtain major and minor copy number values if allele counts provided
    if allele_counts_bed_path == None:
        print("     ... No phased hetSNPs allele counts provided. Only total absolute copy number will be estimated. Consider providing hetSNPs allele counts if possible. See documentation for instructions on how to generate these.")
        fitted_purity, fitted_ploidy = final_fit[0], final_fit[1]
        abs_copy_number_segments = copy.deepcopy(rel_copy_number_segments)
        for x in abs_copy_number_segments:
            acn = round(cnfitter.relative_to_absolute_CN(x[-1], fitted_purity, fitted_ploidy),4)
            x[-1] = acn if acn > 0 else 0
        # set negative values to 0
        for x in abs_copy_number_segments:
            if x[-1] < 0:
                x[-1] = 0

    elif allele_counts_bed_path != None:
        print("     ... Allele counts for phased hetSNPs provided. Minor and total absolute copy number are being estimated ...")
        fitted_purity, fitted_ploidy = final_fit[0], final_fit[1]
        # prepare for multiprocessing
        chr_names = list(dict.fromkeys([x[0] for x in rel_copy_number_segments]))
        # only use multiprocessing if more than 1 thread available/being used.
        if threads == 1:
            # loop through chromosomes
            print("multithreading skipped.")
            abs_copy_number_segments = []
            for chrom in chr_names:
                cn_chr = [x for x in rel_copy_number_segments if x[0] == chrom]
                ac_chr = [x for x in allele_counts if x[0] == chrom]
                cur_chr = cnfitter.relative_to_absolute_minor_total_CN(chrom, cn_chr, ac_chr, fitted_purity, fitted_ploidy)
                abs_copy_number_segments.append(cur_chr)
            abs_copy_number_segments = [x for xs in abs_copy_number_segments for x in xs]

        else:
            print(f"multithreading using {threads} threads.")
            args_in = [[chrom, [x for x in rel_copy_number_segments if x[0] == chrom], [x for x in allele_counts if x[0] == chrom], fitted_purity, fitted_ploidy] for chrom in chr_names]
            with Pool(processes=threads) as pool:
                abs_copy_number_segments = [x for xs in list(pool.starmap(cnfitter.relative_to_absolute_minor_total_CN, args_in)) for x in xs]


    # Prepare and write out results
    ## ranked solutions, final fit, and converted segmented absolute copy number

    outfile1 = open(f"{outdir}/{sample}_ranked_solutions.tsv", "w")
    header=['purity','ploidy','distance','rank']
    outfile1.write('\t'.join(header)+'\n')
    for r in solutions_ranked:
        Line = '\t'.join(str(e) for e in r) + '\n'
        outfile1.write(Line)
    outfile1.close()

    outfile2 = open(f"{outdir}/{sample}_fitted_purity_ploidy.tsv", "w")
    outfile2.write('\t'.join(header)+'\n')
    Line = '\t'.join(str(e) for e in final_fit) + '\n'
    outfile2.write(Line)
    outfile2.close()

    outfile3 = open(f"{outdir}/{sample}_segmented_absolute_copy_number.tsv", "w")
    if allele_counts_bed_path == None:
        header=['chromosome','start','end','segment_id', 'bin_count', 'sum_of_bin_lengths', 'weight', 'copyNumber']
    elif allele_counts_bed_path != None:
        header=['chromosome','start','end','segment_id', 'bin_count', 'sum_of_bin_lengths', 'weight', 'copyNumber', 'minorAlleleCopyNumber']
    outfile3.write('\t'.join(header)+'\n')
    for r in abs_copy_number_segments:
        Line = '\t'.join(str(e) for e in r) + '\n'
        outfile3.write(Line)
    outfile3.close()

    # TODO: add this back in - either manually or one step up
    """
    with open(f"{outdir}/PARAMS_out.tsv", 'w') as params_out:
            params_out.write('PARAMETERS:\n')
            for key, value in vars(args).items():
                    params_out.write(f'{key}: {value}\n')
    """

    ###############
    # Add in code to redefine segment breaks based on rounded acn?
    ###############