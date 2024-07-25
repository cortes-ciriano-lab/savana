"""
Script to to count heterozygous SNPs
(required prior to copy number fitting for purity estimation)
Created: 24/06/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import argparse
import timeit
import cyvcf2
import pysam
import pybedtools
import math
import os

def extract_hets(phased_vcf):
    '''
    Extracts heterozygous SNPs from phased.vcf (e.g. from matched normal bam).
    '''
    print(f"    Extracting phased heterozygous SNPs from {phased_vcf} ...")
    vcf_reader = cyvcf2.VCF(phased_vcf)
    hets = []
    # iterate through variants
    for variant in vcf_reader:
        # check if variant is a SNP
        if variant.is_snp:
            # iterate through each genotype
            for s_idx, gt in enumerate(variant.genotypes):
                # only get heterozygous snps
                if gt[0] != gt[1]:
                    # print(s_idx, gt)
                    # sample = vcf_reader.samples[s_idx]
                    ps = int(variant.format('PS')[s_idx]) if 'PS' in variant.FORMAT else None
                    gt_str = f"{gt[0]}|{gt[1]}" if gt[2] == 1 else f"{gt[0]}/{gt[1]}"
                    # dp = int(variant.format('DP')[s_idx]) if 'DP' in variant.FORMAT else None
                    # ad = int(variant.format('AD')[s_idx]) if 'AD' in variant.FORMAT else None
                    var_out = [variant.CHROM,str(variant.POS),str(variant.POS),variant.REF,variant.ALT[0],gt_str,str(ps)]
                    hets.append(var_out)
    return hets

def chunkify_bed(bed_file, chunk_size):
    '''
    Divides bed file into chunks based on threads available for multiprocessing.
    '''
    chunks = []
    current_chunk = []
    for i, feature in enumerate(bed_file):
        current_chunk.append(feature)
        if (i + 1) % chunk_size == 0:
            chunks.append(pybedtools.BedTool(current_chunk))
            current_chunk = []
    if current_chunk:
        chunks.append(pybedtools.BedTool(current_chunk))
    return chunks

def allele_count(curr_chunk, sites, bam, allele_mapq, allele_min_reads):
    '''
    Count alleles across heterozygous SNP positions from phased normal bam
    '''
    print(f"    Read counting {curr_chunk} ...")
    # open tumour bam
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        list_out = []
        for site in sites:
            # print(site)
            dict_out = {}
            # if len(site[0])==len(chr_now):
            #     chr = site[0];#[3];
            # else:
            #     chr = site[0]#[3:5];
            # if chr != chr_now:
            #     continue
            # chr2    61791980    61791980    G    A    1|0    58907085
            chrom = site[0]
            start = int(site[1])-1
            end = int(site[2])  # start and end cannot be the same position or otherwise the fetching does not work
            ref = site[3]
            alt = site[4]
            GT = site[5]
            PS = site[6]

            reads = [
                read for read in bamfile.fetch(
                    chrom,
                    start,
                    end,
                    multiple_iterators=True
                )
                ]
            # filter reads
            reads = [
                read for read in reads if
                read.mapping_quality >= allele_mapq and not read.is_secondary and not read.is_supplementary
                ]
            # perform base counting
            try:
                if len(reads) >= allele_min_reads:
                    key_now = f"{chr}_{start}_{end}_{ref}_{alt}"
                    # key_now = "{}\t{}\t{}\t{}\t{}".format(chr, start, end, ref ,alt)
                    for read in reads:
                        read_sequence = read.seq

                        aligned_pos = read.get_reference_positions(
                            full_length=True
                        )
                        try:
                            idx = aligned_pos.index(start)
                            #DP+=1
                            snp_read = read_sequence[idx]
                            if key_now in dict_out: #.has_key(key_now):
                                dict_out[key_now][snp_read] = dict_out[key_now][snp_read] + 1
                            else:
                                dict_out[key_now]={"A":0, "C":0, "G":0, "T":0, "N":0}
                                dict_out[key_now][snp_read] = dict_out[key_now][snp_read] + 1
                        except:
                            continue
                    DP = dict_out[key_now]["A"] + dict_out[key_now]["C"] + dict_out[key_now]["G"] + dict_out[key_now]["T"]
                    if GT == "0|1":
                        AF_0 = dict_out[key_now][ref] / DP
                        AF_1 = dict_out[key_now][alt] / DP
                    if GT == "1|0":
                        AF_0 = dict_out[key_now][alt] / DP
                        AF_1 = dict_out[key_now][ref] / DP
                    #print("{} {}".format(AF_0, AF_1))
                    dict_out[key_now]["AF_0"] = AF_0
                    dict_out[key_now]["AF_1"] = AF_1

                    out = [chrom, str(start), str(end), ref, alt, str(dict_out[key_now]["A"]), str(dict_out[key_now]["C"]), str(dict_out[key_now]["G"]), str(dict_out[key_now]["T"]), str(dict_out[key_now]["N"]) , str(dict_out[key_now]["AF_0"]), str(dict_out[key_now]["AF_1"]), GT, str(PS)]

                list_out.append(out)
            except:
                continue
    return(list_out)


#----
# 3. Prepare input and run allele counter across heterozygous SNPs
#----
def perform_allele_counting(outdir, sample, phased_vcf, tumour, allele_mapq, allele_min_reads, threads):
    """ extract heterozygous SNPs """
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Allele counter will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads
    hets_bed_path = f"{outdir}/{sample}_phased_het_snps.bed"
    outfile = open(hets_bed_path, "w")
    het_snps = extract_hets(phased_vcf)
    for r in het_snps:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    ## RUN IN CHUNKS TO SPEED UP
    # Define bed_file chunks
    bed_file =  pybedtools.BedTool(hets_bed_path)
    bed_length = sum(1 for _ in bed_file)
    chunk_size = bed_length // threads + (bed_length % threads > 0)
    print(f"    splitting bed file into n = {math.ceil(bed_length / chunk_size)} chunks ...")
    # Split the BED file into chunks
    chunks = chunkify_bed(bed_file, chunk_size)

    # only use multiprocessing if more than 1 thread available/being used.
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        allele_counts = []
        for idx,bed_chunk in enumerate(chunks):
            curr_chunk = f"chunk {idx+1}"
            allele_counts_chunk = allele_count(curr_chunk, bed_chunk, tumour, allele_mapq, allele_min_reads)
            allele_counts.append(allele_counts_chunk)
        allele_counts = [x for xs in allele_counts for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,tumour,allele_mapq,allele_min_reads] for idx,bed_chunk in enumerate(chunks)]
        with Pool(processes=threads) as pool:
            allele_counts = [x for xs in list(pool.starmap(allele_count, args_in)) for x in xs]

    #----
    # 4. Get results and write out
    #----
    out_path = f"{outdir}/{sample}_allele_counts_hetSNPs.bed"
    outfile = open(out_path, "w")
    for r in allele_counts:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    return out_path

if __name__ == "__main__":
    print("Allele counter")
