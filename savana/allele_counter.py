"""
Script to to count heterozygous SNPs
(required prior to copy number fitting for purity estimation)
Created: 24/06/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import cyvcf2
import pysam
import pybedtools
import os
import glob
import math


def extract_hets(snp_vcf, g1000_vcf, FASTA, window):
    '''
    Extracts heterozygous SNPs from phased.vcf (i.e. from matched normal bam).
    '''
    if snp_vcf != None and g1000_vcf == None:
        print(f"    ... Extracting heterozygous SNPs from {snp_vcf} ... ")
        vcf_reader = cyvcf2.VCF(snp_vcf)
        hets_dict = {}
        # iterate through variants
        for variant in vcf_reader:
            # check if variant is a SNP
            if variant.is_snp:
                # iterate through each genotype
                for s_idx, gt in enumerate(variant.genotypes):
                    # only get heterozygous snps
                    if gt[0] != gt[1]:
                        # ps = int(variant.format('PS')[s_idx]) if 'PS' in variant.FORMAT else None
                        gt_str = f"{gt[0]}|{gt[1]}" if gt[2] == 1 else f"{gt[0]}/{gt[1]}"
                        CHROM=variant.CHROM
                        key=f"{str(variant.POS)}_{str(variant.REF)}"
                        # var_out = [variant.ALT[0],gt_str,str(ps)]
                        # nested dictionary - by chromosome and by position (in 100k increments)
                        pos_cat = math.floor(int(variant.POS)/window) * window
                        var_out = [variant.ALT[0],f'{CHROM}_{pos_cat}']
                        if (CHROM not in hets_dict):
                            hets_dict[CHROM] = {}
                            if (pos_cat not in hets_dict[CHROM]):
                                hets_dict[CHROM][pos_cat] = {}                          
                                hets_dict[CHROM][pos_cat][key] = var_out
                        elif (pos_cat not in hets_dict[CHROM]):
                            hets_dict[CHROM][pos_cat] = {} 
                            hets_dict[CHROM][pos_cat][key] = var_out
                        else: 
                            hets_dict[CHROM][pos_cat][key] = var_out
        print(f"    ... Heterozygous SNPs from {snp_vcf} extracted. Extracting allele counts for heterozygous SNPs ...")
    elif snp_vcf == None and g1000_vcf != None:
        #check for chromosome annotation in FASTA
        inFasta = pysam.FastaFile(FASTA)
        ref_contigs = inFasta.references
        inFasta.close()
        chr_annot = True if 'chr1' in ref_contigs else False
        vcf_dir = os.path.join(os.path.dirname(__file__),'1K_genome_vcf')
        if g1000_vcf == "1000g_hg38":
            #vcf_path = os.path.join(vcf_dir, 'hg38_g1000_biallelic_AF0.35-0.65.vcf.gz')
            vcf_path = os.path.join(vcf_dir, 'hg38_g1000_biallelic_AF0.25-0.75.vcf.gz')
        if g1000_vcf == "1000g_hg19":
            #vcf_path = os.path.join(vcf_dir, 'hg19_g1000_biallelic_AF0.35-0.65.vcf.gz')
            vcf_path = os.path.join(vcf_dir, 'hg19_g1000_biallelic_AF0.25-0.75.vcf.gz')
        if g1000_vcf == "1000g_t2t":
            vcf_path = os.path.join(vcf_dir, 't2t_g1000_biallelic_AF0.25-0.75.vcf.gz')
        vcf_reader = cyvcf2.VCF(vcf_path)
        hets_dict = {}
        # iterate through variants
        for variant in vcf_reader:
            # check if variant is a SNP
            if variant.is_snp:
                CHROM=variant.CHROM
                if chr_annot == True:
                    CHROM=f'chr{CHROM}'
                key=f"{str(variant.POS)}_{str(variant.REF)}"
                # var_out = [variant.ALT[0],gt_str,str(ps)]
                # nested dictionary - by chromosome and by position (in 100k increments)
                pos_cat = math.floor(int(variant.POS)/window) * window
                var_out = [variant.ALT[0],f'{CHROM}_{pos_cat}']
                if (CHROM not in hets_dict):
                    hets_dict[CHROM] = {}
                    if (pos_cat not in hets_dict[CHROM]):
                        hets_dict[CHROM][pos_cat] = {}                          
                        hets_dict[CHROM][pos_cat][key] = var_out
                elif (pos_cat not in hets_dict[CHROM]):
                    hets_dict[CHROM][pos_cat] = {} 
                    hets_dict[CHROM][pos_cat][key] = var_out
                else: 
                    hets_dict[CHROM][pos_cat][key] = var_out
        print(f"    ... Heterozygous SNPs from {g1000_vcf} extracted. Extracting allele counts for heterozygous SNPs ...")
    return hets_dict

def process_allele_counts(allele_counts_path, hets_dict, window):
    '''
    Extracts and process allele counts for heterozygous SNPs dictionary.
    '''
    ac_bed = pybedtools.BedTool(allele_counts_path)
    bed_out = []
    for site in ac_bed:
        chrom,pos,ref,dp= site[0], site[2], site[3], site[9]
        bases = { 'A': site[4], 'C': site[5], 'G': site[6], 'T': site[7], 'N': site[8]}
        key=f"{pos}_{ref}"
        pos_cat = math.floor(int(pos)/window) * window # same key as above
        # try:
        #     hets_dict[chrom][pos_cat][key]
        #     pass
        # except:
        #     continue
        if chrom not in hets_dict.keys():
            continue
        if pos_cat not in hets_dict[chrom].keys():
            continue
        if key not in hets_dict[chrom][pos_cat].keys(): 
            continue
        key_out=hets_dict[chrom][pos_cat][key]
        # alt,gt,ps=key_out[0],key_out[1],key_out[2]
        alt,block_id=key_out[0],key_out[1]
        # estimate AFs
        AF_0 = float(bases[ref]) / float(dp)
        AF_1 = float(bases[alt]) / float(dp)
        # define output list
        out = [chrom, pos, pos, ref, alt, bases['A'],bases['C'],bases['G'],bases['T'],bases['N'], str(AF_0), str(AF_1), block_id]
        bed_out.append(out)
    return bed_out

# the following five functions were adapted and taken from Fran Muyas
def collect_result(result):
    if (result[1] != ''):
        VARIANTS = result[1]
        OUT = result[0]
        out = open(OUT,'w')
        out.write(VARIANTS)
        out.close()

def concatenate_sort_temp_files_and_write(out_file, tmp_dir):
    # Get file paths of temp files
    print("    ... Concatenating temp files into interim allele count results ... ")
    all_files = glob.glob(tmp_dir + '/*.hetsnps_allelecounts.temp')
    # Load temp files
    if (len(all_files) > 0):
        # Organise files in dictionaries of chromosomes and start positions
        Dictionary_of_files = {}
        for filename in all_files:
            basename = os.path.basename(filename)
            basename = basename.split(".")
            # positions
            coordinates = basename[-3]
            CHROM, START, END = coordinates.split("__")
            START = int(START)
            if (CHROM not in Dictionary_of_files):
                Dictionary_of_files[CHROM] = {}
                Dictionary_of_files[CHROM][START] = filename
            else:
                Dictionary_of_files[CHROM][START] = filename
        ## Write to final output file
        out = open(out_file,'w')
        # Header = ['#CHROM', 'Start', 'End','REF', 'A', 'C', 'G', 'T', 'N', 'DP']
        # Header = '\t'.join(Header) + '\n'
        # out.write(Header)
        # Move through filenames by sorted coordinates
        for chrom in sorted(Dictionary_of_files.keys()):
            for start in sorted(Dictionary_of_files[chrom].keys()):
                filename = Dictionary_of_files[chrom][start]
                with open(filename, 'r') as f:
                    out.write(f.read())
                    out.write('\n')
                # Remove temp file
                os.remove(filename)
        out.close()
    else:
        # If temp files not found, print message
        print ('No temporary files found')

def MakeWindows(CONTIG, FASTA, window):
    inFasta = pysam.FastaFile(FASTA)	
    # Generate bed file based on all coordenates from reference fasta file
    # only allow canonical contigs -- might update to parameters later... 
    contigs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    ref_contigs = inFasta.references
    if 'chr1' in ref_contigs:
        contigs = [f'chr{x}' for x in contigs]   
    if (CONTIG != 'all'):
        CONTIG_Names = [contigs[(int(x)-1)] for x in CONTIG]
    else:
        CONTIG_Names = contigs
    # print(CONTIG_Names)
    LIST = [ (x, 1, inFasta.get_reference_length(x)) for x in CONTIG_Names]
    # print(LIST)
    a = pybedtools.BedTool(LIST)
    # final bed to be used for allele counting with windows
    final_bed = a.window_maker(a,w=window)
    inFasta.close()
    return(final_bed)

def BaseCount(LIST):
    Bases=['A','C','G','T','N']
    # Dictinary with base counts
    NUCLEOTIDES = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    # Update dictionary counts
    for x in LIST:
        if x.upper() in Bases:
            NUCLEOTIDES[x.upper()] += 1
    return (NUCLEOTIDES)

def run_interval(interval, BAM, FASTA, MQ, MIN_COV, tmp_dir):
    # Coordinates to analyse
    CHROM  =  interval[0]
    START = int(interval[1])
    END = int(interval[2])
    # Get pileup read counts from coordinates
    bam = pysam.AlignmentFile(BAM)
    i = bam.pileup(CHROM, START, END, min_mapping_quality = MQ, ignore_overlaps = False)
    # Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
    inFasta = pysam.FastaFile(FASTA)
    # Run it for each position in pileup
    POSITIONS = []
    for p in i:
        POS=p.pos
        if POS >= START and POS < END:
            # Get reference base from fasta file
            ref_base = inFasta.fetch(CHROM, POS, POS+1).upper()
            # Get coverage
            DP = p.get_num_aligned()
            # Run only if coverage is more than minimum (arg parameter)
            # if (DP >= MIN_COV and ref_base != 'N'):
            if (DP >= MIN_COV and ref_base in ['A','C','G','T']):
                # Get pileup list
                PILEUP_LIST = p.get_query_sequences(mark_matches=True, add_indels=True)
                BASE_COUNTS = BaseCount(PILEUP_LIST)
                # only include output if ALT != REF
                total_count = sum(BASE_COUNTS.values())
                if BASE_COUNTS[ref_base] + BASE_COUNTS['N'] < total_count:
                    # Pysam provides a 0-based coordinate. We must sum 1 to this value
                    POS_print = POS + 1
                    # Lines to print
                    LINE_0 = '\t'.join([str(BASE_COUNTS['A']), str(BASE_COUNTS['C']), str(BASE_COUNTS['G']), str(BASE_COUNTS['T']), str(BASE_COUNTS['N'])])  
                    LINE_1 = '\t'.join([str(CHROM), str(POS_print), str(POS_print), str(ref_base), LINE_0 , str(total_count)])
                    # Save results in list of positions
                    POSITIONS.append(LINE_1)   
    inFasta.close()
    bam.close()
    # Return list of positions
    ID = '__'.join([str(CHROM), str(START), str(END)])
    out_temp = tmp_dir + '/' + ID + '.hetsnps_allelecounts.temp'
    return([out_temp,'\n'.join(POSITIONS)])

#----
# Prepare input and run allele counter across heterozygous SNPs
#----
def perform_allele_counting(outdir, sample, contigs, fasta_file_path, snp_vcf, g1000_vcf, tumour, ac_window, allele_mapq, allele_min_reads, tmp_dir, threads):
    """ extract allele counts for heterozygous SNPs """
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Allele counter will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads
    #----
    # 1. Create bed file and windows for allele counting
    #----
    print(f"    ... splitting helper bed into windows for allele counting ...")
    bed = MakeWindows(contigs, fasta_file_path, ac_window)
    #----
    # 2. Multiprocess and run bed bins in parallel
    #----
    print(f"    ... Allele counting ...")
    if (threads > 1):
        pool = Pool(threads)
        # Use loop to parallelize
        for row in bed:
            # This funtion writes out temp files
            pool.apply_async(run_interval, args=(row, tumour, fasta_file_path, allele_mapq, allele_min_reads, tmp_dir), callback=collect_result)
        # Close Pool and let all the processes complete    
        pool.close()
        pool.join()
        print("         ... interim finished...")
    else:
        for row in bed:
            # This funtion writes in temp files the results
            collect_result(run_interval(row, tumour, fasta_file_path, allele_mapq, allele_min_reads, tmp_dir))
    #----
    # 3. Get interim results and write out
    #----
    # print(f"    ... Concatenating temp files into interim allele count results ... ")
    out_path_interim = f"{outdir}/{sample}_allele_counts_hetSNPs_INTERIM.bed"
    concatenate_sort_temp_files_and_write(out_path_interim, tmp_dir)
    #----
    # 4. Generate dictionary of heterozygous SNPs from normal VCF
    #----
    # print(f"    ... Extracting heterozygous SNPs from {phased_vcf} ... ")
    het_snps = extract_hets(snp_vcf, g1000_vcf, fasta_file_path, ac_window)
    # print(f"    ... Heterozygous SNPs from {phased_vcf} extracted. Extracting allele counts for heterozygous SNPs ...")
    #----
    # 5. Extract allele counts for hetSNPs
    #----
    allele_counts_final = process_allele_counts(out_path_interim, het_snps, ac_window)
    #----
    # 6. Get results and write out
    #----
    out_path = f"{outdir}/{sample}_allele_counts_hetSNPs.bed"
    outfile = open(out_path, "w")
    for row in allele_counts_final:
        Line = '\t'.join(row) + '\n'
        outfile.write(Line)
    outfile.close()
    # remove interim result bed file
    os.remove(out_path_interim)

    return out_path

if __name__ == "__main__":
    print("Allele counter")
