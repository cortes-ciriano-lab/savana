"""
Script to generate binned reference for copy number estimation
Created: 31/08/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import pybedtools
import pysam
import gzip

#----
# Define functions
#----
def overlaps(a, b):
    '''
    Function to estimate overlaps between two regions
    '''
    return min(a[1], b[1]) - max(a[0], b[0]) + 1

def parse_vcf(breakpoints_file, bp_dict):
    '''
    Generate dictionary from Savana output vcf with SV breakpoint positions for each chromosomes
    '''
    # Open the Savana VCF file containing SVs. (If it is gzipped, use gzip.open, otherwise use open)
    open_vcf = gzip.open if breakpoints_file.endswith('.gz') else open
    with open_vcf(breakpoints_file, 'rt') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            fields = line.split('\t')
            chr = fields[0]  # Chromosome
            pos = int(fields[1])  # Start position
            if chr in bp_dict:
                bp_dict[chr].append(pos)
    return bp_dict

def process_chromosome(chrom, chr_length, fasta_file_path, bin_size, blacklist, breakpoints = []):
    '''
    Function to process each chromosome of interest (either in a loop or in parallel) to generate bins based on bin_size with annotated blacklist and unknown bases values. 
    If a vcf file (e.g. with SV breakpoints from Savana) is provided, bins will be adjusted to those breakpoints.
    '''
    print(f"    Processing {chrom} ...")

    fasta_file = pysam.FastaFile(fasta_file_path)
    chr_starts = list(range(1, chr_length + 1, bin_size))

    # append chr_starts with vcf breakpoint positions if provided
    if breakpoints is not None:
        print(f"    adding vcf breakpoints for {chrom} ...")
        for bp in breakpoints:
            chr_starts.append(bp)
        # make sure all bin starts are unique (in case breakpoints = chr_start position) to avoid bins of 0kbp bin size.
        chr_starts = list(set(chr_starts))
        chr_starts.sort()

    # define bins
    chr_bins = []
    for i in range(len(chr_starts)-1):
        bin_start,bin_end = chr_starts[i], chr_starts[i+1]-1
        chr_bins.append([bin_start,bin_end])
    # add last bin
    chr_bins.append([chr_starts[-1],chr_length])

    # estimate overlap with blacklist if blacklist is provided
    if (blacklist is not None):
        blacklist_chr = pybedtools.BedTool(blacklist).filter(lambda b: b.chrom == chrom)
        blacklist_chr_list = [[int(row[1]), int(row[2])] for row in blacklist_chr]
    else:
        pass

    # Generate output
    List_of_bins = []
    for bin in chr_bins:
        start, end = bin
        cur_bin_size = float(end)-float(start)+1

        # print(chrom,bin,cur_bin_size)
        # if cur_bin_size == 0:
        #     print(f"    ...ERROR! cur_bin_size = {cur_bin_size} ")

        # Estimate overlap with blacklist for each bin
        # Using overlap function
        if blacklist is not None:
            overlap_list = [overlaps(bin,row) for row in blacklist_chr_list]
            overlap_list = [x for x in overlap_list if x > 0]
            overlap_pct = sum(overlap_list) / cur_bin_size * 100
            # print(overlap_pct)

        # Estimate overlap manually
        # bl_bin_bases = 0
        # for row in blacklist_chr_list:
        #     bed_start, bed_end = row
        #     if bed_start >= start and bed_start <= end:
        #         bl_length = min(bed_end, end) - bed_start + 1
        #         bl_bin_bases += bl_length
        #     elif bed_end >= start and bed_end <= end:
        #         bl_length = bed_end - min(bed_start,start) + 1
        #         bl_bin_bases += bl_length
        # overlap_pct = bl_bin_bases / float(cur_bin_size) * 100

        # Estimate number/fraction of known bases for each bin
        chr_bin_seq = fasta_file.fetch(chrom, start, end)
        acgt_len = 0
        for base in chr_bin_seq:
            if base.upper() in 'ACTG':
                acgt_len += 1
        chr_bases = acgt_len / cur_bin_size * 100

        if blacklist is not None:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases), str(overlap_pct)])
        else:
            List_of_bins.append([chrom, str(start), str(end), str(chr_bases)])

    fasta_file.close()
    return(List_of_bins)

def generate_bins(outdir, sample, ref, chromosomes, bin_size, blacklist, breakpoints, threads):
    '''
    Main function to process and bin fasta file
    '''
    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Bin generator will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads
    bin_size = bin_size * int(1000) # in kb
    # check if vcf file was provided
    if breakpoints is not None:
        binmode = "withVCF"
        print(f"... using vcf file: '{breakpoints}' in addition to binning...")
    else:
        binmode = "standard"
    # a. Extract chormosome info and pass fasta file into chr_in list
    fasta = pysam.FastaFile(ref)
    if chromosomes != 'all':
        chr_names = [fasta.references[(int(x)-1)] for x in chromosomes]
    else:
        chr_names = fasta.references[0:24]
    # b. Extract breakpoints into dictionary from savana vcf to integrate into CN analysis if provided, and define input for main processing step
    if (breakpoints is not None):
        in_dict = { chrom : list([]) for chrom in chr_names }
        bp_dict = parse_vcf(breakpoints, in_dict)
        chr_in = [[chrom,fasta.get_reference_length(chrom),ref,bin_size,blacklist,bp_dict[chrom]] for chrom in chr_names]
    else:
        chr_in = [[chrom,fasta.get_reference_length(chrom),ref,bin_size,blacklist] for chrom in chr_names]
    fasta.close()

    # c. Generate bins and count bases per bins for each chrom in parallel (run function using multithreader)
    ## only use multiprocessing if more than 1 thread available/being used.
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        chrData = []
        for contig in chr_in:
            chrom, chr_length, fasta_file_path = contig[0], contig[1], contig[2]
            if (breakpoints is not None):
                binned_chr = process_chromosome(chrom, chr_length, fasta_file_path, bin_size, blacklist, bp_dict[chrom])
            elif (breakpoints is None):
                binned_chr = process_chromosome(chrom, chr_length, fasta_file_path, bin_size, blacklist)
            chrData.append(binned_chr)
    else:
        print(f"multithreading using {threads} threads.")
        # with ThreadPoolExecutor(max_workers=24) as executor:
        # chrData = list(executor.map(process_chromosome, chr_in))
        with Pool(processes=threads) as pool:
            chrData = list(pool.starmap(process_chromosome, chr_in))

    # c. Concat results into a single bed file
    if chromosomes != 'all':
        outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_subset_{sample}{binmode}.bed"
    else:
        outfile_name = f"{outdir}/{int(bin_size/1000)}kbp_bin_ref_all_{sample}{binmode}.bed"

    outfile = open(outfile_name, "w")

    for obj in chrData:
        for r in obj:
            Line = '\t'.join(r) + '\n'
            outfile.write(Line)

    outfile.close()

    return outfile_name

if __name__ == "__main__":
    print("Bin generator")