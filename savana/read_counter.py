"""
Script to count reads copy number estimation
Created: 12/09/2023
Python 3.9.7
Carolin Sauer
"""

from multiprocessing import Pool
from multiprocessing import cpu_count
import pysam
import pybedtools
import copy
import statistics
import math


def count_reads_in_curr_bin(bam, chrom, start, end, readcount_mapq):
    chunk_read_count = 0
    for read in bam.fetch(chrom, start, end, multiple_iterators = True):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < readcount_mapq:
            continue
        else: 
            read_start = read.reference_start
            read_end = read.reference_end
            if read_start >= start and read_start <= end:
                    chunk_read_count += 1
            elif read_end >= start and read_end <= end:
                    chunk_read_count += 1
    return chunk_read_count

def binned_read_counting(curr_chunk, bed_chunk, alns, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq):
    print(f"    Read counting {curr_chunk} ...")

    # Open the BAM file using pysam
    bam_T = pysam.AlignmentFile(alns['tumour'], "rb")
    bam_N = pysam.AlignmentFile(alns['normal'], "rb") if nmode == "mnorm" and len(alns) == 2 else None

    # Read the BED file to get regions and start/end positions for each chromosome
    # bed_file =  pybedtools.BedTool(bed).filter(lambda b: b.chrom == chr)

    chr_read_counts = []
    for bin in bed_chunk:
        chrom, start, end, bases = bin[0], int(bin[1]), int(bin[2]), float(bin[3])
        bin_name = f"{chrom}:{start}_{end}"
        # Filtering of bins by assigning each bin to boolean 'use' variable
        if blacklisting == False:
            if bases_filter == False:
                use = True
            elif bases_filter == True:
                if bases >= bases_threshold:
                    use = True
                else:
                    use = False

        elif blacklisting == True:
            blacklist = float(bin[4])
            if bases_filter == False:
                if blacklist <= bl_threshold:
                    use = True
                else:
                    use = False
            elif bases_filter == True:
                if blacklist <= bl_threshold and bases >= bases_threshold:
                    use = True
                else:
                    use = False
        # counting reads in tumour (and if provided) normal bam files in each bin
        chunk_read_count_T = count_reads_in_curr_bin(bam_T, chrom, start, end, readcount_mapq)
        chunk_read_count_N = count_reads_in_curr_bin(bam_N, chrom, start, end, readcount_mapq) if bam_N else None

        # if nmode == "mnorm" and len(aln_files) == 2:
        #     chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        # elif nmode != "mnorm" and len(aln_files) == 1:
        #     chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T)])
        if blacklisting == True:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(blacklist), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])
        else:
            chr_read_counts.append([bin_name, chrom, str(start), str(end), str(bases), str(use), str(chunk_read_count_T), str(chunk_read_count_N)])

    # Close BAM file and return out
    bam_T.close()
    if bam_N is not None:
        bam_N.close()
    return(chr_read_counts)

def filter_and_normalise(nmode, countData):
    """ perform filtering, normalisation and transformation """
    ## filtering: Remove bins with 0 reads (no sequencing data in this region). - need to be excluded for segmentation and log2 transformation
    if nmode == "self":
        filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0]
        med_self = statistics.median([int(x[-2]) for x in filtered_counts]) #estimate genome wide median for selfnormalisation
        # Normalise and log2 transform
        normalised_counts = copy.deepcopy(filtered_counts)
        for r in normalised_counts:
            r[-2] = str(math.log2(int(r[-2])/med_self)) # median normalise readcounts
            del r[-1]
    # elif nmode == "pon":
    #     print("This function is yet to be implemented...")
    elif nmode == "mnorm":
        filtered_counts = [x for x in countData if x[-3] == 'True' and int(x[-2]) != 0 and int(x[-1]) != 0]
        cov_scaler = statistics.median([math.log2(int(x[-2])/int(x[-1])) for x in filtered_counts])
        normalised_counts = copy.deepcopy(filtered_counts)
        for r in normalised_counts:
            n = str(math.log2(int(r[-2])/int(r[-1])) - cov_scaler)
            del r[-2:]
            r.append(n)
    # return out
    return filtered_counts, normalised_counts

def chunkify_bed(bed_file, chunk_size):
    '''
    Divides bed files into chunks based on threads available for multiprocessing.
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

def count_reads(outdir, tumour, normal, sample, bin_annotations_path, readcount_mapq, blacklisting, bl_threshold, bases_filter, bases_threshold, threads):
    """ Perform binned read counting on bam file/files per chromosome """

    # check and define threads
    new_threads = min(threads, cpu_count())
    print(f"... Bin read counter will use threads = {new_threads}. (threads = {threads} defined; threads = {cpu_count()} available) ...")
    threads = new_threads

    if normal is not None:
        nmode = "mnorm" #matched normal will be used for logR normalisation
        aln_files = {
                'tumour': tumour,
                'normal': normal
            }
    # elif panel_of_normals is not None:
    #     nmode = "pon" #panel of normals will be used for logR normalisation
    #     aln_files = {
    #             'tumour': tumour
    #         }
    elif normal is None: #and panel_of_normals is None:
        nmode = "self"
        aln_files = {
                'tumour': tumour
            }

    print(f"Normalisation mode: {nmode}")

    # Define bed_file chunks
    bed_file =  pybedtools.BedTool(bin_annotations_path)
    bed_length = sum(1 for _ in bed_file)
    chunk_size = bed_length // threads + (bed_length % threads > 0)
    print(f"    splitting bed file into n = {math.ceil(bed_length / chunk_size)} chunks ...")
    # Split the BED file into chunks
    chunks = chunkify_bed(bed_file, chunk_size)

    # only use multiprocessing if more than 1 thread available/being used.
    if threads == 1:
        # loop through chromosomes
        print("multithreading skipped.")
        countData = []
        for idx,bed_chunk in enumerate(chunks):
            curr_chunk = f"chunk {idx+1}"
            counts_binned_chr = binned_read_counting(curr_chunk, bed_chunk, aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq)
            countData.append(counts_binned_chr)
        countData = [x for xs in countData for x in xs]

    else:
        print(f"multithreading using {threads} threads.")
        args_in = [[str(f"chunk {idx+1}"),bed_chunk,aln_files, nmode, blacklisting, bl_threshold, bases_filter, bases_threshold, readcount_mapq] for idx,bed_chunk in enumerate(chunks)]
        # print(args_in)
        with Pool(processes=threads) as pool:
            countData = [x for xs in list(pool.starmap(binned_read_counting, args_in)) for x in xs]

    filtered_counts, normalised_counts = filter_and_normalise(nmode=nmode, countData=countData)

    #----
    # 4. Get results and write out
    #----
    outfile = open(f"{outdir}/{sample}_raw_read_counts.tsv", "w")
    if blacklisting == True:
        header=['bin', 'chromosome','start','end','perc_known_bases', 'use_bin', 'tumour_read_count', 'normal_read_count']
    else: 
        header=['bin', 'chromosome','start','end','known_bases', 'overlap_blacklist', 'use_bin', 'tumour_read_count', 'normal_read_count']
    outfile.write('\t'.join(header)+'\n')
    for r in countData:
        Line = '\t'.join(r) + '\n'
        outfile.write(Line)
    outfile.close()

    
    
    # outfile2 = open(f"{outdir}/{sample}_read_counts_filtered.tsv", "w")
    # for r in filtered_counts:
    #     Line = '\t'.join(r) + '\n'
    #     outfile2.write(Line)
    # outfile2.close()

    log2_ratio_readcounts_path = f"{outdir}/{sample}_read_counts_{nmode}_log2r.tsv"
    outfile3 = open(log2_ratio_readcounts_path, "w")
    if blacklisting == True:
        header=['bin', 'chromosome','start','end','perc_known_bases', 'use_bin', 'log2r_copynumber']
    else: 
        header=['bin', 'chromosome','start','end','known_bases', 'overlap_blacklist', 'use_bin', 'log2r_copynumber']
    outfile3.write('\t'.join(header)+'\n')
    for r in normalised_counts:
        Line = '\t'.join(r) + '\n'
        outfile3.write(Line)
    outfile3.close()

    return log2_ratio_readcounts_path

if __name__ == "__main__":
    print("Read bin counter")