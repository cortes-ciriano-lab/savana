"""
Module containing misc. useful functions for SAVANA
Created: 26/01/2022
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import re
import os
import csv

from time import time
from datetime import datetime

__version__ = "0.3.0"

samflag_desc_to_number = {
	"BAM_CMATCH": 0, # M
	"BAM_CINS": 1, # I
	"BAM_CDEL": 2, # D
	"BAM_CREF_SKIP": 3, # N
	"BAM_CSOFT_CLIP": 4, # S
	"BAM_CHARD_CLIP": 5, # H
	"BAM_CPAD": 6, # P
	"BAM_CEQUAL": 7, # =
	"BAM_CDIFF": 8, # X
	"BAM_CBACK": 9 # B
}
samflag_number_to_desc = {n: d for d, n in samflag_desc_to_number.items()}

samflag_number_to_letter = {
	0: "M",
	1: "I",
	2: "D",
	3: "N",
	4: "S",
	5: "H",
	6: "P",
	7: "=",
	8: "X",
	9: "B"
}
samflag_letter_to_number = {l: n for n, l in samflag_number_to_letter.items()}

"""
                                                             Consumes	Consumes
Op  BAM Description                                             query	reference
M   0   alignment match (can be a sequence match or mismatch)   yes		yes
I   1   insertion to the reference                              yes		no
D   2   deletion from the reference                             no		yes
N   3   skipped region from the reference                       no		yes
S   4   soft clipping (clipped sequences present in SEQ)        yes		no
H   5   hard clipping (clipped sequences NOT present in SEQ)    no		no
P   6   padding (silent deletion from padded reference)         no		no
=   7   sequence match                                          yes		yes
X   8   sequence mismatch                                       yes		yes
"""
consumes_query = {
	0: True,
	1: True,
	2: False,
	3: False,
	4: True,
	5: False,
	6: False,
	7: True,
	8: True
}
consumes_reference = {
	0: True,
	1: False,
	2: True,
	3: True,
	4: False,
	5: False,
	6: False,
	7: True,
	8: True
}

def reverse_complement(sequence):
	""" dna reverse complement of bases """
	bases = {'A':'T', 'T':'A', 'C': 'G', 'G': 'C', 'N': 'N'}
	rev_comp = [bases[b] for b in sequence[::-1]]
	return ''.join(rev_comp)

breakend_type = {
	'N[', # extends right after N
	'[N', # extends right before N
	'N]', # extends left after N
	']N' # extends left before N
}

def is_int(string):
	""" return true if string can be converted to int, False otherwise """
	try:
		_ = int(string)
		return True
	except ValueError as e:
		return False

def flatten(list_of_lists):
	""" given a list of lists, flatten to a 1d list """
	return [item for sublist in list_of_lists for item in sublist]

def sum_cigar(cigarstring):
	""" add up ALL lengths in a CIGAR string (not just ones that consume the query) """
	cigar_split = re.split('[MIDNSHP=X]', cigarstring)[:-1] # split on characters
	cigar_split = [int(x) for x in cigar_split]
	return sum(cigar_split)

def sum_consumed_query(cigarstring):
	""" add up the lengths in a CIGAR string that consume the query """
	cigar_split = re.split('([MIDNSHP=X])', cigarstring)[:-1]
	sum_consumed = 0
	for i in range(0,len(cigar_split)-1,2):
		if consumes_query[samflag_letter_to_number[cigar_split[i+1]]]:
			sum_consumed += int(cigar_split[i])

	return sum_consumed

def sum_consumed_reference(cigarstring):
	""" add up the lengths in a CIGAR string that consume the query """
	cigar_split = re.split('([MIDNSHP=X])', cigarstring)[:-1]
	sum_consumed = 0
	for i in range(0,len(cigar_split)-1,2):
		if consumes_reference[samflag_letter_to_number[cigar_split[i+1]]]:
			sum_consumed += int(cigar_split[i])

	return sum_consumed

def get_cigartuples(cigarstring):
	""" given a cigarstring, return a list of tuples [(operation, length)] """
	cigar_tuples = []
	cigar_split = re.split('([MIDNSHP=X])', cigarstring)[:-1]
	for i in range(0,len(cigar_split)-1,2):
		cigar_tuples.append((samflag_letter_to_number[cigar_split[i+1]], int(cigar_split[i])))

	return cigar_tuples

def trim_supplementary(cigarstring):
	""" trim supplementary edges of a CIGAR string if there are any """
	# e.g.) 10S200M33D42M15S -> 200M33D42M
	trimmed_cigar = ''
	cigar_split = re.split('([MIDNSHP=X])', cigarstring)[:-1]
	for i in range(0, len(cigar_split), 2):
		if cigar_split[i+1] != 'S':
			trimmed_cigar += ''.join(cigar_split[i:i+2])

	return trimmed_cigar

def get_read_boundaries(read):
	""" given a read get its start and end of both aligned and softclipped region.
	|~~~~~|--------|~~~~~|
	rs    qs      qe     re
	ranges:
		rs:(rs + |softclipped bases|)
		qs:qe
		(re - |softclipped bases|):re
	For now, just using qs and qe
	"""
	boundaries = {'chr': read.reference_name, 'range': [int(read.reference_start), int(read.reference_end)]}
	return boundaries

def get_clipping(cigar_tuples, is_reverse):
	""" return direction-adjusted clipping of an alignment from the cigar tuple """

	clipping = {'left_softclip': 0, 'right_softclip': 0}
	if cigar_tuples[0][0] == samflag_desc_to_number["BAM_CSOFT_CLIP"]:
		if is_reverse:
			clipping['right_softclip'] = cigar_tuples[0][1]
		else:
			clipping['left_softclip'] = cigar_tuples[0][1]
	if cigar_tuples[-1][0] == samflag_desc_to_number["BAM_CSOFT_CLIP"]:
		if is_reverse:
			clipping['left_softclip'] = cigar_tuples[-1][1]
		else:
			clipping['right_softclip'] = cigar_tuples[-1][1]

	return clipping

def get_chimeric_regions(read, mapq_filter):
	""" get info for the chimeric regions of a read and store in a dict """
	chimeric_regions = []
	sa_keys = ['chrom', 'pos', 'strand', 'CIGAR', 'mapQ', 'NM']
	for tag, value in read.get_tags():
		if tag == "SA":
			# (chrom, pos, strand, CIGAR, mapQ, NM)
			supp_alignments = [v.split(",") for v in value[:-1].split(";")]
			for supp_alignment in supp_alignments:
				chimeric_region = dict(zip(sa_keys, supp_alignment))
				if int(chimeric_region['mapQ']) < mapq_filter:
					# only consider those above quality threshold
					continue
				cigar_split = re.split('([MIDNSHP=X])', chimeric_region['CIGAR'])[:-1]
				if chimeric_region['strand'] == "-":
					# reverse CIGAR if the direction of the supplementary doesn't match the primary
					cigar_split.reverse()
					cigar_split = [cigar_split[f(x)] for x in range(0,len(cigar_split),2) for f in (lambda x: x+1, lambda x: x)]
				# calculate the number of bases from the SA CIGAR string that are consuming the query and reference
				alignment_sum = sum_consumed_query(trim_supplementary(chimeric_region['CIGAR']))
				reference_sum = sum_consumed_reference(trim_supplementary(chimeric_region['CIGAR']))
				left_softclip = int(cigar_split[0]) if cigar_split[1] == 'S' else 0
				right_softclip = int(cigar_split[-2]) if cigar_split[-1] == 'S' else 0
				# track values in chimeric region dicts
				chimeric_region['left_softclip'] = left_softclip
				chimeric_region['right_softclip'] = right_softclip
				chimeric_region['consumed_query'] = alignment_sum
				chimeric_region['consumed_reference'] = reference_sum
				chimeric_region['seen'] = False
				chimeric_regions.append(chimeric_region)
	return chimeric_regions

def get_contigs(contig_file, ref_index):
	""" given a file of contigs to consider, return them in a list """
	if contig_file:
		with open(contig_file, encoding="utf-8") as f:
			contigs = f.readlines()
			contigs = [contig.rstrip() for contig in contigs]
			return contigs
	elif ref_index:
		# use the fai to get the contig names
		with open(ref_index, encoding="utf-8") as f:
			tab_reader = csv.reader(f, delimiter='\t')
			contigs = []
			for line in tab_reader:
				contig = line[0]
				contigs.append(contig)
			return contigs
	return None

def generate_vcf_header(ref_fasta, ref_fasta_index, tumour_file, example_breakpoint):
	""" given a fasta file and index, generate the VCF header """
	vcf_header_str = []
	vcf_header_str.extend([
		"##fileformat=VCFv4.2",
		f'##fileDate={datetime.now().strftime("%Y%m%d")}',
		"##source=SAVANA.Beta",
		f'##reference={ref_fasta}',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
		'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
		'##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">',
		'##INFO=<ID=NORMAL_SUPPORT,Number=1,Type=Float,Description="Number of variant supporting normal reads">',
		'##INFO=<ID=TUMOUR_SUPPORT,Number=1,Type=Float,Description="Number of variant supporting tumour reads">',
		'##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
		'##INFO=<ID=ORIGINATING_CLUSTER,Number=.,Type=String,Description="SAVANA internal originating cluster id supporting variant">',
		'##INFO=<ID=END_CLUSTER,Number=.,Type=String,Description="SAVANA internal end cluster id supporting variant">',
		'##INFO=<ID=BP_NOTATION,Number=1,Type=String,Description="+- notation format of variant (same for paired breakpoints)">'
	])
	breakpoint_stats_origin = example_breakpoint.originating_cluster.get_stats().keys()
	breakpoint_stats_end = example_breakpoint.end_cluster.get_stats().keys()
	for stat in breakpoint_stats_origin:
		vcf_header_str.append(f'##INFO=<ID=ORIGIN_{stat.upper()},Number=1,Type=Float,Description="Originating cluster value for {stat}">')
	for stat in breakpoint_stats_end:
		vcf_header_str.append(f'##INFO=<ID=END_{stat.upper()},Number=1,Type=Float,Description="End cluster value for {stat}">')
	assembly_name = os.path.basename(ref_fasta)
	with open(ref_fasta_index) as f:
		reader = csv.reader(f, delimiter='\t')
		for line in list(reader):
			contig = line[0]
			length = line[1]
			vcf_header_str.append(f'##contig=<ID={contig},length={length},assembly={assembly_name}>')
	sample_name = os.path.splitext(os.path.basename(tumour_file))[0]
	vcf_header_str.append("#"+"\t".join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name]))

	return "\n".join(vcf_header_str)+"\n"

def time_function(desc, checkpoints, time_str, final=False):
	""" prints the number of seconds elapsed compared to previous checkpoint """
	checkpoints.append(time())
	if not final:
		formatted_time = f'{desc:<40}{round(checkpoints[-1] - checkpoints[-2], 2)} seconds'
	else:
		formatted_time = f'{desc:<40}{round(checkpoints[-1] - checkpoints[0], 2)} seconds'
	time_str.append(formatted_time)
	print(formatted_time)
	return

if __name__ == "__main__":
	print("Helper functions for SAVANA")
