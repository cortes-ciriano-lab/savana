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
import sys

from time import time
from datetime import datetime
import argparse

__version__ = "1.2.4-dev"

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

# for developer debugging
def conditionally_decorate(dec, condition=False):
	""" whether to decorate a function (False by default)"""
	def decorator(func):
		print(f'{func}: {condition}')
		if not condition:
			# Return the function unchanged, not decorated.
			return func
		return dec(func)
	return decorator

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
		__ = int(string)
		return True
	except ValueError as _:
		return False
	except TypeError as _:
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

def get_chimeric_regions(read):
	""" get info for the chimeric regions of a read and store in a dict """
	if read.is_supplementary:
		return None
	chimeric_regions = []
	try:
		value = read.get_tag("SA")
	except KeyError as _:
		return chimeric_regions
	# parse the SA tag if it exists
	supp_alignments = [v.split(",") for v in value[:-1].split(";")]
	sa_keys = ['chrom', 'pos', 'strand', 'CIGAR', 'mapQ', 'NM']
	for supp_alignment in supp_alignments:
		chimeric_region = dict(zip(sa_keys, supp_alignment))
		chimeric_region['mapQ'] = int(chimeric_region['mapQ'])
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
	""" use the contigs file to return contigs and lengths - otherwise use index """
	if contig_file:
		with open(contig_file, encoding="utf-8") as f:
			contigs = f.readlines()
			contigs = [contig.rstrip() for contig in contigs]
			return contigs
	# otherwise, use the fai to get the contig names
	contigs = []
	with open(ref_index, encoding="utf-8") as f:
		tab_reader = csv.reader(f, delimiter='\t')
		for line in tab_reader:
			contig = line[0]
			contigs.append(contig)

	return contigs

def get_contig_lengths(ref_index):
	""" get the contig lengths from the reference """
	contig_lengths = {}
	with open(ref_index, encoding="utf-8") as f:
		tab_reader = csv.reader(f, delimiter='\t')
		for line in tab_reader:
			contig = line[0]
			length = line[1]
			contig_lengths[contig] = int(length)

	return contig_lengths

def generate_vcf_header(args, example_breakpoint):
	""" given a fasta file, index, and example breakpoint generate the VCF header """
	vcf_header_str = []
	vcf_header_str.extend([
		"##fileformat=VCFv4.2",
		f'##fileDate={datetime.now().strftime("%Y%m%d")}',
		f'##source=SAVANAv{__version__}'
	])
	# add contigs
	assembly_name = os.path.basename(args.ref)
	with open(args.ref_index) as f:
		reader = csv.reader(f, delimiter='\t')
		for line in list(reader):
			contig = line[0]
			length = line[1]
			vcf_header_str.append(f'##contig=<ID={contig},length={length},assembly={assembly_name}>')
	# generate command line args string
	cmd_string = '##savana_args="'
	for arg, value in vars(args).items():
		if value and arg != "func":
			cmd_string+=f' --{arg} {value}'
	cmd_string+='"'
	# add info fields
	# TODO: tumour only - there has to be a better way to do this
	vcf_header_str.extend([
		cmd_string,
		f'##reference={args.ref}',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
		'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
		'##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">'
		])
	if not args.tumour_only:
		vcf_header_str.extend([
			'##INFO=<ID=NORMAL_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting normal reads">'
			])
	vcf_header_str.extend([
		'##INFO=<ID=TUMOUR_READ_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumour reads">'
	])
	if not args.tumour_only:
		vcf_header_str.extend([
			'##INFO=<ID=NORMAL_ALN_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting normal alignments">'
		])
	vcf_header_str.extend([
		'##INFO=<ID=TUMOUR_ALN_SUPPORT,Number=1,Type=Integer,Description="Number of SV supporting tumour alignments">',
		'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
		'##INFO=<ID=TUMOUR_DP_BEFORE,Number=2,Type=Integer,Description="Local tumour depth in bin before the breakpoint(s) of an SV">',
		'##INFO=<ID=TUMOUR_DP_AT,Number=2,Type=Integer,Description="Local tumour depth in bin at the breakpoint(s) of an SV">',
		'##INFO=<ID=TUMOUR_DP_AFTER,Number=2,Type=Integer,Description="Local tumour depth in bin after the breakpoint(s) of an SV">'
	])
	if not args.tumour_only:
		vcf_header_str.extend([
			'##INFO=<ID=NORMAL_DP_BEFORE,Number=2,Type=Integer,Description="Local normal depth in bin before the breakpoint(s) of an SV">',
			'##INFO=<ID=NORMAL_DP_AT,Number=2,Type=Integer,Description="Local normal depth in bin at the breakpoint(s) of an SV">',
			'##INFO=<ID=NORMAL_DP_AFTER,Number=2,Type=Integer,Description="Local normal depth in bin after the breakpoint(s) of an SV">'
		])
	vcf_header_str.extend([
		'##INFO=<ID=TUMOUR_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of tumour variant-supporting reads to tumour read depth (DP) at breakpoint">'
	])
	if not args.tumour_only:
		vcf_header_str.extend([
		'##INFO=<ID=NORMAL_AF,Number=2,Type=Float,Description="Allele-fraction (AF) of normal variant-supporting reads to normal read depth (DP) at breakpoint">'
		])
	vcf_header_str.extend([
		'##INFO=<ID=BP_NOTATION,Number=1,Type=String,Description="+- notation format of variant (same for paired breakpoints)">',
		'##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of evidence for a breakpoint - CIGAR (INS, DEL, SOFTCLIP), SUPPLEMENTARY or mixture">',
		'##INFO=<ID=CLUSTERED_READS_TUMOUR,Number=1,Type=Integer,Description="Total number of tumour reads clustered at this location of any SV type">'
		])
	if not args.tumour_only:
		vcf_header_str.extend([
		'##INFO=<ID=CLUSTERED_READS_NORMAL,Number=1,Type=Integer,Description="Total number of normal reads clustered at this location of any SV type">'
		])
	vcf_header_str.extend([
		'##INFO=<ID=TUMOUR_ALT_HP,Number=3,Type=Integer,Description="Counts of SV-supporting reads belonging to each haplotype in the tumour sample (1/2/NA)">',
		'##INFO=<ID=TUMOUR_PS,Number=.,Type=String,Description="List of unique phase sets from the tumour supporting reads">'
	])
	if not args.tumour_only:
		vcf_header_str.extend([
			'##INFO=<ID=NORMAL_ALT_HP,Number=3,Type=Integer,Description="Counts of reads belonging to each haplotype in the normal sample (1/2/NA)">',
			'##INFO=<ID=NORMAL_PS,Number=.,Type=String,Description="List of unique phase sets from the normal supporting reads">'
		])
	vcf_header_str.extend([
		'##INFO=<ID=TUMOUR_TOTAL_HP_AT,Number=3,Type=Integer,Description="Counts of all reads at SV location belonging to each haplotype in the tumour sample (1/2/NA)">'
	])
	if not args.tumour_only:
		vcf_header_str.extend([
		'##INFO=<ID=NORMAL_TOTAL_HP_AT,Number=3,Type=Integer,Description="Counts of all reads at SV location belonging to each haplotype in the normal sample (1/2/NA)">'
	])
	# add the stat info fields
	breakpoint_stats_origin = example_breakpoint.originating_cluster.get_stats().keys()
	breakpoint_stats_end = example_breakpoint.end_cluster.get_stats().keys()
	for stat in breakpoint_stats_origin:
		vcf_header_str.append(f'##INFO=<ID=ORIGIN_{stat.upper()},Number=1,Type=Float,Description="Originating cluster value for {stat}">')
	for stat in breakpoint_stats_end:
		vcf_header_str.append(f'##INFO=<ID=END_{stat.upper()},Number=1,Type=Float,Description="End cluster value for {stat}">')
	# add the final header line
	sample_name = os.path.splitext(os.path.basename(args.tumour))[0]
	vcf_header_str.append("#"+"\t".join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name]))

	return "\n".join(vcf_header_str)+"\n"

def time_function(desc, checkpoints, time_str, final=False):
	""" prints the number of seconds elapsed compared to previous checkpoint """
	checkpoints.append(time())
	if not final:
		formatted_time = f'{desc:<60}{round(checkpoints[-1] - checkpoints[-2], 3)} seconds'
	else:
		formatted_time = f'{desc:<60}{round(checkpoints[-1] - checkpoints[0], 3)} seconds\n'
	time_str.append(formatted_time)
	print(formatted_time)
	return

def check_outdir(args_outdir, args_overwrite, illegal=None):
	# create output dir if it doesn't exist
	outdir = os.path.join(os.getcwd(), args_outdir)
	if not os.path.exists(outdir):
		print(f'Creating directory {outdir} to store results')
		os.mkdir(outdir)
	if args_overwrite:
		# don't check for files, overwrite them
		return outdir
	if not illegal:
		# throw error if ANY files present
		if os.listdir(outdir):
			sys.exit(f'Output directory "{outdir}" already exists and contains files. Please remove the files or supply a different directory name.')
	else:
		# check if illegal files exist in outdir
		for f in os.listdir(outdir):
			if f.endswith(illegal):
				sys.exit(f'Output directory "{outdir}" already exists and contains {illegal} files which may be overwritten. Please remove the files or supply a different directory name.')

	return outdir

def check_tmpdir(args_tmpdir, outdir, args_overwrite, illegal=None):
	# create output dir if it doesn't exist
	print(outdir)
	print(args_overwrite)
	tmpdir = os.path.join(outdir, args_tmpdir)
	if not os.path.exists(tmpdir):
		print(f'Creating directory {tmpdir} to store temp files during hetSNP allele counting')
		os.mkdir(tmpdir)
	if args_overwrite:
		# don't check for files, overwrite them
		return tmpdir
	if not illegal:
		# throw error if ANY files present
		if os.listdir(tmpdir):
			sys.exit(f'Temp directory "{tmpdir}" already exists and contains files. Please remove the files or supply a different directory name.')
	else:
		# check if illegal files exist in outdir
		for f in os.listdir(tmpdir):
			if f.endswith(illegal):
				sys.exit(f'Temp directory "{tmpdir}" already exists and contains {illegal} files which may be overwritten. Please remove the files or supply a different directory name.')

	return tmpdir

def clean_tmpdir(args_tmpdir, outdir):
	tmpdir = os.path.join(outdir, args_tmpdir)
	# remove tmpdir if empty
	if not os.listdir(tmpdir):
		os.rmdir(tmpdir)

# credit to stackoverflow: https://stackoverflow.com/questions/55324449/how-to-specify-a-minimum-or-maximum-float-value-with-argparse
def float_range(mini,maxi):
	"""Return function handle of an argument type function for ArgumentParser checking a float range: mini <= arg <= maxi
		mini - minimum acceptable argument
		maxi - maximum acceptable argument
	"""

	# Define the function with default arguments
	def float_range_checker(arg):
		"""New Type function for argparse - a float within predefined range."""

		try:
			f = float(arg)
		except ValueError:
			raise argparse.ArgumentTypeError("must be a floating point number")
		if f < mini or f > maxi:
			raise argparse.ArgumentTypeError("must be in range [" + str(mini) + " .. " + str(maxi)+"]")
		return f

	# Return function handle to checking function
	return float_range_checker

if __name__ == "__main__":
	print("Helper functions for SAVANA")
