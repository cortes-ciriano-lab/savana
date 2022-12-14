"""
SAVANA strucural variant caller for long-read data - main program
Created: 21/09/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import sys
import os
import argparse

from time import time
from math import ceil
from multiprocessing import Pool, cpu_count

import pysam

import savana.helper as helper
import savana.validation as validation
from savana.breakpoints import *
from savana.clusters import *

logo = """
███████  █████  ██    ██  █████  ███    ██  █████
██      ██   ██ ██    ██ ██   ██ ████   ██ ██   ██
███████ ███████ ██    ██ ███████ ██ ██  ██ ███████
     ██ ██   ██  ██  ██  ██   ██ ██  ██ ██ ██   ██
███████ ██   ██   ████   ██   ██ ██   ████ ██   ██
"""
__version__ = "0.2.0"

def pool_get_potential_breakpoints(bam_files, args):
	""" split the genome into 500kBp chunks and identify PotentialBreakpoints """
	chunk_size = 500000 # 5.0e5 or (half a million)
	pool_potential = Pool(processes=args.threads)
	pool_potential_args = []
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	for label, bam_file in bam_files.items():
		for contig in bam_file.get_index_statistics():
			if contig.contig not in contigs_to_consider:
				if args.debug:
					print(f'Skipping reads aligned to {contig.contig} - not in contigs file')
				continue
			if contig.mapped == 0:
				continue
			chrom_length = int(bam_file.get_reference_length(contig.contig))
			if chrom_length > chunk_size:
				# split the chrom into parts
				num_intervals = ceil(chrom_length/chunk_size) + 1
				start_pos = 0
				for i in range(1, num_intervals):
					end_pos = start_pos + chunk_size
					end_pos = chrom_length if end_pos > chrom_length else end_pos # don't extend past end
					pool_potential_args.append((bam_file.filename, args, label, contigs_to_consider, contig.contig, start_pos, end_pos))
					start_pos = end_pos + 1
			else:
				pool_potential_args.append((bam_file.filename, args, label, contigs_to_consider, contig.contig))
	potential_breakpoints_results = pool_potential.starmap(get_potential_breakpoints, pool_potential_args)
	pool_potential.close()
	pool_potential.join()
	return potential_breakpoints_results

def pool_cluster_breakpoints(args, chrom_potential_breakpoints):
	""" perform initial clustering of Potential Breakpoints """
	pool_clustering = Pool(processes=args.threads)
	pool_clustering_args = []
	for breakpoints in chrom_potential_breakpoints.values():
		pool_clustering_args.append((breakpoints, args))
	clustering_results = pool_clustering.starmap(cluster_breakpoints, pool_clustering_args)
	pool_clustering.close()
	pool_clustering.join()
	clusters = {} # collect results
	for result in clustering_results:
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			clusters.setdefault(bp_type, []).extend(result[bp_type])
	return clusters

def pool_output_clusters(args, clusters, outdir):
	""" output trimmed fastqs of the reads in each cluster """
	pool_output = Pool(processes=args.threads)
	pool_output_args = []
	# split list into equal chunks from https://stackoverflow.com/a/2135920
	quotient, remainder = divmod(len(clusters), args.threads)
	clusters_split = (clusters[i*quotient+min(i, remainder):(i+1)*quotient+min(i+1, remainder)] for i in range(args.threads))
	for split in clusters_split:
		pool_output_args.append((split, outdir))
	pool_output.starmap(output_clusters, pool_output_args)
	pool_output.close()
	pool_output.join()

def apply_somatic_filters(validated_breakpoints):
	""" use heuristics to eliminate noise and germline variants """
	somatic_breakpoints_lenient = []
	for bp in validated_breakpoints:
		if bp.support['tumour'] == 0:
			continue
		elif bp.support['normal']/bp.support['tumour'] < 0.1:
			originating_cluster_stats = bp.originating_cluster.get_stats()
			end_cluster_stats = bp.end_cluster.get_stats()
			if originating_cluster_stats['starts_std_dev'] < 150 and originating_cluster_stats['event_heuristic'] < 3:
				if bp.breakpoint_notation == "<INS>" and bp.support['tumour'] > 25:
					somatic_breakpoints_lenient.append(bp)
				elif bp.breakpoint_notation != "<INS>" and bp.support['tumour'] > 5:
					somatic_breakpoints_lenient.append(bp)

	somatic_breakpoints_strict = []
	for bp in validated_breakpoints:
		if bp.support['normal'] > 0:
			continue
		if bp.support['tumour'] > 7:
			originating_cluster_stats = bp.originating_cluster.get_stats()
			end_cluster_stats = bp.end_cluster.get_stats()
			if bp.support['tumour'] > 12 and originating_cluster_stats['uncertainty'] <= 15:
				if originating_cluster_stats['event_heuristic'] <= 0.025:
					somatic_breakpoints_strict.append(bp)
					continue
			elif bp.support['tumour'] > 12 and end_cluster_stats['uncertainty'] <= 30:
				somatic_breakpoints_strict.append(bp)
				continue
			elif end_cluster_stats['uncertainty'] <= 10:
				somatic_breakpoints_strict.append(bp)
				continue

	return somatic_breakpoints_lenient, somatic_breakpoints_strict

def spawn_processes(args, bam_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using multiprocessing with {args.threads} threads\n')
	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints_results = pool_get_potential_breakpoints(bam_files, args)
	if args.debug:
		time_function("Identified potential breakpoints", checkpoints, time_str)
	# collect results per chrom
	chrom_potential_breakpoints = {}
	for result in potential_breakpoints_results:
		for chrom, breakpoints in result.items():
			chrom_potential_breakpoints.setdefault(chrom,[]).extend(breakpoints)

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clusters = pool_cluster_breakpoints(args, chrom_potential_breakpoints)
	if args.debug:
		time_function("Clustered potential breakpoints", checkpoints, time_str)

	if args.debug:
		# 2.1) OUTPUT CLUSTERS
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			pool_output_clusters(args, clusters[bp_type], outdir)
		time_function("Output originating clusters", checkpoints, time_str)

	# 3) CALL BREAKPOINTS
	validated_breakpoints = call_breakpoints(clusters, args.buffer)
	if args.debug:
		time_function("Called consensus breakpoints", checkpoints, time_str)

	# 3.1) OUTPUT BREAKPOINTS
	bedpe_string = ''
	vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, validated_breakpoints[0])
	read_support_string = ''
	variant_stats_string, variant_stats_cols = helper.generate_variant_stats_header(validated_breakpoints[0])
	ref_fasta = pysam.FastaFile(args.ref)
	for count, bp in enumerate(validated_breakpoints):
		bedpe_string += bp.as_bedpe(count)
		vcf_string += bp.as_vcf(ref_fasta)
		read_support_string += bp.as_read_support(count)
		variant_stats_string += bp.as_variant_stats(count, variant_stats_cols)
	with open(os.path.join(outdir, "sv_breakpoints.vcf"), 'w') as output:
		output.write(vcf_string)
	with open(os.path.join(outdir, "sv_breakpoints.bedpe"), 'w') as output:
		output.write(bedpe_string)
	with open(os.path.join(outdir, "sv_breakpoints_read_support.tsv"), 'w') as output:
		output.write(read_support_string)
	with open(os.path.join(outdir, "variant.stats"), 'w') as output:
		output.write(variant_stats_string)
	if args.debug:
		time_function("Output consensus breakpoints", checkpoints, time_str)

	# 4) CATEGORIZE BREAKPOINTS
	somatic_breakpoints_lenient, somatic_breakpoints_strict = apply_somatic_filters(validated_breakpoints)
	# output lenient vcf
	lenient_vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, validated_breakpoints[0])
	for bp in somatic_breakpoints_lenient:
		lenient_vcf_string += bp.as_vcf(ref_fasta)
	with open(os.path.join(outdir, "somatic.sv_breakpoints.lenient.vcf"), 'w') as output:
		output.write(lenient_vcf_string)
	# output strict vcf
	strict_vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, validated_breakpoints[0])
	for bp in somatic_breakpoints_strict:
		strict_vcf_string += bp.as_vcf(ref_fasta)
	with open(os.path.join(outdir, "somatic.sv_breakpoints.strict.vcf"), 'w') as output:
		output.write(strict_vcf_string)

	if args.debug:
		time_function("Applied somatic filters", checkpoints, time_str)

	return clusters, checkpoints, time_str

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

def main():
	""" main function for SAVANA - collects command line arguments and executes algorithm """
	parser = argparse.ArgumentParser(description="SAVANA - somatic SV caller")
	parser.add_argument('--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
	parser.add_argument('--normal', nargs='?', type=str, required=True, help='Normal BAM file (must have index)')
	parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
	parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
	parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")
	parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
	parser.add_argument('--mapq', nargs='?', type=int, default=5, help='MAPQ filter on reads which are considered (default=5)')
	parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer to add when clustering adjacent potential breakpoints (default=10)')
	parser.add_argument('--depth', nargs='?', type=int, default=3, help='Threshold for number of supporting reads (default=3)')
	parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
	parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
	parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
	parser.add_argument('--validation', nargs='?', type=str, required=False, help='VCF file to use as validation (optional)')
	args = parser.parse_args()

	print(logo)
	print(f'Version {__version__} - beta')
	src_location = __file__
	print(f'Source: {src_location}\n')

	# create output dir if it doesn't exist
	outdir = os.path.join(os.getcwd(), args.outdir)
	if not os.path.exists(outdir):
		print(f'Creating directory {outdir} to store results')
		os.mkdir(outdir)
	elif os.listdir(outdir):
		sys.exit(f'Output directory "{outdir}" already exists and contains files. Please remove the files or supply a different directory name.')

	# set number of threads to cpu count if none set
	if not args.threads:
		args.threads = cpu_count()

	# read bam files (must have bai)
	bam_files = {
		'tumour': pysam.AlignmentFile(args.tumour, "rb"),
		'normal': pysam.AlignmentFile(args.normal, "rb")
	}

	# confirm ref and ref fasta index exist
	if not os.path.exists(args.ref):
		sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
	elif args.ref_index and not os.path.exists(args.ref_index):
		sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
	elif not os.path.exists(f'{args.ref}.fai'):
		sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
	else:
		args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
		print(f'Using {args.ref_index} as reference fasta index')

	# initialize timing
	checkpoints = [time()]
	time_str = []

	# run SAVANA processes
	consensus_clusters, checkpoints, time_str = spawn_processes(args, bam_files, checkpoints, time_str, outdir)

	if args.debug:
		write_cluster_bed(consensus_clusters, outdir)
		calculate_cluster_stats(consensus_clusters, outdir)

	# validate strict
	output_vcf = os.path.join(outdir, 'somatic.sv_breakpoints.strict.vcf')
	if args.validation:
		try:
			validation.validate_vcf(outdir, output_vcf, args.validation, 'strict')
		except Exception as e:
			print(f'\nWARNING: Validation of breakpoints against {args.validation} failed due to "{str(e)}"')
			print(f'You can retry by running "python savana/validation.py --outdir testing --input {output_vcf} --validation {args.validation}"')

	# validate lenient
	output_vcf = os.path.join(outdir, 'somatic.sv_breakpoints.lenient.vcf')
	if args.validation:
		try:
			validation.validate_vcf(outdir, output_vcf, args.validation, 'lenient')
		except Exception as e:
			print(f'\nWARNING: Validation of breakpoints against {args.validation} failed due to "{str(e)}"')
			print(f'You can retry by running "python savana/validation.py --outdir testing --input {output_vcf} --validation {args.validation}"')

	time_function("Total time", checkpoints, time_str, final=True)
	f = open(os.path.join(outdir, 'time.log'), "w+")
	f.write("\n".join(time_str))
	f.write("\n")
	f.close()

if __name__ == "__main__":
	main()
