"""
SAVANA strucural variant caller for long-read data - run
Created: 17/02/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import os

from math import ceil, floor
from multiprocessing import Pool

import numpy as np
import pysam
import pysam.bcftools as bcftools
import pybedtools

import savana.helper as helper
from savana.breakpoints import get_potential_breakpoints, call_breakpoints, compute_depth
from savana.clusters import cluster_breakpoints, output_clusters

"""
# developer dependencies
from memory_profiler import profile
from pympler import muppy, summary, refbrowser
import objgraph
"""

def multithreading_get_potential_breakpoints(aln_files, args):
	""" split the genome into chunks and identify PotentialBreakpoints via multithreading """
	from concurrent.futures import ThreadPoolExecutor
	from threading import Lock

	executor = ThreadPoolExecutor(max_workers=args.threads)
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	get_potential_breakpoints_args = {
		'aln_filename': [],
		'is_cram': [],
		'ref': [],
		'length': [],
		'mapq': [],
		'label': [],
		'contigs_to_consider': [],
		'contig': [],
		'start_pos': [],
		'end_pos': []
	}
	for label, aln_file in aln_files.items():
		for contig, contig_length in contig_lengths.items():
			if contig not in contigs_to_consider:
				continue
			if contig_length > args.chunksize:
				# split the chrom into parts
				num_intervals = floor(contig_length/args.chunksize) + 1
				start_pos = 0
				for i in range(1, num_intervals):
					end_pos = start_pos + args.chunksize
					end_pos = contig_length if end_pos > contig_length else end_pos # don't extend past end
					get_potential_breakpoints_args['aln_filename'].append(aln_file.filename)
					get_potential_breakpoints_args['is_cram'].append(args.is_cram)
					get_potential_breakpoints_args['ref'].append(args.ref)
					get_potential_breakpoints_args['length'].append(args.length)
					get_potential_breakpoints_args['mapq'].append(args.mapq)
					get_potential_breakpoints_args['label'].append(label)
					get_potential_breakpoints_args['contigs_to_consider'].append(contigs_to_consider)
					get_potential_breakpoints_args['contig'].append(contig)
					get_potential_breakpoints_args['start_pos'].append(start_pos)
					get_potential_breakpoints_args['end_pos'].append(end_pos)
					start_pos = end_pos + 1
			else:
				get_potential_breakpoints_args['aln_filename'].append(aln_file.filename)
				get_potential_breakpoints_args['is_cram'].append(args.is_cram)
				get_potential_breakpoints_args['ref'].append(args.ref)
				get_potential_breakpoints_args['length'].append(args.length)
				get_potential_breakpoints_args['mapq'].append(args.mapq)
				get_potential_breakpoints_args['label'].append(label)
				get_potential_breakpoints_args['contigs_to_consider'].append(contigs_to_consider)
				get_potential_breakpoints_args['contig'].append(contig)
				get_potential_breakpoints_args['start_pos'].append(0)
				get_potential_breakpoints_args['end_pos'].append(contig_length)
	results = executor.map(
		get_potential_breakpoints,
		get_potential_breakpoints_args['aln_filename'],
		get_potential_breakpoints_args['is_cram'],
		get_potential_breakpoints_args['ref'],
		get_potential_breakpoints_args['length'],
		get_potential_breakpoints_args['mapq'],
		get_potential_breakpoints_args['label'],
		get_potential_breakpoints_args['contigs_to_consider'],
		get_potential_breakpoints_args['contig'],
		get_potential_breakpoints_args['start_pos'],
		get_potential_breakpoints_args['end_pos']
	)

	return results


def pool_get_potential_breakpoints(aln_files, args):
	""" split the genome into chunks and identify PotentialBreakpoints """
	#pool_potential = Pool(processes=args.threads, maxtasksperchild=args.maxtasksperchild) if args.maxtasksperchild else Pool(processes=args.threads)
	pool_potential_args = []
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	if args.debug:
		print(f' > Setting chunksize for split to {args.chunksize}')
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	for label, aln_file in aln_files.items():
		for contig, contig_length in contig_lengths.items():
			if contig not in contigs_to_consider:
				continue
			if contig_length > args.chunksize:
				# split the chrom into parts
				num_intervals = floor(contig_length/args.chunksize) + 1
				start_pos = 0
				for i in range(1, num_intervals):
					end_pos = start_pos + args.chunksize
					end_pos = contig_length if end_pos > contig_length else end_pos # don't extend past end
					pool_potential_args.append((aln_file.filename, args.is_cram, args.ref, args.length, args.mapq, label, contigs_to_consider, contig, start_pos, end_pos))
					start_pos = end_pos + 1
			else:
				pool_potential_args.append((aln_file.filename, args.is_cram, args.ref, args.length, args.mapq, label, contigs_to_consider, contig, 0, contig_length))
	if args.debug:
		print(f' > Submitting {len(pool_potential_args)} "get_potential_breakpoints" tasks to {args.threads} worker threads ({args.maxtasksperchild} maxtasksperchild)')

	potential_breakpoints_results = pool_potential.starmap(get_potential_breakpoints, pool_potential_args)
	pool_potential.close()
	pool_potential.join()
	return potential_breakpoints_results

def pool_cluster_breakpoints(threads, buffer, ins_buffer, chrom_potential_breakpoints):
	""" perform initial clustering of Potential Breakpoints """
	pool_clustering = Pool(processes=threads)
	pool_clustering_args = []
	for chrom, breakpoints in chrom_potential_breakpoints.items():
		pool_clustering_args.append((chrom, breakpoints, buffer, ins_buffer))
	clustering_results = pool_clustering.starmap(cluster_breakpoints, pool_clustering_args)
	pool_clustering.close()
	pool_clustering.join()
	clusters = {}
	""" collect & store results in multi-level dict
	clusters = {
		'chr1': {
			{
				'+-': [cluster, cluster, ...],
				'--': [cluster, cluster, ...]
			}
		}
		...
	} """
	for chrom, result in clustering_results:
		if chrom not in clusters:
			clusters[chrom] = result

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

def pool_call_breakpoints(threads, buffer, length, depth, clusters, debug):
	""" parallelise the identification of consensus breakpoints """
	pool_calling = Pool(processes=threads)
	pool_calling_args = []

	for chrom, chrom_clusters in clusters.items():
		pool_calling_args.append((chrom_clusters, buffer, length, depth, chrom))
	calling_results = pool_calling.starmap(call_breakpoints, pool_calling_args)

	breakpoint_dict_chrom = {}
	seen_cluster_uids = {}
	pruned_clusters = {} if debug else None
	for result_breakpoints, result_pruned_clusters, result_chrom in calling_results:
		# collect breakpoint calling results
		breakpoint_dict_chrom[result_chrom] = result_breakpoints
		if debug:
			for bp_type in result_pruned_clusters.keys():
				for cluster in result_pruned_clusters[bp_type]:
					if cluster.uid not in seen_cluster_uids:
						pruned_clusters.setdefault(bp_type, []).append(cluster)
						seen_cluster_uids[cluster.uid] = True

	return breakpoint_dict_chrom, pruned_clusters

def multithreading_compute_depth(threads, breakpoint_dict_chrom, contig_coverages_merged, debug):
	""" computes the depth of breakpoints using coverage arrays """
	from concurrent.futures import ThreadPoolExecutor
	from threading import Lock

	executor = ThreadPoolExecutor(max_workers=threads+6)
	lock = Lock()
	results = executor.map(
		compute_depth,
		[breakpoints for breakpoints in breakpoint_dict_chrom.values()],
		[contig_coverages_merged]*len(breakpoint_dict_chrom.keys()),
		[lock]*len(breakpoint_dict_chrom.keys())
	)
	result_collector = []
	for i, result in enumerate(results):
		result_collector.append((i, len(result)))
	if debug:
		print(f' > Finished computing depth, work split across {len(result_collector)} units in total')

	return breakpoint_dict_chrom

def spawn_processes(args, aln_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using {args.threads} thread(s)\n')

	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints_results = pool_get_potential_breakpoints(aln_files, args)
	helper.time_function("Identified potential breakpoints", checkpoints, time_str)

	# collect results per chrom
	chrom_potential_breakpoints = {}
	contig_coverages_merged = {}
	for result in potential_breakpoints_results:
		#potential_breakpoints_dict = result
		potential_breakpoints_dict = result[0]
		contig_coverages = result[1]
		result_chrom = None
		for chrom, potential_breakpoints in potential_breakpoints_dict.items():
			result_chrom = chrom if not result_chrom else result_chrom
			chrom_potential_breakpoints.setdefault(chrom,[]).extend(potential_breakpoints)
		contig_coverages_merged.setdefault(contig_coverages.pop('contig'), []).append(contig_coverages)
	# get rid (heavy memory footprint - no longer needed)
	del potential_breakpoints_results

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clusters = pool_cluster_breakpoints(args.threads, args.buffer, args.insertion_buffer, chrom_potential_breakpoints)
	helper.time_function("Clustered potential breakpoints", checkpoints, time_str)

	# 3) CALL BREAKPOINTS FROM CLUSTERS
	breakpoint_dict_chrom, pruned_clusters = pool_call_breakpoints(args.threads, args.buffer, args.length, args.depth, clusters, args.debug)
	helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	if args.debug:
		total_breakpoints = 0
		for _, b in breakpoint_dict_chrom.items():
			total_breakpoints+=len(b)
		print(f' > Total number of breakpoints post-calling: {total_breakpoints}')

	#TODO: UNCOMMENT LATER - SKIP THIS DEBUGGING STEP FOR NOW
	"""
	if args.debug:
		# 3.1) OUTPUT CLUSTERS
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			if bp_type in pruned_clusters:
				pool_output_clusters(args, pruned_clusters[bp_type], outdir)
		helper.time_function("Output pruned clusters", checkpoints, time_str)
	"""

	# 4) COMPUTE LOCAL DEPTH
	multithreading_compute_depth(args.threads, breakpoint_dict_chrom, contig_coverages_merged, args.debug)
	helper.time_function("Computed local depth for breakpoints", checkpoints, time_str)

	# 5) OUTPUT BREAKPOINTS
	# define filenames
	vcf_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.vcf')
	bedpe_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.bedpe')
	tsv_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints_read_support.tsv')
	# build strings
	ref_fasta = pysam.FastaFile(args.ref)
	bedpe_string = ''
	vcf_string = helper.generate_vcf_header(args, breakpoint_dict_chrom[list(breakpoint_dict_chrom.keys())[0]][0])
	read_support_string = 'VARIANT_ID\tTUMOUR_SUPPORTING_READS\tNORMAL_SUPPORTING_READS\n'
	count = 0
	for chrom, chrom_breakpoints in breakpoint_dict_chrom.items():
		for bp in chrom_breakpoints:
			bedpe_string += bp.as_bedpe(count)
			vcf_string += bp.as_vcf(ref_fasta)
			read_support_string += bp.as_read_support(count)
			count+=1
	# write output files
	with open(vcf_file, 'w') as output:
		output.write(vcf_string)
	with open(bedpe_file, 'w') as output:
		output.write(bedpe_string)
	with open(tsv_file, 'w') as output:
		output.write(read_support_string)
	# sort vcf
	bcftools.sort('-o', vcf_file, vcf_file, catch_stdout=False)

	helper.time_function("Output consensus breakpoints", checkpoints, time_str)

	return checkpoints, time_str

if __name__ == "__main__":
	print("Functions to run SAVANA")
