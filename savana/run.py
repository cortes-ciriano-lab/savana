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

import pysam
import pysam.bcftools as bcftools
import pybedtools

import savana.helper as helper
from savana.breakpoints import get_potential_breakpoints, call_breakpoints, add_local_depth
from savana.clusters import cluster_breakpoints, output_clusters

# developer dependencies
"""
from memory_profiler import profile
from pympler import muppy, summary, refbrowser
import objgraph
"""


def pool_get_potential_breakpoints(aln_files, args):
	""" split the genome into chunks and identify PotentialBreakpoints """
	pool_potential = Pool(processes=args.threads)
	pool_potential_args = []
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	if not args.is_cram:
		# calculate how to split contigs based on total mapped reads
		total_num_mapped_reads = 0
		for label, aln_file in aln_files.items():
			for contig in aln_file.get_index_statistics():
				if contig.contig in contigs_to_consider:
					total_num_mapped_reads+=contig.mapped
		ideal_reads_per_thread = ceil(total_num_mapped_reads/args.threads)
		# balance approx. number of reads per worker thread
		for label, aln_file in aln_files.items():
			for contig in aln_file.get_index_statistics():
				if contig.contig not in contigs_to_consider:
					if args.debug:
						print(f'Skipping reads aligned to {contig.contig} - not in contigs file')
					continue
				if contig.mapped == 0:
					continue
				mapped_reads = contig.mapped
				chrom_length = int(aln_file.get_reference_length(contig.contig))
				num_chunks = max(floor(mapped_reads/ideal_reads_per_thread), 1)
				if num_chunks > 1:
					start_pos = 0
					chunk_size = ceil(chrom_length/num_chunks)
					for i in range(1, num_chunks+1):
						end_pos = start_pos + chunk_size
						end_pos = chrom_length if end_pos > chrom_length else end_pos # don't extend past end
						pool_potential_args.append((aln_file.filename, args, label, contigs_to_consider, contig.contig, start_pos, end_pos))
						start_pos = end_pos + 1
				else:
					pool_potential_args.append((aln_file.filename, args, label, contigs_to_consider, contig.contig))
	else:
		# parallelize by contig (unable to see num. mapped reads per contig with cram)
		chunk_size = 60000000 # 60 million
		contig_lengths = helper.get_contig_lengths(args.ref_index)
		for label, aln_file in aln_files.items():
			for contig, contig_length in contig_lengths.items():
				if contig not in contigs_to_consider:
					continue
				if contig_length > chunk_size:
					# split the chrom into parts
					num_intervals = ceil(contig_length/chunk_size) + 1
					start_pos = 0
					for i in range(1, num_intervals):
						end_pos = start_pos + chunk_size
						end_pos = contig_length if end_pos > contig_length else end_pos # don't extend past end
						pool_potential_args.append((aln_file.filename, args, label, contigs_to_consider, contig, start_pos, end_pos))
						start_pos = end_pos + 1
				else:
					pool_potential_args.append((aln_file.filename, args, label, contigs_to_consider, contig))

	print(f'Submitting {len(pool_potential_args)} "get_potential_breakpoints" tasks to {args.threads} worker threads')

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

def pool_add_local_depth(threads, sorted_bed, breakpoint_dict_chrom, aln_files, is_cram=False, ref=False):
	""" """
	from itertools import groupby

	intervals_by_chrom = sorted([list(intervals) for _chrom, intervals in groupby(sorted_bed, lambda x: x[0])], key=len, reverse=True)
	total_length = sum([len(c) for c in intervals_by_chrom])
	redistributed_intervals = []
	#ideal_binsize = floor(total_length/(threads-2))
	#ideal_binsize = floor(total_length/(threads*2))
	ideal_binsize = max(floor(total_length/(threads*threads)),1)
	for chrom_chunk in intervals_by_chrom:
		if len(chrom_chunk) > 2*ideal_binsize:
			num_subchunks = floor(len(chrom_chunk)/ideal_binsize)
			# split list into equal chunks from https://stackoverflow.com/a/2135920
			quotient, remainder = divmod(len(chrom_chunk), num_subchunks)
			chunk_split = (chrom_chunk[i*quotient+min(i, remainder):(i+1)*quotient+min(i+1, remainder)] for i in range(num_subchunks))
			redistributed_intervals.extend(chunk_split)
		else:
			# don't bother splitting
			redistributed_intervals.append(chrom_chunk)
	print(f'Using {ideal_binsize} as binsize, there are {len(redistributed_intervals)} redistributed intervals')
	if not redistributed_intervals:
		import sys
		sys.exit('Issue calculating redistributed_intervals. Check input parameters')
	max_bin = max([len(c) for c in redistributed_intervals])
	min_bin = min([len(c) for c in redistributed_intervals])
	print(f'Max binsize {max_bin}, min binsize {min_bin}')

	# calculate max tasksperchild (max num)
	max_total_intervals_per_child = max(10000, max_bin+1) # figure this out by experiements
	max_tasks = floor(max_total_intervals_per_child/max_bin)
	print(f'Setting maxtasksperchild to {max_tasks}')

	pool_local_depth = Pool(processes=threads, maxtasksperchild=max_tasks)
	pool_local_depth_args = []
	# convert aln_files into filenames (rather than objects - breaks parallelization)
	for label in aln_files.keys():
		aln_files[label] = aln_files[label].filename
	for chrom_split in redistributed_intervals:
		pool_local_depth_args.append((chrom_split, aln_files, is_cram, ref))
	local_depth_results = pool_local_depth.starmap(add_local_depth, pool_local_depth_args)

	uid_dp_dict = {}
	"""
	dictionary structured as so:
	{
		'uid': {'TUMOUR_DP': [0,0], 'NORMAL_DP': [0,0]}
	}
	"""
	for result in local_depth_results:
		for uid, counts in result.items():
			if uid not in uid_dp_dict:
				uid_dp_dict[uid] = {}
			for aln_file, values in counts.items():
				if aln_file not in uid_dp_dict[uid]:
					uid_dp_dict[uid][aln_file] = [None, None]
				for i, dp in enumerate(values):
					uid_dp_dict[uid][aln_file][i] = dp if dp else uid_dp_dict[uid][aln_file][i]

	for breakpoints in breakpoint_dict_chrom.values():
		for bp in breakpoints:
			if bp.breakpoint_notation == "<INS>":
				# remove None values for insertions
				reformatted_dp = {}
				for t,d in uid_dp_dict[bp.uid].items():
					reformatted_dp[t] = [d[0]]
				bp.local_depths = reformatted_dp
			else:
				bp.local_depths = uid_dp_dict[bp.uid]

	return breakpoints

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

def spawn_processes(args, aln_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using multiprocessing with {args.threads} threads\n')
	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints_results = pool_get_potential_breakpoints(aln_files, args)
	helper.time_function("Identified potential breakpoints", checkpoints, time_str)
	# collect results per chrom
	chrom_potential_breakpoints = {}
	for result in potential_breakpoints_results:
		for chrom, potential_breakpoints in result.items():
			chrom_potential_breakpoints.setdefault(chrom,[]).extend(potential_breakpoints)

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clusters = pool_cluster_breakpoints(args.threads, args.buffer, args.insertion_buffer, chrom_potential_breakpoints)
	helper.time_function("Clustered potential breakpoints", checkpoints, time_str)

	# 3) CALL BREAKPOINTS FROM CLUSTERS
	breakpoint_dict_chrom, pruned_clusters = pool_call_breakpoints(args.threads, args.buffer, args.length, args.depth, clusters, args.debug)
	helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	total_breakpoints = 0
	for c, b in breakpoint_dict_chrom.items():
		total_breakpoints+=len(b)
	print(f'Length after: {total_breakpoints}')

	if args.debug:
		# 3.1) OUTPUT CLUSTERS
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			if bp_type in pruned_clusters:
				pool_output_clusters(args, pruned_clusters[bp_type], outdir)
		helper.time_function("Output pruned clusters", checkpoints, time_str)

	# 4) ADD LOCAL DEPTH
	# generate interval files
	bed_string = ''
	total_num_breakpoints = 0
	total_num_insertions = 0
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	for chrom, chrom_breakpoints in breakpoint_dict_chrom.items():
		total_num_breakpoints+=len(chrom_breakpoints)
		for bp in chrom_breakpoints:
			if bp.breakpoint_notation == "<INS>":
				total_num_insertions += 1
			bed_string += bp.as_bed(contig_lengths)
	sorted_bed = pybedtools.BedTool(bed_string, from_string=True).sort(faidx=args.ref_index)
	print(f'Total breakpoints: {total_num_breakpoints} ({total_num_insertions} insertions)')
	pool_add_local_depth(args.threads, sorted_bed, breakpoint_dict_chrom, aln_files, args.is_cram, args.ref)
	helper.time_function("Added local depth to breakpoints", checkpoints, time_str)

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
