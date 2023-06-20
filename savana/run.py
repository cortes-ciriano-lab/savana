"""
SAVANA strucural variant caller for long-read data - run
Created: 17/02/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import os

from math import ceil
from multiprocessing import Pool

import pysam
import pysam.bcftools as bcftools
import pybedtools

import savana.helper as helper
from savana.breakpoints import get_potential_breakpoints, call_breakpoints, add_local_depth, add_local_depth_old
from savana.clusters import cluster_breakpoints

from memory_profiler import profile
from pympler import muppy, summary, refbrowser
import requests
import sys
import objgraph

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

def pool_add_local_depth_old(threads, breakpoint_dict_chrom, bam_files):
	pool_local_depth = Pool(processes=threads)
	pool_local_depth_args = []
	# convert bam_files into filenames (rather than objects - breaks parallelization)
	for label in bam_files.keys():
		bam_files[label] = bam_files[label].filename
	for chrom, breakpoints in breakpoint_dict_chrom.items():
		pool_local_depth_args.append((breakpoints, bam_files))
		local_depth_results = pool_local_depth.starmap(add_local_depth_old, pool_local_depth_args)
	pool_local_depth.close()
	pool_local_depth.join()
	breakpoint_by_chrom = {}
	for result in local_depth_results:
		chrom = result[0].start_chr
		if chrom not in breakpoint_by_chrom:
			breakpoint_by_chrom[chrom] = result
		else:
			breakpoint_by_chrom[chrom].extend(result)

	return breakpoint_by_chrom

def single_add_local_depth(intervals, bam_filenames):
	""" SINGLE THREAD OPTION
	given intervals and uids, get the local depth for each interval """
	if True:
		uid_dp_dict = {}
		chrom = intervals[0][0]
		start = int(intervals[0][1]) # first start
		end = int(intervals[-1][2]) # last end
		read_stats = {}
		for bam_type, bam_filename in bam_filenames.items():
			with pysam.AlignmentFile(bam_filename, "rb") as bam_file:
				bam_file = pysam.AlignmentFile(bam_filename, "rb")
				read_stats[bam_type] = []
				for read in bam_file.fetch(chrom, start, end):
					continue
				for read in bam_file.fetch(chrom, start, end):
					if read.mapping_quality == 0 or read.is_duplicate:
						continue
					read_stats[bam_type].append([int(read.reference_start), int(read.reference_length)])
				bam_file.close()
		#print(f'Calculating DP for {len(intervals)} Intervals')
		for i in intervals:
			uid = i[3]
			interval_end = int(i[2])
			edge = int(i[4])
			for bam_type, reads in read_stats.items():
				subtraction = [((interval_end-r[0])+r[1]) for r in reads]
				dp = sum(1 for x in subtraction if x >= 0)
				del subtraction
				if uid not in uid_dp_dict:
					uid_dp_dict[uid] = {}
				if bam_type not in uid_dp_dict[uid]:
					uid_dp_dict[uid][bam_type] = [None, None]
				uid_dp_dict[uid][bam_type][edge] = str(dp)
		#print(f'Done calculating DP for {len(intervals)} Intervals')

	else:
		# ALTERNATELY
		uid_dp_dict = {}
		for bam_type, bam_filename in bam_filenames.items():
			bam_file = pysam.AlignmentFile(bam_filename, "rb")
			for i in intervals:
				i_chrom, i_start, i_end, uid, edge = i
				reads = [read for read in bam_file.fetch(i_chrom, int(i_start), int(i_end))]
				reads = [read for read in reads if not read.is_duplicate and read.mapping_quality >= 0]
				if uid not in uid_dp_dict:
					uid_dp_dict[uid] = {}
				if bam_type not in uid_dp_dict[uid]:
					uid_dp_dict[uid][bam_type] = [None, None]
				uid_dp_dict[uid][bam_type][int(edge)] = str(len(reads))
			bam_file.close()

	return uid_dp_dict

def pool_add_local_depth(threads, sorted_bed, breakpoint_dict_chrom, bam_files):
	""" """
	from itertools import groupby
	from math import floor

	# OPTION FOR THREADS = 1
	if threads == 1:
		intervals_by_chrom = sorted([list(intervals) for _chrom, intervals in groupby(sorted_bed, lambda x: x[0])], key=len)
		local_depth_results = []
		for label in bam_files.keys():
			bam_files[label] = bam_files[label].filename
		for chrom_split in intervals_by_chrom:
			local_depth_results.append(single_add_local_depth(chrom_split, bam_files))
	else:
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
				# don't both splitting
				redistributed_intervals.append(chrom_chunk)
		print(f'Using {ideal_binsize} as binsize, there are {len(redistributed_intervals)} redistributed intervals')
		max_bin = max([len(c) for c in redistributed_intervals])
		min_bin = min([len(c) for c in redistributed_intervals])
		print(f'Max binsize {max_bin}, min binsize {min_bin}')

		# calculate max tasksperchild (max num)
		max_total_intervals_per_child = max(10000, max_bin+1) # figure this out by experiements
		max_tasks = floor(max_total_intervals_per_child/max_bin)
		print(f'Setting maxtasksperchild to {max_tasks}')

		pool_local_depth = Pool(processes=threads, maxtasksperchild=max_tasks)
		pool_local_depth_args = []
		# convert bam_files into filenames (rather than objects - breaks parallelization)
		for label in bam_files.keys():
			bam_files[label] = bam_files[label].filename
		for chrom_split in redistributed_intervals:
			pool_local_depth_args.append((chrom_split, bam_files))
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
			for bam_file, values in counts.items():
				if bam_file not in uid_dp_dict[uid]:
					uid_dp_dict[uid][bam_file] = [None, None]
				for i, dp in enumerate(values):
					uid_dp_dict[uid][bam_file][i] = dp if dp else uid_dp_dict[uid][bam_file][i]

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

def pool_call_breakpoints(threads, buffer, length, depth, clusters):
	""" parallelise the identification of consensus breakpoints """
	pool_calling = Pool(processes=threads)
	pool_calling_args = []

	for chrom, chrom_clusters in clusters.items():
		pool_calling_args.append((chrom_clusters, buffer, length, depth, chrom))
	calling_results = pool_calling.starmap(call_breakpoints, pool_calling_args)

	breakpoint_dict_chrom = {}
	seen_cluster_uids = {}
	pruned_clusters = {}
	for result_breakpoints, result_pruned_clusters, result_chrom in calling_results:
		# collect breakpoint calling results
		breakpoint_dict_chrom[result_chrom] = result_breakpoints
		# for debug
		#TODO: add a debug flag here, pass from args
		for bp_type in result_pruned_clusters.keys():
			for cluster in result_pruned_clusters[bp_type]:
				if cluster.uid not in seen_cluster_uids:
					pruned_clusters.setdefault(bp_type, []).append(cluster)
					seen_cluster_uids[cluster.uid] = True

	return breakpoint_dict_chrom, pruned_clusters

def profiling_experiment(intervals, bam_filenames):
	""" to identify and characterise potential memory leak in pysam """
	import gc
	import sys
	chrom = intervals[0][0]
	start = int(intervals[0][1]) # first start
	end = int(intervals[-1][2]) # last end
	# INVENTORY OF OBJECTS BEFORE
	objects_before = muppy.get_objects()
	print(f'Num. objects before: {len(objects_before)}')
	summary_before = summary.summarize(objects_before)
	objgraph.show_growth(limit=3)
	for bam_type, bam_filename in bam_filenames.items():
		with pysam.AlignmentFile(bam_filename, "rb") as bam_file:
			bam_file = pysam.AlignmentFile(bam_filename, "rb")
			# ITERATE THROUGH FETCHED READS
			count_reads = 0
			for read in bam_file.fetch(chrom, start, end):
				count_reads+=1
				continue
			#del read
			print(f'Num reads: {count_reads}')
		#del bam_file
	gc.collect()
	objgraph.show_growth(limit=3)
	# INVENTORY OF OBJECTS AFTER
	objects_after = muppy.get_objects()
	print(f'Num. objects after: {len(objects_after)}')
	summary_after = summary.summarize(objects_after)
	summary.print_(summary_after)
	print('Diff:')
	diff = summary.get_diff(summary_before, summary_after)
	summary.print_(diff)
	pysam_objs = muppy.filter(objects_after, Type=(pysam.libcalignedsegment.AlignedSegment))
	print(f'Num Aligned Segment Objects: {len(pysam_objs)}')
	for o in pysam_objs:
		cb = refbrowser.ConsoleBrowser(o, maxdepth=2, str_func=lambda x: str(type(x)))
		cb.print_tree()

def spawn_processes(args, bam_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using multiprocessing with {args.threads} threads\n')
	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints_results = pool_get_potential_breakpoints(bam_files, args)
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
	breakpoint_dict_chrom, pruned_clusters = pool_call_breakpoints(args.threads, args.buffer, args.length, args.depth, clusters)
	helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	total_breakpoints = 0
	for c, b in breakpoint_dict_chrom.items():
		total_breakpoints+=len(b)
	print(f'Length after: {total_breakpoints}')

	# skip for now for testing purposes
	if args.debug and False:
		# 3.1) OUTPUT CLUSTERS
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			pool_output_clusters(args, pruned_clusters[bp_type], outdir)
		helper.time_function("Output pruned clusters", checkpoints, time_str)

	# 5) ADD LOCAL DEPTH
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
	pool_add_local_depth(args.threads, sorted_bed, breakpoint_dict_chrom, bam_files)
	helper.time_function("Added local depth to breakpoints", checkpoints, time_str)

	# 6) OUTPUT BREAKPOINTS
	# define filenames
	vcf_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.vcf')
	bedpe_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.bedpe')
	tsv_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints_read_support.tsv')
	# build strings
	ref_fasta = pysam.FastaFile(args.ref)
	bedpe_string = ''
	vcf_string = helper.generate_vcf_header(args, breakpoint_dict_chrom[list(breakpoint_dict_chrom.keys())[0]][0])
	read_support_string = ''
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

	return pruned_clusters, breakpoint_dict_chrom, checkpoints, time_str

if __name__ == "__main__":
	print("Functions to run SAVANA")
