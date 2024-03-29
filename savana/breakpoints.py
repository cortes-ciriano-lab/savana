"""
Module containing functions related to Breakpoint classes for SAVANA
Created: 06/09/2022
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

from statistics import median
import pysam
import numpy as np

from math import floor, ceil

import savana.helper as helper
from savana.core import PotentialBreakpoint, ConsensusBreakpoint, Cluster

def get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order):
	""" reconstruct the breakpoints from the supplementary alignments """
	primary_clipping = helper.get_clipping(cigar_tuples, read.is_reverse)
	breakpoint_pairs = []
	# sort the chimeric regions by the left soft clip pos (ascending)
	chimeric_regions = sorted(chimeric_regions, key=lambda d: d['left_softclip'])
	left_index = 0
	# in the left softclip
	while left_index < len(chimeric_regions) and chimeric_regions[left_index]['left_softclip'] < primary_clipping['left_softclip']:
		if left_index != 0:
			# for all others except first
			# add the second edge of the supp. breakpoint
			bp_end = int(chimeric_regions[left_index]['pos'])
			if chimeric_regions[left_index]['strand'] == '-':
				bp_end += chimeric_regions[left_index]['consumed_reference']
				bp_notation = ("+")
			else:
				bp_notation = ("-")
			# add to existing breakpoint
			breakpoint_pairs[-1].append({
				'chr': chimeric_regions[left_index]['chrom'],
				'loc': bp_end,
				'bp_notation': bp_notation
			})
		# create the first edge of each supp. breakpoint in the left softclip
		bp_start = int(chimeric_regions[left_index]['pos'])
		if chimeric_regions[left_index]['strand'] == '+':
			bp_start += chimeric_regions[left_index]['consumed_reference']
			bp_notation = ("+")
		else:
			bp_notation = ("-")
		breakpoint_pairs.append([{
			'chr': chimeric_regions[left_index]['chrom'],
			'loc': bp_start,
			'bp_notation': bp_notation
		}])
		left_index+=1
	if breakpoint_pairs:
		# add the primary alignment as the second edge on the last breakpoint
		breakpoint_pairs[-1].append({
			'chr': read.reference_name,
			'loc': read.reference_end if read.is_reverse else read.reference_start,
			'bp_notation': ("+") if read.is_reverse else ("-")
		})
	# if still have chimeric regions, need to place in right softclip
	if left_index < len(chimeric_regions):
		# add the primary as the first edge
		breakpoint_pairs.append([{
			'chr': read.reference_name,
			'loc': read.reference_start if read.is_reverse else read.reference_end,
			'bp_notation': ("-") if read.is_reverse else ("+")
		}])
		# sort remaining chimeric regions by the right soft clip (descending)
		right_chimeric = sorted(chimeric_regions[left_index:], key=lambda d: d['right_softclip'], reverse=True)
		right_index = 0
		while right_index < len(right_chimeric) and right_chimeric[right_index]['right_softclip'] < primary_clipping['right_softclip']:
			# supp as second edge
			bp_end = int(right_chimeric[right_index]['pos'])
			if right_chimeric[right_index]['strand'] == '-':
				bp_end += right_chimeric[right_index]['consumed_reference']
				bp_notation = ("+")
			else:
				bp_notation = ("-")
			# add to existing breakpoint
			breakpoint_pairs[-1].append({
				'chr': right_chimeric[right_index]['chrom'],
				'loc': bp_end,
				'bp_notation': bp_notation
			})
			if right_index != len(right_chimeric) - 1:
				# first edge for supp
				# (as long as not last in right chimeric regions)
				bp_start = int(right_chimeric[right_index]['pos'])
				if right_chimeric[right_index]['strand'] == "+":
					bp_start += right_chimeric[right_index]['consumed_reference']
					bp_notation = ("+")
				else:
					bp_notation = ("-")
				breakpoint_pairs.append([{
					'chr': right_chimeric[right_index]['chrom'],
					'loc': bp_start,
					'bp_notation': bp_notation
				}])
			right_index+=1
	# once all pairs completed, create breakpoints for each edge
	supplementary_breakpoints = []
	for start, end in breakpoint_pairs:
		if start['chr'] not in contig_order or end['chr'] not in contig_order:
			continue
		if start['chr'] == end['chr']:
			if start['loc'] < end['loc']:
				supplementary_breakpoints.append(PotentialBreakpoint([start, end], "SUPP", read.query_name, read.mapping_quality, label, "".join((start['bp_notation'], end['bp_notation']))))
			else:
				supplementary_breakpoints.append(PotentialBreakpoint([end, start], "SUPP", read.query_name, read.mapping_quality, label, "".join((end['bp_notation'], start['bp_notation']))))
		elif contig_order.index(start['chr']) <= contig_order.index(end['chr']):
			supplementary_breakpoints.append(PotentialBreakpoint([start, end], "SUPP", read.query_name, read.mapping_quality, label, "".join((start['bp_notation'], end['bp_notation']))))
		elif start['loc'] < end['loc']:
			supplementary_breakpoints.append(PotentialBreakpoint([end, start], "SUPP", read.query_name, read.mapping_quality, label, "".join((end['bp_notation'], start['bp_notation']))))
	return supplementary_breakpoints

def count_num_labels(source_breakpoints):
	""" given a list of unique breakpoints, return the counts for each label """
	label_counts = {}
	seen_reads = {}
	for bp in source_breakpoints:
		if bp.read_name in seen_reads:
			continue
		if bp.label not in label_counts:
			label_counts[bp.label] = [bp.read_name]
		else:
			label_counts[bp.label].append(bp.read_name)
		seen_reads[bp.read_name] = True

	return label_counts

def get_potential_breakpoints(aln_filename, args, label, contig_order, chrom=None, start=None, end=None):
	""" iterate through alignment file, tracking potential breakpoints and saving relevant reads to fastq """
	potential_breakpoints = {}
	if args.is_cram:
		aln_file = pysam.AlignmentFile(aln_filename, "rc", reference_filename=args.ref)
	else:
		aln_file = pysam.AlignmentFile(aln_filename, "rb")
	# adjust the thresholds depending on sample source
	args_length = max((args.length - floor(args.length/5)), 0) if label == 'normal' else args.length
	mapq = min((args.mapq - ceil(args.mapq/2)), 1) if label == 'normal' else args.mapq
	for read in aln_file.fetch(chrom, start, end):
		if read.is_secondary or read.is_supplementary:
			continue # only consider primary
		if read.mapping_quality < mapq:
			continue # discard if mapping quality lower than threshold
		curr_pos = {
			'cigar': 0,
			'reference': read.reference_start,
			'query': 0
		}
		cigar_tuples = read.cigartuples
		chimeric_regions = helper.get_chimeric_regions(read, mapq)
		if chimeric_regions:
			chimeric_breakpoints = get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order)
			for bp in chimeric_breakpoints:
				potential_breakpoints.setdefault(bp.start_chr,[]).append(bp)
		# look for insertions and deletions in the CIGAR
		curr_chrom = read.reference_name
		prev_deletion = []
		for sam_flag, length in cigar_tuples:
			if sam_flag == helper.samflag_desc_to_number["BAM_CINS"] and length > args_length:
				# start and end location for insertion breakpoint are the same
				location = [
					{'chr': curr_chrom, 'loc': curr_pos['reference']},
					{'chr': curr_chrom, 'loc': curr_pos['reference']}
				]
				if prev_deletion:
					# record and clear the previously tracked deletion
					potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "DEL", read.query_name, read.mapping_quality, label, "+-"))
					prev_deletion = []
				# record the insertion
				inserted_sequence = read.query_sequence[curr_pos['query']:(curr_pos['query']+length)]
				potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(location, "INS", read.query_name, read.mapping_quality, label, "<INS>", inserted_sequence))
			elif sam_flag == helper.samflag_desc_to_number["BAM_CDEL"] and length > args_length:
				# deletion has one breakpoint (read->read)
				if prev_deletion:
					# expand the deletion end location to end of current deletion
					prev_deletion[1]['loc'] = curr_pos['reference']+length
				else:
					prev_deletion = [
						{'chr': curr_chrom, 'loc': curr_pos['reference']},
						{'chr': curr_chrom, 'loc': curr_pos['reference']+length}
					]
			elif prev_deletion and length > args_length:
				# record and clear the previously tracked deletion
				potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "DEL", read.query_name, read.mapping_quality, label, "+-"))
				prev_deletion = []
			# increment values
			curr_pos['cigar'] += length
			if helper.consumes_query[sam_flag]:
				curr_pos['query'] += length
			if helper.consumes_reference[sam_flag]:
				curr_pos['reference'] += length
		if prev_deletion:
			# if reached end of string and no chance to expand deletion, add it
			potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "DEL", read.query_name, read.mapping_quality, label, "+-"))

	aln_file.close()

	return potential_breakpoints

def add_local_depth(intervals, aln_filenames, is_cram, ref):
	""" given intervals and uids, get the local depth for each interval """
	uid_dp_dict = {}
	chrom = intervals[0][0]
	start = max(int(intervals[0][1])-1, 0) # first start
	end = int(intervals[-1][2]) # last end
	read_stats = {}
	for file_type, aln_filename in aln_filenames.items():
		if is_cram:
			aln_file = pysam.AlignmentFile(aln_filename, "rc", reference_filename=ref)
		else:
			aln_file = pysam.AlignmentFile(aln_filename, "rb")
		read_stats[file_type] = []
		for read in aln_file.fetch(chrom, start, end):
			if read.mapping_quality == 0 or read.is_duplicate:
				continue
			read_stats[file_type].append([int(read.reference_start), int(read.reference_end), read.query_name])
			del read
		aln_file.close()
		del aln_file
	for i in intervals:
		uid = i[3]
		interval_start = int(i[1])
		interval_end = int(i[2])
		edge = int(i[4])
		for file_type, reads in read_stats.items():
			comparison = [[(interval_start - r[1]), (r[0] - interval_end)] for r in reads]
			dp = sum(1 for x,y in comparison if x <= 0 and y <= 0)
			del comparison
			# for some reason, these methods aren't faster
			'''
			comparison = [[(interval_start <= r[1]), (r[0] <= interval_end)] for r in reads]
			dp = sum(1 for x,y in comparison if x and y)
			del comparison
			'''
			# nor
			#dp = sum(1 for r in reads if (interval_start - r[1]) <= 0 and (r[0] - interval_end) <= 0)
			if uid not in uid_dp_dict:
				uid_dp_dict[uid] = {}
			if file_type not in uid_dp_dict[uid]:
				uid_dp_dict[uid][file_type] = [None, None]
			uid_dp_dict[uid][file_type][edge] = str(dp)

	return uid_dp_dict

def call_breakpoints(clusters, buffer, min_length, min_depth, chrom):
	""" identify consensus breakpoints from list of clusters """
	# N.B. all breakpoints in a cluster must be from same chromosome!
	final_breakpoints = []
	pruned_clusters = {}
	num_insertions = 0
	for bp_type in clusters.keys():
		for cluster in clusters[bp_type]:
			if bp_type == "<INS>":
				# call validated insertions
				# average out the start/end, keep longest insertion sequence
				starts, ends = [], []
				longest_seq = ''
				for bp in cluster.breakpoints:
					starts.append(bp.start_loc)
					ends.append(bp.end_loc)
					if len(bp.inserted_sequence) > len(longest_seq):
						longest_seq = bp.inserted_sequence
				source_breakpoints = cluster.breakpoints
				label_counts = count_num_labels(source_breakpoints)
				if max([len(v) for v in label_counts.values()]) >= min_depth and len(longest_seq) > min_length:
					num_insertions += 1
					final_breakpoints.append(ConsensusBreakpoint(
						[{'chr': cluster.chr, 'loc': median(starts)}, {'chr': cluster.chr, 'loc': median(ends)}],
						"INS", cluster, None, label_counts, bp_type, longest_seq))
					pruned_clusters.setdefault(bp_type, []).append(cluster)
			else:
				# call all other types
				per_end_chrom = {}
				# separate breakpoints by end chrom
				for bp in cluster.breakpoints:
					if bp.end_chr not in per_end_chrom:
						per_end_chrom[bp.end_chr] = {
							'starts': [bp.start_loc],
							'ends': [bp.end_loc],
							'breakpoints': [bp],
							'originating_cluster': cluster
						}
					else:
						per_end_chrom[bp.end_chr]['starts'].append(bp.start_loc)
						per_end_chrom[bp.end_chr]['ends'].append(bp.end_loc)
						per_end_chrom[bp.end_chr]['breakpoints'].append(bp)
				# cluster by end location
				for _, end_chrom_info in per_end_chrom.items():
					# flip original breakpoints
					source_breakpoints = [reversed(b) for b in end_chrom_info['breakpoints']]
					if len(source_breakpoints) > 1:
						# sort by "start" - which is actually end
						source_breakpoints.sort()
					# cluster sorted breakpoints
					cluster_stack = []
					for bp in source_breakpoints:
						new_cluster = False
						if len(cluster_stack) == 0 or not abs(cluster_stack[-1].start - bp.start_loc) < buffer:
							# put a new cluster onto top of stack
							new_cluster = Cluster(bp)
							cluster_stack.append(new_cluster)
							new_cluster = True
						else:
							# append to cluster on top of stack
							cluster_stack[-1].add(bp)
					for end_cluster in cluster_stack:
						# create ConsensusBreakpoint per end_cluster
						end_cluster_breakpoints = end_cluster.breakpoints
						label_counts = count_num_labels(end_cluster_breakpoints)
						median_start = median([bp.end_loc for bp in end_cluster_breakpoints])
						median_end = median([bp.start_loc for bp in end_cluster_breakpoints])
						# need to create a new "originating cluster" with only the used breakpoints
						# this ensures that the statistics are calculated correctly
						new_start_cluster = None
						for bp in end_cluster_breakpoints:
							if not new_start_cluster:
								new_start_cluster = Cluster(reversed(bp))
							else:
								new_start_cluster.add(reversed(bp))
						if max([len(v) for v in label_counts.values()]) >= min_depth:
							new_breakpoint = ConsensusBreakpoint(
								[{'chr': cluster.chr, 'loc': median_start}, {'chr': end_cluster.chr, 'loc': median_end}],
								bp_type, new_start_cluster, end_cluster, label_counts, bp_type)
							if new_breakpoint.sv_length >= min_length or new_breakpoint.sv_length == 0:
								final_breakpoints.append(new_breakpoint)
								pruned_clusters.setdefault(bp_type, []).append(new_start_cluster)
								pruned_clusters.setdefault(bp_type, []).append(end_cluster)

	return final_breakpoints, pruned_clusters, chrom

if __name__ == "__main__":
	print("Breakpoint Functions")
