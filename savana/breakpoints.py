"""
Module containing functions related to Breakpoint classes for SAVANA
Created: 06/09/2022
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import pysam
import numpy as np

from math import floor, ceil
from statistics import median

import savana.helper as helper
from savana.core import PotentialBreakpoint, ConsensusBreakpoint, Cluster

"""
# developer dependencies
from memory_profiler import profile
from pympler import muppy, summary, refbrowser

# decorate functions like so:
#helper.conditionally_decorate(profile, True)
"""

#TODO: investigate/test if passing relevant read attributes to this function (rather than entire object), improves speed/mem
def get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order, single_bnd, single_bnd_min_length, single_bnd_max_mapq):
	""" reconstruct the breakpoints from the supplementary alignments """
	breakpoint_pairs = []
	primary_clipping = helper.get_clipping(cigar_tuples, read.is_reverse)
	#	primary_clipping = helper.get_clipping(cigar_tuples, False) sort the chimeric regions by the left soft clip pos (ascending)
	chimeric_regions = sorted(chimeric_regions, key=lambda d: d['left_softclip'])
	left_index = 0
	# track the position in the read in order to extract sbnd bases
	curr_pos_query = min(chimeric_regions[left_index]['left_softclip'], primary_clipping['left_softclip'])
	while left_index < len(chimeric_regions) and chimeric_regions[left_index]['left_softclip'] < primary_clipping['left_softclip']:
		# in the left softclip
		region_mapped = False if single_bnd and (chimeric_regions[left_index]['mapQ'] <= single_bnd_max_mapq or chimeric_regions[left_index]['chrom'] not in contig_order) else True
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
				'bp_notation': bp_notation,
				'region_mapped': region_mapped,
				'query_pos': chimeric_regions[left_index]['left_softclip'] + chimeric_regions[left_index-1]['consumed_query']
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
			'bp_notation': bp_notation,
			'region_mapped': region_mapped,
			'query_pos': chimeric_regions[left_index]['left_softclip'] # curr_pos_query
		}])
		# increment position in the query
		curr_pos_query += chimeric_regions[left_index]['consumed_query']
		left_index+=1
	if breakpoint_pairs:
		# add the primary alignment as the second edge on the last breakpoint
		breakpoint_pairs[-1].append({
			'chr': read.reference_name,
			'loc': read.reference_end if read.is_reverse else read.reference_start,
			'bp_notation': ("+") if read.is_reverse else ("-"),
			'region_mapped': False if single_bnd and read.mapq <= single_bnd_max_mapq else True,
			'query_pos': primary_clipping['left_softclip'] + read.query_alignment_length # curr_pos_query
		})
	# if still have chimeric regions, need to place in right softclip
	if left_index < len(chimeric_regions):
		# add the primary as the first edge
		breakpoint_pairs.append([{
			'chr': read.reference_name,
			'loc': read.reference_start if read.is_reverse else read.reference_end,
			'bp_notation': ("-") if read.is_reverse else ("+"),
			'region_mapped': False if single_bnd and read.mapq <= single_bnd_max_mapq else True,
			'query_pos': primary_clipping['left_softclip'] # curr_pos_query
		}])
		#curr_pos_query += read.query_alignment_length
		curr_pos_query += helper.sum_consumed_query(helper.trim_supplementary(read.cigarstring))
		# sort remaining chimeric regions by the right soft clip (descending)
		# TODO: don't need to resort this can just use old one
		#right_chimeric = sorted(chimeric_regions[left_index:], key=lambda d: d['right_softclip'], reverse=True)
		right_chimeric = sorted(chimeric_regions[left_index:], key=lambda d: d['left_softclip'])
		right_index = 0
		while right_index < len(right_chimeric) and right_chimeric[right_index]['right_softclip'] < primary_clipping['right_softclip']:
			region_mapped = False if single_bnd and (right_chimeric[right_index]['mapQ'] <= single_bnd_max_mapq or right_chimeric[right_index]['chrom'] not in contig_order) else True
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
				'bp_notation': bp_notation,
				'region_mapped': region_mapped,
				'query_pos': right_chimeric[right_index]['left_softclip'] + right_chimeric[right_index]['consumed_query'] # curr_pos_query
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
					'bp_notation': bp_notation,
					'region_mapped': region_mapped,
					'query_pos': right_chimeric[right_index]['left_softclip'] # curr_pos_query
				}])
				# increment position in the query
				curr_pos_query += right_chimeric[right_index]['consumed_query']
			right_index+=1

	# once all pairs completed, create breakpoints for each edge
	return create_breakpoint_objects(read, label, contig_order, breakpoint_pairs, chimeric_regions)

def create_breakpoint_objects(read, label, contig_order, breakpoint_pairs, chimeric_regions):
	supplementary_breakpoints = []
	index = 0 # track manually to allow skipping/merging
	while index < len(breakpoint_pairs):
		curr_start, curr_end = breakpoint_pairs[index]
		# don't create breakpoints for irrelevant contigs
		if curr_start['chr'] not in contig_order or curr_end['chr'] not in contig_order:
			index += 1
			continue
		if not curr_end['region_mapped']:
			# consume until mapped region found (merging if multiple unmapped)
			found_mapped = False
			next_index = index + 1
			while next_index < len(breakpoint_pairs) and not found_mapped:
				next_start, next_end = breakpoint_pairs[next_index]
				if next_start['region_mapped']:
					# segment mapped, create single breakend
					found_mapped = True
					# get sequence for query - from curr_start to next_start
					sbnd_seq = read.query_sequence[curr_start['query_pos']:(next_start['query_pos'])]
					# only one location - next start since it's mapped
					location = [{'chr': next_start['chr'], 'loc': next_start['loc']}, {'chr': next_start['chr'], 'loc': next_start['loc']}]
					if len(sbnd_seq) > 0:
						supplementary_breakpoints.append(PotentialBreakpoint(
							location,
							"SUPPLEMENTARY",
							read.query_name,
							0, # TODO: what to set for MAPQ?
							label,
							"<SBND>",
							sbnd_seq
						))
					else:
						print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos"]}:{next_start["query_pos"]} (curr_start - next_start)')
						print(breakpoint_pairs)
						print(chimeric_regions)
				elif next_end['region_mapped']:
					# segment mapped, create single breakend
					found_mapped = True
					# get sequence for query - from curr_start to next_end
					sbnd_seq = read.query_sequence[curr_start['query_pos']:(next_end['query_pos'])]
					location = [{'chr': next_end['chr'], 'loc': next_end['loc']}, {'chr': next_end['chr'], 'loc': next_end['loc']}]
					if len(sbnd_seq) > 0:
						supplementary_breakpoints.append(PotentialBreakpoint(
							location,
							"SUPPLEMENTARY",
							read.query_name,
							0, # TODO: what to set for MAPQ?
							label,
							"<SBND>",
							sbnd_seq
						))
					else:
						print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos"]}:{next_end["query_pos"]} (curr_start - next_end)')
						print(breakpoint_pairs)
						print(chimeric_regions)
				next_index += 1
			index = next_index - 1
			# got to end and didn't find map, means we should use the curr_start to last end
			if not found_mapped:
				_, last_end = breakpoint_pairs[-1]
				# get sequence for query - from curr_start to next_start
				sbnd_seq = read.query_sequence[curr_start['query_pos']:(last_end['query_pos'])]
				# only one location - can't be end since isn't mapped
				location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']}, {'chr': curr_start['chr'], 'loc': curr_start['loc']}]
				if len(sbnd_seq) > 0:
					supplementary_breakpoints.append(PotentialBreakpoint(
						location,
						"SUPPLEMENTARY",
						read.query_name,
						0, # TODO: what to set for MAPQ?
						label,
						"<SBND>",
						sbnd_seq
					))
				else:
					print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos"]}:{last_end["query_pos"]} (curr_start - last_end)')
					print(breakpoint_pairs)
					print(chimeric_regions)
		elif not curr_start['region_mapped']:
			# end is mapped, no need to merge/consume next pair
			sbnd_seq = read.query_sequence[curr_start['query_pos']:(curr_end['query_pos'])]
			# only one location - the mapped end
			location = [{'chr': curr_end['chr'], 'loc': curr_end['loc']}, {'chr': curr_end['chr'], 'loc': curr_end['loc']}]
			if len(sbnd_seq) > 0:
				supplementary_breakpoints.append(PotentialBreakpoint(
					location,
					"SUPPLEMENTARY",
					read.query_name,
					0, # TODO: what to set for MAPQ?
					label,
					"<SBND>",
					sbnd_seq
				))
			else:
				print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos"]}:{curr_end["query_pos"]} (curr_start - curr_end)')
				print(breakpoint_pairs)
				print(chimeric_regions)
		else:
			# boring region, everything mapped
			location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']}, {'chr': curr_end['chr'], 'loc': curr_end['loc']}]
			if curr_start['chr'] == curr_end['chr']:
				if curr_start['loc'] < curr_end['loc']:
					supplementary_breakpoints.append(PotentialBreakpoint(
						location,
						"SUPPLEMENTARY",
						read.query_name,
						read.mapping_quality,
						label,
						"".join((curr_start['bp_notation'],curr_end['bp_notation']))
					))
				else:
					supplementary_breakpoints.append(PotentialBreakpoint(
						list(reversed(location)),
						"SUPPLEMENTARY",
						read.query_name,
						read.mapping_quality,
						label,
						"".join((curr_end['bp_notation'], curr_start['bp_notation']))
					))
			elif contig_order.index(curr_start['chr']) <= contig_order.index(curr_end['chr']):
				supplementary_breakpoints.append(PotentialBreakpoint(
					location,
					"SUPPLEMENTARY",
					read.query_name,
					read.mapping_quality,
					label,
					"".join((curr_start['bp_notation'], curr_end['bp_notation']))
				))
			elif curr_start['loc'] < curr_end['loc']:
				supplementary_breakpoints.append(PotentialBreakpoint(
					list(reversed(location)),
					"SUPPLEMENTARY",
					read.query_name,
					read.mapping_quality,
					label,
					"".join((curr_end['bp_notation'], curr_start['bp_notation']))
				))
		index += 1

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

def get_potential_breakpoints(aln_filename, is_cram, ref, length, mapq, label, contig_order, contig, start, end, coverage_binsize, contig_coverage_array):
	""" iterate through alignment file, tracking potential breakpoints and saving relevant reads to fastq """
	potential_breakpoints = {}
	aln_file = pysam.AlignmentFile(aln_filename, "rc", reference_filename=ref) if is_cram else pysam.AlignmentFile(aln_filename, "rb")
	# adjust the thresholds depending on sample source
	args_length = max((length - floor(length/5)), 0) if label == 'normal' else length
	mapq = min((mapq - ceil(mapq/2)), 0) if label == 'normal' else mapq
	#TODO: make these args
	single_bnd = True
	single_bnd_min_length = 100
	single_bnd_max_mapq = 20
	####
	for read in aln_file.fetch(contig, start, end):
		if read.is_secondary or read.is_unmapped:
			continue
		try:
			# record start/end in read incrementer
			contig_coverage_array[floor((read.reference_start-1)/coverage_binsize)]+=1
			contig_coverage_array[floor((read.reference_end-1)/coverage_binsize)]-=1
		except IndexError as e:
			print(f'Unable to update coverage for contig {contig}')
			print(f'Attempting to update bins {floor((read.reference_start-1)/coverage_binsize)} and {floor((read.reference_end-1)/coverage_binsize)}')
			print(f'Length of array {len(contig_coverage_array)}, Read {read.reference_start} to {read.reference_end}')
		if read.mapping_quality < mapq:
			continue # discard if mapping quality lower than threshold
		curr_pos = {
			'cigar': 0,
			'reference': read.reference_start,
			'query': 0
		}
		cigar_tuples = read.cigartuples
		chimeric_regions = helper.get_chimeric_regions(read)
		if chimeric_regions:
			chimeric_breakpoints = get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order, single_bnd, single_bnd_min_length, single_bnd_max_mapq)
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
					potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-"))
					prev_deletion = []
				# record the insertion
				inserted_sequence = read.query_sequence[curr_pos['query']:(curr_pos['query']+length)]
				potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(location, "CIGAR", read.query_name, read.mapping_quality, label, "<INS>", inserted_sequence))
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
			elif sam_flag == helper.samflag_desc_to_number["BAM_CSOFT_CLIP"] and single_bnd and not chimeric_regions and length > single_bnd_min_length:
				# start and end location for single breakend are the same (for clustering purposes)
				location = [
					{'chr': curr_chrom, 'loc': curr_pos['reference']},
					{'chr': curr_chrom, 'loc': curr_pos['reference']}
				]
				# record the softclipped bases
				softclipped_sequence = read.query_sequence[curr_pos['query']:(curr_pos['query']+length)]
				if len(softclipped_sequence) > 0:
					potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(location, "CIGAR", read.query_name, read.mapping_quality, label, "<SBND>", softclipped_sequence))
				else:
					print(f'softclipped sequence empty! (from {curr_pos["query"]}-{curr_pos["query"]+length})')
			elif prev_deletion and length > args_length:
				# record and clear the previously tracked deletion
				potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-"))
				prev_deletion = []
			# increment values
			curr_pos['cigar'] += length
			if helper.consumes_query[sam_flag]:
				curr_pos['query'] += length
			if helper.consumes_reference[sam_flag]:
				curr_pos['reference'] += length
		if prev_deletion:
			# if reached end of string and no chance to expand deletion, add it
			potential_breakpoints.setdefault(curr_chrom,[]).append(PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-"))

	aln_file.close()
	del aln_file

	return potential_breakpoints

def call_breakpoints(clusters, end_buffer, min_length, min_support, chrom):
	""" identify consensus breakpoints from list of clusters """
	# N.B. all breakpoints in a cluster must be from same chromosome!
	final_breakpoints = []
	for bp_type in clusters.keys():
		for cluster in clusters[bp_type]:
			if bp_type == "<INS>":
				# call validated insertions
				# average out the start/end
				starts, ends = [], []
				inserts = []
				for bp in cluster.breakpoints:
					starts.append(bp.start_loc)
					ends.append(bp.end_loc)
					inserts.append(bp.inserted_sequence)
				source_breakpoints = cluster.breakpoints
				label_counts = count_num_labels(source_breakpoints)
				if max([len(v) for v in label_counts.values()]) >= min_support:
					final_breakpoints.append(ConsensusBreakpoint(
						[{'chr': cluster.chr, 'loc': median(starts)}, {'chr': cluster.chr, 'loc': median(ends)}],
						"INS", cluster, None, label_counts, bp_type, inserts))
					pruned_clusters.setdefault(bp_type, []).append(cluster)
			elif bp_type == "<SBND>":
				starts, ends = [], []
				inserts = []
				for bp in cluster.breakpoints:
					starts.append(bp.start_loc)
					ends.append(bp.end_loc)
					inserts.append(bp.inserted_sequence)
				source_breakpoints = cluster.breakpoints
				label_counts = count_num_labels(source_breakpoints)
				if max([len(v) for v in label_counts.values()]) >= min_support:
					starts, inserts, sources = [], [], []
					for bp in cluster_breakpoints:
						starts.append(bp.start_loc)
						inserts.append(bp.inserted_sequence)
						sources.append(bp.source)
					median_start = median(starts)
					# can be from softclip or supplementary
					sources = set(sources)
					consensus_source = "CIGAR/SUPPLEMENTARY" if len(sources) > 1 else sources.pop()
					final_breakpoints.append(ConsensusBreakpoint(
						[{'chr': cluster.chr, 'loc': median_start}, {'chr': cluster.chr, 'loc': median_start}],
						consensus_source, cluster, None, label_counts, bp_type, inserts))
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
						if len(cluster_stack) == 0 or not abs(cluster_stack[-1].start - bp.start_loc) < end_buffer:
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
						if max([len(v) for v in label_counts.values()]) >= min_support:
							# need to create a new "originating cluster" with only the used breakpoints
							# this ensures that the statistics are calculated correctly
							new_start_cluster = None
							consensus_source = [] if bp_type == "+-" else end_cluster_breakpoints[0].source
							for bp in end_cluster_breakpoints:
								if bp_type == "+-":
									consensus_source.append(bp.source)
								if not new_start_cluster:
									new_start_cluster = Cluster(reversed(bp))
								else:
									new_start_cluster.add(reversed(bp))
							if bp_type == "+-":
								# can be from CIGAR or SUPPLEMENTARY
								consensus_source = set(consensus_source)
								consensus_source = "CIGAR/SUPPLEMENTARY" if len(consensus_source) > 1 else consensus_source.pop()
							median_start = median([bp.end_loc for bp in end_cluster_breakpoints])
							median_end = median([bp.start_loc for bp in end_cluster_breakpoints])
							new_breakpoint = ConsensusBreakpoint(
								[{'chr': cluster.chr, 'loc': median_start}, {'chr': end_cluster.chr, 'loc': median_end}],
								consensus_source, new_start_cluster, end_cluster, label_counts, bp_type)
							if new_breakpoint.sv_length >= min_length or new_breakpoint.sv_length == 0:
								final_breakpoints.append(new_breakpoint)

	return final_breakpoints, chrom

def compute_depth(breakpoints, shared_cov_arrays, coverage_binsize):
	""" use the contig coverages to annotate the depth of breakpoints (same method as mosdepth) """
	for bp in breakpoints:
		bp.local_depths = {
			'tumour': [[0,0],[0,0],[0,0]],
			'normal': [[0,0],[0,0],[0,0]],
		}
		for label in bp.local_depths.keys():
			for i, (chrom, loc) in enumerate([(bp.start_chr, bp.start_loc),(bp.end_chr, bp.end_loc)]):
				centre_bin = floor(loc/coverage_binsize)
				if centre_bin > 0:
					bp.local_depths[label][0][i] = np.sum(shared_cov_arrays[label][chrom][0:centre_bin-1])
				else:
					# unlikely case it's the first bin
					bp.local_depths[label][0][i] = 0
				bp.local_depths[label][1][i] = bp.local_depths[label][0][i] + shared_cov_arrays[label][chrom][centre_bin]
				try:
					# handle extremely unlikely case that there's no next bin
					bp.local_depths[label][2][i] = bp.local_depths[label][1][i] + shared_cov_arrays[label][chrom][centre_bin+1]
				except IndexError as e:
					bp.local_depths[label][2][i] = 0

	return breakpoints

if __name__ == "__main__":
	print("Breakpoint Functions")
