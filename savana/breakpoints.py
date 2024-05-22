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
from savana.clusters import cluster_breakpoints
from savana.core import PotentialBreakpoint, ConsensusBreakpoint, Cluster

"""
# developer dependencies
from memory_profiler import profile
from pympler import muppy, summary, refbrowser

# decorate functions like so:
#helper.conditionally_decorate(profile, True)
"""

#TODO: investigate/test if passing relevant read attributes to this function (rather than entire object), improves speed/mem
def get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order, single_bnd, single_bnd_max_mapq):
	""" reconstruct the breakpoints from the supplementary alignments """
	breakpoint_pairs = []
	primary_clipping = helper.get_clipping(cigar_tuples, read.is_reverse)
	primary_consumed_query = helper.sum_consumed_query(helper.trim_supplementary(read.cigarstring))
	# sort the chimeric regions by the left soft clip pos (ascending)
	chimeric_regions = sorted(chimeric_regions, key=lambda d: d['left_softclip'])
	left_index = 0
	# track the position in the read in order to extract sbnd bases
	#TODO: curr_pos_query isn't used, remove
	curr_pos_query = min(chimeric_regions[left_index]['left_softclip'], primary_clipping['left_softclip'])
	while left_index < len(chimeric_regions) and chimeric_regions[left_index]['left_softclip'] < primary_clipping['left_softclip']:
		# in the left softclip
		region_mapped = False if single_bnd and (chimeric_regions[left_index]['mapQ'] < single_bnd_max_mapq or chimeric_regions[left_index]['chrom'] not in contig_order) else True
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
				'primary': False,
				'query_pos_start': chimeric_regions[left_index]['left_softclip'],
				'query_pos_end': chimeric_regions[left_index]['left_softclip'] + chimeric_regions[left_index]['consumed_query']
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
			'primary': False,
			'query_pos_start': chimeric_regions[left_index]['left_softclip'],
			'query_pos_end': chimeric_regions[left_index]['left_softclip'] + chimeric_regions[left_index]['consumed_query']
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
			'region_mapped': False if single_bnd and read.mapq < single_bnd_max_mapq else True,
			'primary': True,
			'query_pos_start': primary_clipping['left_softclip'],
			'query_pos_end': primary_clipping['left_softclip'] + primary_consumed_query
		})
	# if still have chimeric regions, need to place in right softclip
	if left_index < len(chimeric_regions):
		# add the primary as the first edge
		breakpoint_pairs.append([{
			'chr': read.reference_name,
			'loc': read.reference_start if read.is_reverse else read.reference_end,
			'bp_notation': ("-") if read.is_reverse else ("+"),
			'region_mapped': False if single_bnd and read.mapq < single_bnd_max_mapq else True,
			'primary': True,
			'query_pos_start': primary_clipping['left_softclip'],
			'query_pos_end': primary_clipping['left_softclip'] + primary_consumed_query
		}])
		# sort remaining chimeric regions by the right soft clip (descending)
		# TODO: don't need to resort this can just use old one?
		right_chimeric = sorted(chimeric_regions[left_index:], key=lambda d: d['left_softclip'])
		right_index = 0
		while right_index < len(right_chimeric) and right_chimeric[right_index]['right_softclip'] < primary_clipping['right_softclip']:
			region_mapped = False if single_bnd and (right_chimeric[right_index]['mapQ'] < single_bnd_max_mapq or right_chimeric[right_index]['chrom'] not in contig_order) else True
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
				'primary': False,
				'query_pos_start': right_chimeric[right_index]['left_softclip'],
				'query_pos_end': right_chimeric[right_index]['left_softclip'] + right_chimeric[right_index]['consumed_query']
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
					'primary': False,
					'query_pos_start': right_chimeric[right_index]['left_softclip'],
					'query_pos_end': right_chimeric[right_index]['left_softclip'] + right_chimeric[right_index]['consumed_query']
				}])
				# increment position in the query
				curr_pos_query += right_chimeric[right_index]['consumed_query']
			right_index+=1

	# trust the query_start (based on clipping) over the end (based on consumed math from CIGAR)
	for i, (curr_start, curr_end) in enumerate(breakpoint_pairs):
		if (curr_start['query_pos_end'] > curr_end['query_pos_start']) and (curr_start['query_pos_start'] < curr_end['query_pos_start']):
			# only update if it doesn't cause the curr_start to be of 0 length
			breakpoint_pairs[i][0]['query_pos_end'] = curr_end['query_pos_start']
		elif (curr_start['query_pos_end'] > curr_end['query_pos_start']) and (curr_start['query_pos_start'] >= curr_end['query_pos_start']):
			# if it does cause a 0 length for curr_end then update curr_end's start to curr_start's end
			breakpoint_pairs[i][1]['query_pos_start'] = curr_start['query_pos_end']

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
					# next segment mapped, create single breakend from unmapped to mapped next_start
					found_mapped = True
					# get sequence for query - from current (unmapped) start to next (mapped) start
					sbnd_seq = read.query_sequence[curr_end['query_pos_start']:(next_start['query_pos_start'])]
					# only one location - next start since it's mapped
					location = [{'chr': next_start['chr'], 'loc': next_start['loc']}, {'chr': next_start['chr'], 'loc': next_start['loc']}]
					if len(sbnd_seq) > 0:
						supplementary_breakpoints.append(PotentialBreakpoint(
							location,
							"SUPPLEMENTARY",
							read.query_name,
							0, # TODO: what to set for MAPQ?
							label,
							next_start['bp_notation'], #SBND,
							sbnd_seq
						))
					else:
						print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_end["query_pos_start"]}:{next_start["query_pos_start"]} (curr_start - next_start)')
						print(breakpoint_pairs)
						print(chimeric_regions)
					# if curr_start was also mapped, need to create a second single-breakend
					if curr_start['region_mapped']:
						# get sequence for query from curr_start (mapped) end to next_start (mapped) start (same as prev unmapped end)
						sbnd_seq = read.query_sequence[curr_start['query_pos_end']:(next_start['query_pos_start'])]
						# only one location
						location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']}, {'chr': curr_start['chr'], 'loc': curr_start['loc']}]
						if len(sbnd_seq) > 0:
							supplementary_breakpoints.append(PotentialBreakpoint(
								location,
								"SUPPLEMENTARY",
								read.query_name,
								0, # TODO: what to set for MAPQ?
								label,
								curr_start['bp_notation'], #SBND
								sbnd_seq
							))
						else:
							print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos_end"]}:{next_start["query_pos_start"]} (curr_start - next_start)')
							print(breakpoint_pairs)
							print(chimeric_regions)
				elif next_end['region_mapped']:
					# segment mapped, create single breakend
					found_mapped = True
					# get sequence for query - from current (unmapped) end to next (mapped) end
					sbnd_seq = read.query_sequence[curr_start['query_pos_end']:(next_end['query_pos_end'])]
					location = [{'chr': next_end['chr'], 'loc': next_end['loc']}, {'chr': next_end['chr'], 'loc': next_end['loc']}]
					if len(sbnd_seq) > 0:
						supplementary_breakpoints.append(PotentialBreakpoint(
							location,
							"SUPPLEMENTARY",
							read.query_name,
							0, # TODO: what to set for MAPQ?
							label,
							next_end['bp_notation'], #SBND
							sbnd_seq
						))
					else:
						print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos_end"]}:{next_end["query_pos_end"]} (curr_start - next_end)')
						print(breakpoint_pairs)
						print(chimeric_regions)
					# if curr_start was also mapped, need to create a second single-breakend
					if curr_start['region_mapped']:
						# get sequence for query from curr_start (mapped) end to next_end (mapped) start (same as prev unmapped end)
						sbnd_seq = read.query_sequence[curr_start['query_pos_end']:(next_end['query_pos_start'])]
						# only one location
						location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']}, {'chr': curr_start['chr'], 'loc': curr_start['loc']}]
						if len(sbnd_seq) > 0:
							supplementary_breakpoints.append(PotentialBreakpoint(
								location,
								"SUPPLEMENTARY",
								read.query_name,
								0, # TODO: what to set for MAPQ?
								label,
								curr_start['bp_notation'], #SBND
								sbnd_seq
							))
						else:
							print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos_end"]}:{next_end["query_pos_start"]} (curr_start - next_start)')
							print(breakpoint_pairs)
							print(chimeric_regions)
				next_index += 1
			index = next_index - 1
			# got to end and didn't find map, means we should use the curr_start to last end
			if not found_mapped:
				_, last_end = breakpoint_pairs[-1]
				# get sequence for query - from curr_start end to last_end end
				sbnd_seq = read.query_sequence[curr_start['query_pos_end']:(last_end['query_pos_end'])]
				# only one location - can't be end since isn't mapped
				location = [{'chr': curr_start['chr'], 'loc': curr_start['loc']}, {'chr': curr_start['chr'], 'loc': curr_start['loc']}]
				if len(sbnd_seq) > 0:
					supplementary_breakpoints.append(PotentialBreakpoint(
						location,
						"SUPPLEMENTARY",
						read.query_name,
						0, # TODO: what to set for MAPQ?
						label,
						curr_start['bp_notation'], #SBND
						sbnd_seq
					))
				else:
					print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos_end"]}:{last_end["query_pos_end"]} (curr_end - last_end)')
					print(breakpoint_pairs)
					print(chimeric_regions)
		elif not curr_start['region_mapped']:
			# end is mapped, no need to merge/consume next pair
			sbnd_seq = read.query_sequence[curr_start['query_pos_start']:(curr_end['query_pos_start'])]
			# only one location - the mapped end
			location = [{'chr': curr_end['chr'], 'loc': curr_end['loc']}, {'chr': curr_end['chr'], 'loc': curr_end['loc']}]
			if len(sbnd_seq) > 0:
				supplementary_breakpoints.append(PotentialBreakpoint(
					location,
					"SUPPLEMENTARY",
					read.query_name,
					0, # TODO: what to set for MAPQ?
					label,
					curr_end['bp_notation'], #SBND
					sbnd_seq
				))
			else:
				print(f'sbnd_seq empty! {read.reference_name}:{read.reference_start} {read.query_name} {len(read.query_sequence)} {curr_start["query_pos_start"]}:{curr_end["query_pos_start"]} (curr_start - curr_start)')
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
						"".join((curr_start['bp_notation'], curr_end['bp_notation']))
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
			else:
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

def get_potential_breakpoints(aln_filename, is_cram, ref, length, mapq, label, contig_order, contig, start, end, coverage_binsize, contig_coverage_array, single_bnd=False, single_bnd_min_length=None, single_bnd_max_mapq=None):
	""" iterate through alignment file, tracking potential breakpoints and saving relevant reads to fastq """
	potential_breakpoints = {}
	aln_file = pysam.AlignmentFile(aln_filename, "rc", reference_filename=ref) if is_cram else pysam.AlignmentFile(aln_filename, "rb")
	# adjust the thresholds depending on sample source
	args_length = max((length - floor(length/5)), 0) if label == 'normal' else length
	mapq = min((mapq - ceil(mapq/2)), 0) if label == 'normal' else mapq
	for read in aln_file.fetch(contig, start, end):
		if read.is_secondary or read.is_unmapped:
			continue
		try:
			# record start/end in read incrementer
			contig_coverage_array[floor((read.reference_start-1)/coverage_binsize)]+=1
			contig_coverage_array[floor((read.reference_end-1)/coverage_binsize)]-=1
		except IndexError as _:
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
			chimeric_breakpoints = get_supplementary_breakpoints(read, cigar_tuples, chimeric_regions, label, contig_order, single_bnd, single_bnd_max_mapq)
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
					potential_breakpoints.setdefault(curr_chrom,[]).append(
						PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-")
						)
					prev_deletion = []
				# record the insertion
				inserted_sequence = read.query_sequence[curr_pos['query']:(curr_pos['query']+length)]
				potential_breakpoints.setdefault(curr_chrom,[]).append(
					PotentialBreakpoint(location, "CIGAR", read.query_name, read.mapping_quality, label, "<INS>", inserted_sequence)
					)
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
			elif sam_flag == helper.samflag_desc_to_number["BAM_CSOFT_CLIP"] and single_bnd and not chimeric_regions and not read.is_supplementary and length >= single_bnd_min_length:
				# start and end location for single breakend are the same (for clustering purposes)
				location = [
					{'chr': curr_chrom, 'loc': curr_pos['reference']},
					{'chr': curr_chrom, 'loc': curr_pos['reference']}
				]
				# record the softclipped bases
				softclipped_sequence = read.query_sequence[curr_pos['query']:(curr_pos['query']+length)]
				breakpoint_notation = "-" if curr_pos['query'] == 0 else "+"
				if len(softclipped_sequence) > 0:
					potential_breakpoints.setdefault(curr_chrom,[]).append(
						PotentialBreakpoint(location, "SOFTCLIP", read.query_name, read.mapping_quality, label, breakpoint_notation, softclipped_sequence)
						)
				else:
					print(f'softclipped sequence empty! (from {curr_pos["query"]}-{curr_pos["query"]+length})')
			elif prev_deletion and length > args_length:
				# record and clear the previously tracked deletion
				potential_breakpoints.setdefault(curr_chrom,[]).append(
					PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-")
					)
				prev_deletion = []
			# increment values
			curr_pos['cigar'] += length
			if helper.consumes_query[sam_flag]:
				curr_pos['query'] += length
			if helper.consumes_reference[sam_flag]:
				curr_pos['reference'] += length
		if prev_deletion:
			# if reached end of string and no chance to expand deletion, add it
			potential_breakpoints.setdefault(curr_chrom,[]).append(
				PotentialBreakpoint(prev_deletion, "CIGAR", read.query_name, read.mapping_quality, label, "+-")
				)
	aln_file.close()
	del aln_file

	return potential_breakpoints

def cluster_by_insert_length(insertion_like_breakpoints, fraction):
	""" separate insertion breakpoints by their insert lengths - add single-breakends to largest insertion cluster """
	# sort by insert length
	cigar_insertions, sbnds = [], []
	for bp in insertion_like_breakpoints:
		if bp.breakpoint_notation == "<INS>":
			cigar_insertions.append(bp)
		else:
			sbnds.append(bp)
	cigar_insertions.sort(key=lambda bp: len(bp.inserted_sequence))
	stack = []
	if cigar_insertions:
		for bp in cigar_insertions:
			if len(stack) == 0:
				new_cluster = [bp]
				stack.append(new_cluster)
			elif ceil(len(stack[-1][-1].inserted_sequence)*(2-fraction)) >= len(bp.inserted_sequence):
				stack[-1].append(bp)
			else:
				new_cluster = [bp]
				stack.append(new_cluster)
		# now add the single-breakends to the largest cluster on the stack (last one)
		stack[-1].extend(sbnds)
	else:
		stack.append(sbnds)

	return stack

def call_breakpoints(clusters, end_buffer, min_length, min_support, chrom):
	""" identify consensus breakpoints from list of clusters """
	# N.B. all breakpoints in a cluster must be from same start chromosome!
	#TODO: refactor. this works but it's ugly and overcomplicated
	final_breakpoints = []
	for cluster in clusters:
		# separate breakpoints by their end chromosome (ignoring type)
		breakpoints_by_end_chrom = {}
		for bp in cluster.breakpoints:
			breakpoints_by_end_chrom.setdefault(bp.end_chr, []).append(bp)
		for end_chrom, end_chrom_breakpoints in breakpoints_by_end_chrom.items():
			breakpoints_for_end_chrom = []
			if len(end_chrom_breakpoints) >= min_support:
				# TODO: this might not be necessary if we decide to scrap in model
				read_counts = {"tumour": 0, "normal": 0}
				for bp in end_chrom_breakpoints:
					read_counts[bp.label] += 1
				# first deal with insertion/single-breakends
				insertion_like_breakpoints = [b for b in end_chrom_breakpoints if b.breakpoint_notation in ["<INS>", "+", "-"]]
				if len(insertion_like_breakpoints) >= min_support:
					insertion_like_clusters = cluster_by_insert_length(insertion_like_breakpoints, 0.75)
					for ins_cluster in insertion_like_clusters:
						ins_label_counts = count_num_labels(ins_cluster)
						if max(len(count) for count in ins_label_counts.values()) >= min_support:
							num_insertions = sum(b.breakpoint_notation == "<INS>" for b in ins_cluster)
							if num_insertions >= 2:
								# require at least two cigar insertions of similar lengths to call ins
								# create new "originating cluster" with only used breakpoints
								# (ensures statistics are calculated correctly)
								event_info = {'starts':[], 'inserts': [], 'sources': {}}
								start_cluster = None
								for bp in ins_cluster:
									event_info['starts'].append(bp.start_loc)
									event_info['inserts'].append(bp.inserted_sequence)
									event_info['sources'].setdefault(bp.source, True)
									if not start_cluster:
										start_cluster = Cluster(bp)
									else:
										start_cluster.add(bp)
								median_start = median(event_info['starts'])
								consensus_source = "/".join(sorted(event_info['sources'].keys()))
								breakpoints_for_end_chrom.append(ConsensusBreakpoint(
									[{'chr': cluster.chr, 'loc': median_start}, {'chr': cluster.chr, 'loc': median_start}],
									consensus_source, start_cluster, None, ins_label_counts, "<INS>", read_counts, event_info['inserts']))
							else:
								# separate here by bp_notation
								sorted_breakpoints = {}
								for bp in ins_cluster:
									sorted_breakpoints.setdefault(bp.breakpoint_notation, []).append(bp)
								for notation_type, breakpoints in sorted_breakpoints.items():
									sbnd_label_counts = count_num_labels(breakpoints)
									event_info = {'starts':[], 'inserts': [], 'sources': {}}
									start_cluster = None
									for bp in breakpoints:
										event_info['starts'].append(bp.start_loc)
										event_info['inserts'].append(bp.inserted_sequence)
										event_info['sources'].setdefault(bp.source, True)
										if not start_cluster:
											start_cluster = Cluster(bp)
										else:
											start_cluster.add(bp)
									median_start = median(event_info['starts'])
									consensus_source = "/".join(sorted(event_info['sources'].keys()))
									if max(len(count) for count in sbnd_label_counts.values()) >= min_support:
										breakpoints_for_end_chrom.append(ConsensusBreakpoint(
											[{'chr': cluster.chr, 'loc': median_start}, {'chr': cluster.chr, 'loc': median_start}],
											consensus_source, start_cluster, None, sbnd_label_counts, notation_type, read_counts, event_info['inserts']))
				# now deal with other breakpoint notation types
				# reverse start/end in order to re-cluster
				flipped_breakpoints = [reversed(b) for b in end_chrom_breakpoints if b.breakpoint_notation not in ["<INS>", "+", "-"]]
				# sort by new start
				flipped_breakpoints.sort()
				# cheeky reclustering with reversed breakpoints
				_, cluster_stack = cluster_breakpoints(end_chrom, flipped_breakpoints, end_buffer)
				for end_chrom_cluster in cluster_stack:
					# sort breakpoints by their breakpoint_notation
					sorted_breakpoints = {}
					for bp in end_chrom_cluster.breakpoints:
						sorted_breakpoints.setdefault(bp.breakpoint_notation, []).append(bp)
					# go through the breakpoints by notation type
					for bp_type, breakpoints in sorted_breakpoints.items():
						# create a consensus breakpoint for each end cluster that has enough supporting reads
						label_counts = count_num_labels(breakpoints)
						if max([len(count) for count in label_counts.values()]) >= min_support:
							# un-flip the breakpoints
							unflipped_breakpoints = [reversed(b) for b in breakpoints]
							# create a new "originating cluster" with only the used breakpoints
							# (ensures statistics are calculated correctly)
							start_cluster = None
							source = {}
							for bp in unflipped_breakpoints:
								if not start_cluster:
									start_cluster = Cluster(bp)
								else:
									start_cluster.add(bp)
								source.setdefault(bp.source, True)
							median_start = median([bp.start_loc for bp in unflipped_breakpoints])
							median_end = median([bp.end_loc for bp in unflipped_breakpoints])
							consensus_source = "/".join(sorted(source.keys()))
							new_breakpoint = ConsensusBreakpoint(
								[{'chr': start_cluster.chr, 'loc': median_start}, {'chr': end_chrom_cluster.chr, 'loc': median_end}],
								consensus_source, start_cluster, end_chrom_cluster, label_counts, bp_type, read_counts)
							if new_breakpoint.sv_length >= min_length or (new_breakpoint.start_chr != new_breakpoint.end_chr and new_breakpoint.sv_length == 0):
								breakpoints_for_end_chrom.append(new_breakpoint)

			if breakpoints_for_end_chrom:
				# only report single-breakends if no other breakpoint type reported
				counts = {}
				for new_breakpoint in breakpoints_for_end_chrom:
					counts[new_breakpoint.breakpoint_notation] = counts.get(new_breakpoint.breakpoint_notation, 0) + 1
				bp_types = counts.keys()
				if "-" in bp_types or "+" in bp_types:
					major_types = [t for t in bp_types if t not in ("+", "-")]
					if major_types:
						major_types_bps = [b for b in breakpoints_for_end_chrom if b.breakpoint_notation not in ("+","-")]
						final_breakpoints.extend(major_types_bps)
					else:
						# no other types present, add them
						final_breakpoints.extend(breakpoints_for_end_chrom)
				else:
					final_breakpoints.extend(breakpoints_for_end_chrom)

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
				if chrom not in shared_cov_arrays[label]:
					print(f'WARNING: contig {chrom} not in shared_cov_array!')
					continue
				# calculate centre bin (zero-based array, subtract 1)
				centre_bin = floor((loc)/coverage_binsize)
				# ensure is bound by 0 and last element of contig coverage array
				centre_bin = min(max(centre_bin, 0), (len(shared_cov_arrays[label][chrom]) - 1))
				# before
				if centre_bin == 0:
					# no coverage since out of range
					bp.local_depths[label][0][i] = 0
				else:
					bp.local_depths[label][0][i] = np.sum(shared_cov_arrays[label][chrom][0:centre_bin-1])
				# at
				bp.local_depths[label][1][i] = bp.local_depths[label][0][i] + shared_cov_arrays[label][chrom][centre_bin]
				# after
				if centre_bin == (len(shared_cov_arrays[label][chrom]) - 1):
					# no coverage since out of range for contig
					bp.local_depths[label][2][i] = 0
				else:
					bp.local_depths[label][2][i] = bp.local_depths[label][1][i] + shared_cov_arrays[label][chrom][centre_bin+1]

		# now calculate the allele fractions
		for label, [_, dp_at, _] in bp.local_depths.items():
			af = [None, None]
			for i in [0,1]:
				af[i] = round(bp.support[label]/dp_at[i], 3) if dp_at[i] != 0 else 0.0
				# due to sloping edges this can sometimes be inflated
				af[i] = 1.0 if af[i] > 1.0 else af[i]
			bp.allele_fractions[label] = af
	return breakpoints

if __name__ == "__main__":
	print("Breakpoint Functions")
