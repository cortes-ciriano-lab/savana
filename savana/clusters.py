"""
Module containing functions related to the Cluster class for SAVANA
Created: 06/09/2022
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import os
import json
import subprocess
from statistics import median, mean

import pysam
import pybedtools

from savana.core import Cluster

def cluster_breakpoints(chrom, breakpoints, buffer, ins_buffer=None):
	""" given a list of Breakpoints (starting on same chrom) cluster them on location and type """
	stack = []
	breakpoints.sort()
	for bp in breakpoints:
		bp_notation_type = str(bp.breakpoint_notation)
		bp_buffer = ins_buffer if (ins_buffer and bp_notation_type == "<INS>") else buffer
		if len(stack) == 0:
			# put a new cluster on the stack
			new_cluster = Cluster(bp)
			stack.append(new_cluster)
		elif not stack[-1].overlaps(bp, bp_buffer):
			# put a new cluster on the stack
			new_cluster = Cluster(bp)
			stack.append(new_cluster)
		else:
			stack[-1].add(bp)

	# require two supporting reads per cluster
	filtered_stack = [cluster for cluster in stack if len(cluster.supporting_reads) >= 2]

	return chrom, filtered_stack

def write_cluster_bed(clusters, outdir):
	""" store clusters in bed file"""
	cluster_file = os.path.join(outdir, 'cluster.bed')
	cluster_file_compressed = f'{cluster_file}.gz'
	cluster_bed = ''
	for clusters_sv_type in clusters.values():
		for cluster in clusters_sv_type:
			if cluster.start <= cluster.end:
				cluster_bed+="\t".join([cluster.chr, str(cluster.start), str(cluster.end)])
			else:
				cluster_bed+="\t".join([cluster.chr, str(cluster.end), str(cluster.start)])
			cluster_bed+="\n"
	sorted_bed = pybedtools.BedTool(cluster_bed, from_string=True).sort()
	sorted_bed.saveas(cluster_file)
	pysam.tabix_compress(cluster_file, cluster_file_compressed)
	pysam.tabix_index(cluster_file_compressed, preset='bed', keep_original=True)

def calculate_cluster_stats(clusters, outdir):
	""" compute and output statistical information about clusters """
	num_breakpoints = []
	breakpoint_counter = {}
	seen_query_names = {}
	qualities = []
	cluster_sizes = []
	num_clusters = 0
	for clusters_sv_type in clusters.values():
		for cluster in clusters_sv_type:
			num_clusters += 1
			# to prevent double-counting qualities, track which have been seen
			for b in cluster.breakpoints:
				if b.read_name not in seen_query_names:
					qualities.append(b.mapq)
					seen_query_names[b.read_name] = True
			cluster_sizes.append(abs(cluster.start - cluster.end))
			n_bp = len(cluster.breakpoints)
			num_breakpoints.append(n_bp)
			if n_bp not in breakpoint_counter:
				breakpoint_counter[n_bp] = 1
			else:
				breakpoint_counter[n_bp] += 1

	stats = {
		"Number of Clusters": num_clusters,
		"Mean MAPQ": round(mean(qualities), 2),
		"Median MAPQ": median(qualities),
		"Mean Number of Breakpoints": round(mean(num_breakpoints),2),
		"Median Number of Breakpoints": median(num_breakpoints),
		"Max Number of Breakpoints": max(num_breakpoints),
		"Mean Cluster Window Size": round(mean(cluster_sizes),2),
		"Median Cluster Window Size": median(cluster_sizes)
	}

	# print outputs and write to outdir
	f = open(os.path.join(outdir, 'clusters.stats'), "w+")
	f.write(",".join(stats.keys())+"\n")
	f.write(",".join([str(v) for v in stats.values()])+"\n")
	f.close()

if __name__ == "__main__":
	print("Clustering Functions")
