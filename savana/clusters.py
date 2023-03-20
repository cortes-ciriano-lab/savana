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

def cluster_breakpoints(breakpoints, buffer, ins_buffer):
	""" given a list of Breakpoints (starting on same chrom) cluster them on location and type """
	cluster_stacks = {
		"+-": [],
		"++": [],
		"-+": [],
		"--": [],
		"<INS>": []
	}
	breakpoints.sort()
	for bp in breakpoints:
		bp_notation_type = str(bp.breakpoint_notation)
		if len(cluster_stacks[bp_notation_type]) == 0:
			# put a new cluster onto the sv stack
			new_cluster = Cluster(bp)
			cluster_stacks[bp_notation_type].append(new_cluster)
		elif bp_notation_type == "<INS>":
			if not cluster_stacks[bp_notation_type][-1].overlaps(bp, ins_buffer):
				# put a new cluster onto the sv stack
				new_cluster = Cluster(bp)
				cluster_stacks[bp_notation_type].append(new_cluster)
			else:
				# add to cluster on top of stack
				cluster_stacks[bp_notation_type][-1].add(bp)
		elif bp_notation_type != "<INS>" and not cluster_stacks[bp_notation_type][-1].overlaps(bp, buffer):
			# put a new cluster onto the sv stack
			new_cluster = Cluster(bp)
			cluster_stacks[bp_notation_type].append(new_cluster)
		else:
			# add to cluster on top of stack
			cluster_stacks[bp_notation_type][-1].add(bp)
	for bp_notation_type, stack in cluster_stacks.items():
		# can't cluster with only one read - require two
		filtered_cluster_stacks = [c for c in stack if len(c.supporting_reads) >= 2]
		cluster_stacks[bp_notation_type] = filtered_cluster_stacks
	return cluster_stacks

def output_clusters(refined_clusters, outdir):
	""" output the json files of evidence """
	for cluster in refined_clusters:
		cluster_id = str(cluster.uid)
		cluster_outdir = os.path.join(outdir, 'clusters', cluster_id)
		if not os.path.exists(cluster_outdir):
			os.makedirs(cluster_outdir)
		output_json = open(os.path.join(cluster_outdir, f'{cluster_id}.json'), 'w')
		json.dump(cluster.as_dict(), output_json, sort_keys=False, indent=2)
		output_json.close()

def wrap_subprocess(command, outfile=None, wait=False):
	""" error handling for subprocess commands """
	if wait and outfile:
		try:
			with open(outfile, 'w') as out:
				p = subprocess.Popen(command, stdout=out, stderr=subprocess.DEVNULL)
				p.wait()
		except Exception as e:
			print(f'Command {command} failed with error (code {e.returncode}): {e.output}')
	if outfile:
		try:
			with open(outfile, 'w') as out:
				subprocess.run(command, stdout=out, stderr=subprocess.DEVNULL)
		except Exception as e:
			print(f'Command {command} failed with error (code {e.returncode}): {e.output}')

def write_cluster_bed(clusters, outdir):
	""" store clusters in bed file"""
	cluster_file = os.path.join(outdir, 'cluster.bed')
	cluster_file_compressed = f'{cluster_file}.gz'
	cluster_bed = ''
	for clusters_sv_type in clusters.values():
		for cluster in clusters_sv_type:
			cluster_id = str(cluster.uid)
			cluster_bed+="\t".join([cluster.chr, str(cluster.start), str(cluster.end), cluster_id])
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
