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

import savana.helper as helper
from savana.breakpoints import *
from savana.clusters import *

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

def spawn_processes(args, bam_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using multiprocessing with {args.threads} threads\n')
	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints_results = pool_get_potential_breakpoints(bam_files, args)
	if args.debug:
		helper.time_function("Identified potential breakpoints", checkpoints, time_str)
	# collect results per chrom
	chrom_potential_breakpoints = {}
	for result in potential_breakpoints_results:
		for chrom, potential_breakpoints in result.items():
			chrom_potential_breakpoints.setdefault(chrom,[]).extend(potential_breakpoints)

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clusters = pool_cluster_breakpoints(args, chrom_potential_breakpoints)
	if args.debug:
		helper.time_function("Clustered potential breakpoints", checkpoints, time_str)

	if args.debug:
		# 2.1) OUTPUT CLUSTERS
		for bp_type in ["+-", "++", "-+", "--", "<INS>"]:
			pool_output_clusters(args, clusters[bp_type], outdir)
		helper.time_function("Output originating clusters", checkpoints, time_str)

	# 3) CALL BREAKPOINTS
	breakpoints = call_breakpoints(clusters, args.buffer)
	if args.debug:
		helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	# 3.1) OUTPUT BREAKPOINTS
	bedpe_string = ''
	vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, breakpoints[0])
	read_support_string = ''
	ref_fasta = pysam.FastaFile(args.ref)
	for count, bp in enumerate(breakpoints):
		bedpe_string += bp.as_bedpe(count)
		vcf_string += bp.as_vcf(ref_fasta)
		read_support_string += bp.as_read_support(count)
	with open(os.path.join(outdir, f'{args.sample}.sv_breakpoints.vcf'), 'w') as output:
		output.write(vcf_string)
	with open(os.path.join(outdir, f'{args.sample}.sv_breakpoints.bedpe'), 'w') as output:
		output.write(bedpe_string)
	with open(os.path.join(outdir, f'{args.sample}.sv_breakpoints_read_support.tsv'), 'w') as output:
		output.write(read_support_string)
	if args.debug:
		helper.time_function("Output consensus breakpoints", checkpoints, time_str)

	return clusters, breakpoints, checkpoints, time_str

if __name__ == "__main__":
	print("Functions to run SAVANA")
