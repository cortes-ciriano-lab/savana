"""
SAVANA strucural variant caller for long-read data - run
Created: 17/02/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import os
import gc

from math import ceil, floor
from multiprocessing import Pool
from multiprocessing import Array, Pipe, Process

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
import objgraph
from pympler import muppy, summary, refbrowser
"""

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

def execute_get_potential_breakpoint_task(task_arg_dict, contig_coverage_array, task_tracker, conn):
	""" submit the task arguments to the function and send the results through the pipe conn """
	potential_breakpoints = get_potential_breakpoints(
		task_arg_dict['aln_file'],
		task_arg_dict['is_cram'],
		task_arg_dict['ref'],
		task_arg_dict['length'],
		task_arg_dict['mapq'],
		task_arg_dict['label'],
		task_arg_dict['contig_order'],
		task_arg_dict['contig'],
		task_arg_dict['start_pos'],
		task_arg_dict['end_pos'],
		contig_coverage_array
		)
	task_tracker[task_arg_dict['task_id']] = 1
	conn.send(potential_breakpoints)
	# cleanup
	del task_arg_dict
	gc.collect()
	conn.close()

def generate_get_potential_breakpoint_tasks(aln_files, args):
	""" generate get_potential_breakpoint tasks by chunking the genome """
	tasks = []
	task_id_counter = 0
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	for label, aln_file in aln_files.items():
		contigs = helper.get_contigs(args.contigs, args.ref_index)
		if not args.is_cram:
			# if nothing mapped, don't consider contig
			[contigs.remove(c.contig) for c in aln_file.get_index_statistics() if not c.mapped and (c.contig in contigs)]
		for contig, contig_length in contig_lengths.items():
			if contig not in contigs:
				continue
			if contig_length > args.chunksize:
				# split the chrom into parts
				num_intervals = floor(contig_length/args.chunksize) + 1
				start_pos = 0
				for _ in range(0, num_intervals):
					end_pos = start_pos + args.chunksize
					end_pos = contig_length if end_pos > contig_length else end_pos # don't extend past end
					tasks.append({
						'aln_file': aln_file.filename,
						'label': label,
						'is_cram': args.is_cram,
						'ref': args.ref,
						'length': args.length,
						'mapq': args.mapq,
						'contig_order': contigs,
						'contig': contig,
						'start_pos': start_pos,
						'end_pos': end_pos,
						'task_id': task_id_counter
					})
					start_pos = end_pos + 1
					task_id_counter += 1
			else:
				tasks.append({
					'aln_file': aln_file.filename,
					'label': label,
					'is_cram': args.is_cram,
					'ref': args.ref,
					'length': args.length,
					'mapq': args.mapq,
					'contig_order': contigs,
					'contig': contig,
					'start_pos': 0,
					'end_pos': contig_length,
					'task_id': task_id_counter
				})
				task_id_counter += 1
	return tasks

def generate_coverage_arrays(aln_files, args):
	""" generate empty numpy arrays """
	from multiprocessing import sharedctypes
	coverage_arrays = {}
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	for label, _ in aln_files.items():
		for contig, contig_length in contig_lengths.items():
			if contig not in contigs_to_consider:
				continue
			coverage_arrays.setdefault(label, {})[contig] = sharedctypes.RawArray('h', contig_length)
	return coverage_arrays

def execute_call_breakpoints(task_arg_dict, task_tracker, conn):
	""" submit task arguments to the function and send results through pipe """
	breakpoints, _, contig = call_breakpoints(
		task_arg_dict['clusters'],
		task_arg_dict['buffer'],
		task_arg_dict['length'],
		task_arg_dict['depth'],
		task_arg_dict['contig']
	)
	task_tracker[task_arg_dict['task_id']] = 1
	conn.send((breakpoints, contig))
	# cleanup
	del task_arg_dict
	gc.collect()
	conn.close()

def execute_cluster_breakpoints(task_arg_dict, task_tracker, conn):
	""" submit the task arguments to the function and send results through pipe """
	contig, clustered_breakpoints = cluster_breakpoints(
		task_arg_dict['contig'],
		task_arg_dict['breakpoints'],
		task_arg_dict['buffer'],
		task_arg_dict['ins_buffer'],
	)
	task_tracker[task_arg_dict['task_id']] = 1
	conn.send((contig, clustered_breakpoints))
	# cleanup
	del task_arg_dict
	gc.collect()
	conn.close()

def generate_cluster_breakpoints_tasks(potential_breakpoints, args):
	""" """
	tasks = []
	task_id_counter = 0
	for contig, breakpoints in potential_breakpoints.items():
		tasks.append({
			'contig': contig,
			'breakpoints': breakpoints,
			'buffer': args.buffer,
			'ins_buffer': args.insertion_buffer,
			'task_id': task_id_counter
		})
		task_id_counter += 1
	return tasks

def generate_call_breakpoint_tasks(clustered_breakpoints, args):
	""""""
	tasks = []
	task_id_counter = 0
	for contig, clusters in clustered_breakpoints.items():
		tasks.append({
			'contig': contig,
			'clusters': clusters,
			'buffer': args.buffer,
			'length': args.length,
			'depth': args.depth,
			'task_id': task_id_counter
		})
		task_id_counter += 1
	return tasks

def run_call_breakpoints(clustered_breakpoints, args):
	""" generate task list for calling breakpoints and coordinate the execution """
	tasks = generate_call_breakpoint_tasks(clustered_breakpoints, args)
	task_tracker = Array('h', [0]*len(tasks))
	print(f' > {len(tasks)} calling tasks')

	results = []
	pipes = [None] * min(args.threads, len(tasks))
	processes = [None] * min(args.threads, len(tasks))
	while tasks or processes.count(None) != len(processes):
		for i, process in enumerate(processes):
			# check if free process and remaining tasks
			if not process and tasks:
				# pop task off task list
				task_args = tasks.pop(0)
				# create new pipe
				pipes[i] = Pipe(duplex=False)
				# create a new process and replace it (match pipe w process)
				proc = Process(target=execute_call_breakpoints, args=(task_args, task_tracker, pipes[i][1]))
				processes[i] = (proc, task_args['task_id'])
				proc.start()
		for i, process in enumerate(processes):
			if process:
				proc, assigned_task_id = process
				result_conn = pipes[i][0]
				if proc and not proc.is_alive():
					if result_conn.poll():
						result = result_conn.recv()
						proc.join()
						result_conn.close()
						results.append(result)
					del process
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see task is done
					if result_conn.poll(10):
						result = result_conn.recv()
						proc.join(timeout=1) # these hang
						result_conn.close()
						results.append(result)
					else:
						print(f'Timeout polling process {i} (task {assigned_task_id})')
					del process
					processes[i] = None
					pipes[i] = None
	del task_tracker
	called_breakpoints = {}
	for result in results:
		breakpoints, contig = result
		called_breakpoints.setdefault(contig, []).append(breakpoints)

	return called_breakpoints

def run_cluster_breakpoints(potential_breakpoints, args):
	""" generate the task list and coordinate the task execution """
	tasks = generate_cluster_breakpoints_tasks(potential_breakpoints, args)
	task_tracker = Array('h', [0]*len(tasks))
	print(f' > {len(tasks)} clustering tasks')

	results = []
	pipes = [None] * min(args.threads, len(tasks))
	processes = [None] * min(args.threads, len(tasks))
	while tasks or processes.count(None) != len(processes):
		for i, process in enumerate(processes):
			# check if free process and remaining tasks
			if not process and tasks:
				# pop task off the tasks list
				task_args = tasks.pop(0)
				# create new pipe
				pipes[i] = Pipe(duplex=False)
				# create a new process and replace it (matching pipe with its process)
				proc = Process(target=execute_cluster_breakpoints, args=(task_args, task_tracker, pipes[i][1]))
				processes[i] = (proc, task_args['task_id'])
				proc.start()
		for i, process in enumerate(processes):
			if process:
				proc, assigned_task_id = process
				result_conn = pipes[i][0]
				if proc and not proc.is_alive():
					if result_conn.poll():
						result = result_conn.recv()
						proc.join()
						result_conn.close()
						results.append(result)
					del process
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see task is done
					if result_conn.poll(10):
						result = result_conn.recv()
						proc.join(timeout=1) # these processes hang
						result_conn.close()
						results.append(result)
					else:
						print(f'Timeout polling process {i} (task {assigned_task_id})')
					del process
					processes[i] = None
					pipes[i] = None
	del task_tracker
	clustered_breakpoints = {}
	for result in results:
		contig, breakpoints = result
		clustered_breakpoints[contig] = breakpoints

	return clustered_breakpoints

def run_get_potential_breakpoints(aln_files, args):
	""" generate the task list and coordinate the task execution by processes  """
	tasks = generate_get_potential_breakpoint_tasks(aln_files, args)
	task_tracker = Array('h', [0]*len(tasks))
	shared_cov_arrays = generate_coverage_arrays(aln_files, args)

	results = []
	pipes = [None] * args.threads
	processes = [None] * args.threads
	regions_running = [None] * args.threads
	print(f'Total of {len(tasks)} tasks and {len(processes)} processes') if args.debug else None
	while tasks or processes.count(None) != len(processes):
		for i, process in enumerate(processes):
			# check if free process and remaining tasks
			if not process and tasks:
				# pop task off the tasks list
				task_args = tasks.pop(0)
				# retrieve reference to the relevant contig's coverage array
				contig_coverage_array = shared_cov_arrays[task_args['label']][task_args['contig']]
				# create a new pipe
				pipes[i] = Pipe(duplex=False)
				# create a new process and replace it (matching the pipe with its process)
				proc = Process(target=execute_get_potential_breakpoint_task, args=(task_args, contig_coverage_array, task_tracker, pipes[i][1]))
				processes[i] = (proc, task_args['task_id'])
				regions_running[i] = (task_args['label'], task_args['contig'], task_args['start_pos'], task_args['end_pos'])
				proc.start()
		# Check for completed processes/tasks
		for i, process in enumerate(processes):
			if process:
				proc, assigned_task_id = process
				result_conn = pipes[i][0]
				if proc and not proc.is_alive():
					if result_conn.poll():
						result = result_conn.recv()
						proc.join()
						result_conn.close()
						results.append(result)
					del process
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see that the task is done
					if result_conn.poll(10):
						result = result_conn.recv()
						proc.join(timeout=1) # timeout because these processes hang
						result_conn.close()
						results.append(result)
					else:
						print(f'Timeout polling process {i} (task {assigned_task_id})')
						print(f' > label: {regions_running[i][0]}')
						print(f' > region: {regions_running[i][1]}:{regions_running[i][2]}-{regions_running[i][3]}')
					del process
					processes[i] = None
					pipes[i] = None
	# TODO: handle case where some tasks don't complete
	del task_tracker
	contig_potential_breakpoints = {}
	for result in results:
		for contig, potential_breakpoints in result.items():
			if potential_breakpoints:
				contig_potential_breakpoints.setdefault(contig, []).extend(potential_breakpoints)

	return contig_potential_breakpoints, shared_cov_arrays

def spawn_processes(args, aln_files, checkpoints, time_str, outdir):
	""" run main algorithm steps in parallel processes """
	print(f'Using {args.threads} thread(s)\n')

	# 1) GET POTENTIAL BREAKPOINTS
	potential_breakpoints, shared_cov_arrays = run_get_potential_breakpoints(aln_files, args)
	helper.time_function("Identified potential breakpoints", checkpoints, time_str)

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clustered_breakpoints = run_cluster_breakpoints(potential_breakpoints, args)
	helper.time_function("Clustered potential breakpoints", checkpoints, time_str)

	# 3) CALL BREAKPOINTS FROM CLUSTERS
	called_breakpoints = run_call_breakpoints(clustered_breakpoints, args)
	helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	# cleanup
	del potential_breakpoints
	del clustered_breakpoints
	gc.collect()

	#print(sum(shared_cov_arrays['tumour']['chr10'][0:53717136]))
	#print(sum(coverage_arrays['tumour']['chr1'][0:224595093]))

	# 4) COMPUTE LOCAL DEPTH
	#multithreading_compute_depth(args.threads, breakpoint_dict_chrom, contig_coverages_merged, args.debug)
	#helper.time_function("Computed local depth for breakpoints", checkpoints, time_str)

	return checkpoints, time_str

	# 2) CLUSTER POTENTIAL BREAKPOINTS
	clusters = pool_cluster_breakpoints(args.threads, args.buffer, args.insertion_buffer, potential_breakpoints)
	helper.time_function("Clustered potential breakpoints", checkpoints, time_str)

	# 3) CALL BREAKPOINTS FROM CLUSTERS
	breakpoint_dict_chrom, pruned_clusters = pool_call_breakpoints(args.threads, args.buffer, args.length, args.depth, clusters, args.debug)
	helper.time_function("Called consensus breakpoints", checkpoints, time_str)

	if args.debug:
		total_breakpoints = 0
		for _, b in breakpoint_dict_chrom.items():
			total_breakpoints+=len(b)
		print(f' > Total number of breakpoints post-calling: {total_breakpoints}')

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

	# 1) GET POTENTIAL BREAKPOINTS
	#potential_breakpoints_results = pool_get_potential_breakpoints(aln_files, args)
	#potential_breakpoints_results = multiprocessing_get_potential_breakpoints(aln_files, args)
	#potential_breakpoints_results = multithreading_get_potential_breakpoints(aln_files, args)
	potential_breakpoints_results = single_thread_get_potential_breakpoints(aln_files, args)

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
	helper.time_function("Identified potential breakpoints", checkpoints, time_str)

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
