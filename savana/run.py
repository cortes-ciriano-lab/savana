"""
SAVANA strucural variant caller for long-read data - run
Created: 17/02/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import sys
import os
import gc

from math import ceil, floor
from multiprocessing import Array, Pipe, Process, Manager, sharedctypes

import numpy as np
import pysam
import pysam.bcftools as bcftools

import savana.helper as helper
from savana.breakpoints import get_potential_breakpoints, call_breakpoints, compute_depth
from savana.clusters import cluster_breakpoints

"""
# developer dependencies
from memory_profiler import profile
import objgraph
from pympler import muppy, summary, refbrowser
"""

def execute_annotate(task_arg_dict, contig_coverage_array, task_tracker, conn):
	""" submit task arguments to the annotation function and send results through pipe """
	breakpoints = compute_depth(
		task_arg_dict['breakpoints'],
		contig_coverage_array,
		task_arg_dict['coverage_binsize']
	)
	task_tracker[task_arg_dict['task_id']] = 1
	conn.send((breakpoints))
	# cleanup
	del task_arg_dict
	gc.collect()
	conn.close()

def execute_call_breakpoints(task_arg_dict, task_tracker, conn):
	""" submit task arguments to the function and send results through pipe """
	breakpoints, contig = call_breakpoints(
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

def execute_get_potential_breakpoint_task(task_arg_dict, task_tracker, conn):
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
		task_arg_dict['coverage_binsize'],
		task_arg_dict['contig_coverage_array'],
		task_arg_dict['single_bnd'],
		task_arg_dict['single_bnd_min_length'],
		task_arg_dict['single_bnd_max_mapq']
		)
	task_tracker[task_arg_dict['task_id']] = 1
	conn.send(potential_breakpoints)
	# cleanup
	del task_arg_dict
	gc.collect()
	conn.close()

def generate_annotate_depth_tasks(called_breakpoints, args):
	""" generate task list by splitting work into approx. evenly-sized chunks """
	tasks = []
	task_id_counter = 0
	total_breakpoints = sum(len(bps) for bps in called_breakpoints.values())
	split = max(floor(total_breakpoints/(args.threads*10)), 1)
	spare_breakpoints = []
	for _, breakpoints in called_breakpoints.items():
		if len(breakpoints) > split:
			n_chunks = ceil(len(breakpoints)/split)
			for chunk in list(map(lambda x: breakpoints[x * split:x * split + split], list(range(n_chunks)))):
				tasks.append({
					'breakpoints': chunk,
					'coverage_binsize': args.coverage_binsize,
					'task_id': task_id_counter
				})
				task_id_counter += 1
		else:
			spare_breakpoints.extend(breakpoints)
			if len(spare_breakpoints) > split:
				tasks.append({
					'breakpoints': spare_breakpoints,
					'coverage_binsize': args.coverage_binsize,
					'task_id': task_id_counter
				})
				task_id_counter += 1
				spare_breakpoints = []
	# assign remaining spare breakpoints to last task
	if spare_breakpoints:
		tasks.append({
			'breakpoints': spare_breakpoints,
			'coverage_binsize': args.coverage_binsize,
			'task_id': task_id_counter
		})

	return tasks

def generate_call_breakpoint_tasks(clustered_breakpoints, args):
	""""""
	tasks = []
	task_id_counter = 0
	for contig, clusters in clustered_breakpoints.items():
		tasks.append({
			'contig': contig,
			'clusters': clusters,
			'buffer': args.end_buffer,
			'length': args.length,
			'depth': args.min_support,
			'task_id': task_id_counter
		})
		task_id_counter += 1
	return tasks

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

def generate_get_potential_breakpoint_tasks(aln_files, args):
	""" generate get_potential_breakpoint tasks by chunking the genome """
	tasks = []
	task_id_counter = 0
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	for label, aln_file in aln_files.items():
		contigs = helper.get_contigs(args.contigs, args.ref_index)
		if not args.is_cram:
			# if nothing mapped, don't consider contig
			if not args.contigs:
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
						'coverage_binsize': args.coverage_binsize,
						'single_bnd': args.single_bnd,
						'single_bnd_min_length': args.single_bnd_min_length,
						'single_bnd_max_mapq': args.single_bnd_max_mapq,
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
					'coverage_binsize': args.coverage_binsize,
					'single_bnd': args.single_bnd,
					'single_bnd_min_length': args.single_bnd_min_length,
					'single_bnd_max_mapq': args.single_bnd_max_mapq,
					'task_id': task_id_counter
				})
				task_id_counter += 1
	return tasks

def generate_coverage_arrays(aln_files, args):
	""" generate empty numpy arrays """
	coverage_arrays = {}
	contig_lengths = helper.get_contig_lengths(args.ref_index)
	contigs_to_consider = helper.get_contigs(args.contigs, args.ref_index)
	for label, _ in aln_files.items():
		for contig, contig_length in contig_lengths.items():
			if contig not in contigs_to_consider:
				continue
			for haplotype in [1,2,None]:
				coverage_arrays.setdefault(label, {}).setdefault(contig, {})[haplotype] = sharedctypes.RawArray('i', ceil(contig_length/args.coverage_binsize))
	return coverage_arrays

def run_annotate(called_breakpoints, shared_cov_arrays, args):
	""" generate task list for annotating breakpoints and coordinate the execution """
	tasks = generate_annotate_depth_tasks(called_breakpoints, args)
	task_tracker = Array('h', [0]*len(tasks))
	print(f' > {len(tasks)} depth annotation tasks')

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
				proc = Process(target=execute_annotate, args=(task_args, shared_cov_arrays, task_tracker, pipes[i][1]))
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
					del result_conn
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see task is done
					if result_conn.poll(90):
						result = result_conn.recv()
						proc.join(timeout=1) # these hang
						result_conn.close()
						results.append(result)
					else:
						print(f'ERROR: Timeout polling process {i} (task {assigned_task_id})', file=sys.stderr)
						print(f'Unable to annotate depth in given region', file=sys.stderr)
					del process
					del result_conn
					processes[i] = None
					pipes[i] = None
	del task_tracker
	annotated_breakpoints = []
	for result in results:
		annotated_breakpoints.extend(result)

	return annotated_breakpoints

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
					del result_conn
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see task is done
					if result_conn.poll(90):
						result = result_conn.recv()
						proc.join(timeout=1) # these hang
						result_conn.close()
						results.append(result)
					else:
						print(f'ERROR: Timeout polling process {i} (task {assigned_task_id})', file=sys.stderr)
						print(f'Unable to call final breakpoints in given region', file=sys.stderr)
					del process
					del result_conn
					processes[i] = None
					pipes[i] = None
	del task_tracker
	called_breakpoints = {}
	for result in results:
		breakpoints, contig = result
		called_breakpoints.setdefault(contig, []).extend(breakpoints)

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
					del result_conn
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see task is done
					if result_conn.poll(90):
						result = result_conn.recv()
						proc.join(timeout=1) # these processes hang
						result_conn.close()
						results.append(result)
					else:
						print(f'ERROR: Timeout polling process {i} (task {assigned_task_id})', file=sys.stderr)
						print(f'Unable to cluster breakpoints in given region', file=sys.stderr)
					del process
					del result_conn
					processes[i] = None
					pipes[i] = None
	del task_tracker
	clustered_breakpoints = {}
	for result in results:
		contig, clusters = result
		clustered_breakpoints[contig] = clusters

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
				task_args['contig_coverage_array'] = shared_cov_arrays[task_args['label']][task_args['contig']]
				# create a new pipe
				pipes[i] = Pipe(duplex=False)
				# create a new process and replace it (matching the pipe with its process)
				proc = Process(target=execute_get_potential_breakpoint_task, args=(task_args, task_tracker, pipes[i][1]))
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
					del result_conn
					processes[i] = None
					pipes[i] = None
				elif proc.is_alive() and task_tracker[assigned_task_id] != 0:
					# even though process alive, can see that the task is done
					if result_conn.poll(90):
						result = result_conn.recv()
						proc.join(timeout=1) # timeout because these processes hang
						result_conn.close()
						results.append(result)
					else:
						region = f'{regions_running[i][1]}:{regions_running[i][2]}-{regions_running[i][3]}'
						print(f'ERROR: Timeout polling process {i} (task {assigned_task_id})', file=sys.stderr)
						print(f'Unable to fetch results from {regions_running[i][0]} file for region {region}', file=sys.stderr)
						print(f'Please try re-running with a smaller --chunksize argument than {args.chunksize}', file=sys.stderr)
						proc.terminate()
						proc.join()
					del process
					del result_conn
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

	"""
	# TESTING SIZE OF PB OBJECTS
	from pympler import asizeof
	max_size = -1
	max_size_non_insertion = -1
	mean_size = []
	for _, chrom_bps in potential_breakpoints.items():
		for bp in chrom_bps:
			size = asizeof.asizeof(bp)
			mean_size.append(size)
			if size > max_size:
				max_size = size
				print(f'New max size of {round((max_size/1000), 2)}KB {bp.breakpoint_notation}')
			if bp.breakpoint_notation not in ['+', '-', '<INS>'] and size > max_size_non_insertion:
				max_size_non_insertion = size
				print(f'New (non ins) max size of {round((max_size_non_insertion/1000), 2)}KB {bp.breakpoint_notation}')
	print(f'The maximum size of PB objects is {round((max_size/1000), 2)}KB')
	print(f'The maximum size of (non ins) PB objects is {round((max_size_non_insertion/1000), 2)}KB')
	print(f'The mean size is {round((sum(mean_size)/len(mean_size))/1000, 2)}KB for {len(mean_size)} breakpoints')
	"""

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

	# 4) COMPUTE AND ANNOTATE LOCAL DEPTH
	# use converted arrays to annotate depth
	annotated_breakpoints = run_annotate(called_breakpoints, shared_cov_arrays, args)
	helper.time_function("Annotated breakpoints with depth and phasing", checkpoints, time_str)

	print(f'{len(annotated_breakpoints)} annotated breakpoints to output')

	# more cleanup
	del called_breakpoints
	gc.collect()

	# 5) OUTPUT BREAKPOINTS
	# TODO: parallelize this somehow??
	# by chromosome and then concat?
	# define filenames
	vcf_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.vcf')
	bedpe_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints.bedpe')
	tsv_file = os.path.join(outdir, f'{args.sample}.sv_breakpoints_read_support.tsv')
	inserted_seq_fasta_file = os.path.join(outdir, f'{args.sample}.inserted_sequences.fa')
	# build strings
	ref_fasta = pysam.FastaFile(args.ref)
	bedpe_string = ''
	vcf_string = helper.generate_vcf_header(args, annotated_breakpoints[0])
	read_support_string = 'VARIANT_ID\tTUMOUR_SUPPORTING_READS\tNORMAL_SUPPORTING_READS\n'
	insertion_fasta_string = ''
	count = 0
	for bp in annotated_breakpoints:
		bedpe_string += bp.as_bedpe(count)
		vcf_string += bp.as_vcf(ref_fasta)
		read_support_string += bp.as_read_support(count)
		if bp.breakpoint_notation in ["<INS>", "+", "-"]:
			insertion_fasta_string += bp.as_insertion_fasta(count)
		count+=1
	# write output files
	with open(vcf_file, 'w') as output:
		output.write(vcf_string)
	with open(bedpe_file, 'w') as output:
		output.write(bedpe_string)
	with open(tsv_file, 'w') as output:
		output.write(read_support_string)
	with open(inserted_seq_fasta_file, 'w') as output:
		output.write(insertion_fasta_string)
	# sort vcf
	bcftools.sort('-o', vcf_file, vcf_file, catch_stdout=False)
	helper.time_function("Output consensus breakpoints", checkpoints, time_str)

	return checkpoints, time_str

if __name__ == "__main__":
	print("Functions to run SAVANA")
