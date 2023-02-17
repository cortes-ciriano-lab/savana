"""
SAVANA strucural variant caller for long-read data - main program
Created: 21/09/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import sys
import os
import argparse

from time import time
from math import ceil
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pysam

import savana.run as run
import savana.helper as helper
import savana.validation as validation
from savana.breakpoints import *
from savana.clusters import *

logo = """
███████  █████  ██    ██  █████  ███    ██  █████
██      ██   ██ ██    ██ ██   ██ ████   ██ ██   ██
███████ ███████ ██    ██ ███████ ██ ██  ██ ███████
     ██ ██   ██  ██  ██  ██   ██ ██  ██ ██ ██   ██
███████ ██   ██   ████   ██   ██ ██   ████ ██   ██
"""

def savana_run(args):
	""" main function for SAVANA """
	if not args.sample:
		# set sample name to default if req.
		args.sample = os.path.splitext(os.path.basename(args.tumour))[0]
	print(f'Running as sample {args.sample}')
	# create output dir if it doesn't exist
	outdir = os.path.join(os.getcwd(), args.outdir)
	if not os.path.exists(outdir):
		print(f'Creating directory {outdir} to store results')
		os.mkdir(outdir)
	elif os.listdir(outdir):
		sys.exit(f'Output directory "{outdir}" already exists and contains files. Please remove the files or supply a different directory name.')
	# set number of threads to cpu count if none set
	if not args.threads:
		args.threads = cpu_count()
	# read bam files (must have bai)
	bam_files = {
		'tumour': pysam.AlignmentFile(args.tumour, "rb"),
		'normal': pysam.AlignmentFile(args.normal, "rb")
	}
	# confirm ref and ref fasta index exist
	if not os.path.exists(args.ref):
		sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
	elif args.ref_index and not os.path.exists(args.ref_index):
		sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
	elif not os.path.exists(f'{args.ref}.fai'):
		sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
	else:
		args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
		print(f'Using {args.ref_index} as reference fasta index')
	# initialize timing
	checkpoints = [time()]
	time_str = []
	# run SAVANA processes
	consensus_clusters, breakpoints, checkpoints, time_str = run.spawn_processes(args, bam_files, checkpoints, time_str, outdir)
	# write debugging files
	if args.debug:
		run.write_cluster_bed(consensus_clusters, outdir)
		run.calculate_cluster_stats(consensus_clusters, outdir)
	# validate vcf
	output_vcf = os.path.join(outdir, f'{args.sample}.sv_breakpoints.vcf')
	if args.validation:
		try:
			validation.validate_vcf(outdir, output_vcf, args.validation)
		except Exception as e:
			print(f'\nWARNING: Validation of breakpoints against {args.validation} failed due to "{str(e)}"')
			print(f'You can retry by running "python savana/validation.py --outdir testing --input {output_vcf} --validation {args.validation}"')
	# finish timing
	helper.time_function("Total time", checkpoints, time_str, final=True)
	f = open(os.path.join(outdir, 'time.log'), "w+")
	f.write("\n".join(time_str))
	f.write("\n")
	f.close()

def savana_evaluate(args):
	""" main function for savana evaluate """
	# create output dir if it doesn't exist
	outdir = os.path.join(os.getcwd(), args.outdir)
	if not os.path.exists(outdir):
		print(f'Creating directory {outdir} to store results')
		os.mkdir(outdir)
	elif os.listdir(outdir):
		sys.exit(f'Output directory "{outdir}" already exists and contains files. Please remove the files or supply a different directory name.')
	else:
		print(f'Using {outdir} to store output')
	# check input VCFs
	if not os.path.exists(args.input):
		sys.exist(f'Provided input vcf: "{args.input}" does not exist. Please provide full path')
	if not os.path.exists(args.validation):
		sys.exist(f'Provided validation vcf: "{args.validation}" does not exist. Please provide full path')
	# perform validation
	validation.validate_vcf(outdir, args.input, args.validation)
	print("Done.")

def main():
	""" main function for SAVANA - collects command line arguments and executes algorithm """

	# parse arguments - separated into subcommands
	global_parser = argparse.ArgumentParser(prog="savana", description="SAVANA - somatic SV caller")
	global_parser.add_argument('--version', action='version', version=f'SAVANA {helper.__version__}')
	subparsers = global_parser.add_subparsers(title="subcommands", help='SAVANA sub-commands')

	# savana run
	run_parser = subparsers.add_parser("run", help="run SAVANA on tumour and normal long-read BAMs to detect SVs")
	run_parser.add_argument('--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
	run_parser.add_argument('--normal', nargs='?', type=str, required=True, help='Normal BAM file (must have index)')
	run_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
	run_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
	run_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")
	run_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
	run_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='MAPQ filter on reads which are considered (default=5)')
	run_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer to add when clustering adjacent potential breakpoints (default=10)')
	run_parser.add_argument('--depth', nargs='?', type=int, default=3, help='Threshold number of supporting reads required to keep a cluster (default=3)')
	run_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
	run_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
	run_parser.add_argument('--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
	run_parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
	run_parser.add_argument('--validation', nargs='?', type=str, required=False, help='VCF file to use as validation (optional)')
	run_parser.set_defaults(func=savana_run)

	# savana evaluate
	evaluate_parser = subparsers.add_parser("evaluate", help="evaluate and label SAVANA VCF given a VCF to evaluate against")
	evaluate_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
	evaluate_parser.add_argument('--input', nargs='?', type=str, required=False, help='VCF file to evaluate')
	evaluate_parser.add_argument('--validation', nargs='?', type=str, required=False, help='VCF file to evaluate against')
	evaluate_parser.set_defaults(func=savana_evaluate)

	args = global_parser.parse_args()
	print(logo)
	print(f'Version {helper.__version__} - beta')
	src_location = __file__
	print(f'Source: {src_location}\n')
	args.func(args)

if __name__ == "__main__":
	main()
