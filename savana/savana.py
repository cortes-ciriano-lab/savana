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
import json

from time import time
from multiprocessing import cpu_count

import pysam

import savana.run as run
import savana.evaluate as evaluate
import savana.train as train
import savana.classify as classify
import savana.helper as helper

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
	outdir = helper.check_outdir(args.outdir)
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
	# finish timing
	helper.time_function("Total time to call raw variants", checkpoints, time_str, final=True)

def savana_classify(args):
	""" main function for savana classify """
	# initialize timing
	checkpoints = [time()]
	time_str = []
	if args.legacy:
		classify.classify_legacy(args, checkpoints, time_str)
	elif args.params:
		classify.classify_by_params(args, checkpoints, time_str)
	elif args.model:
		classify.classify_by_model(args, checkpoints, time_str)
	# finish timing
	helper.time_function("Total time to classify variants", checkpoints, time_str, final=True)

def savana_evaluate(args):
	""" main function for savana evaluate """
	# check input VCFs
	vcf_string = ''
	if not os.path.exists(args.input):
		sys.exit(f'Provided input vcf: "{args.input}" does not exist. Please provide full path')
	if not os.path.exists(args.somatic):
		sys.exit(f'Provided somatic VCF: "{args.somatic}" does not exist. Please provide full path')
	else:
		vcf_string += f'somatic vcf: "{args.somatic}"'
	if args.germline:
		if not os.path.exists(args.germline):
			sys.exit(f'Provided germline VCF: "{args.germline}" does not exist. Please provide full path')
		else:
			vcf_string += f' and germline vcf: "{args.germline}"'

	# initialize timing
	checkpoints = [time()]
	time_str = []
	# perform validation
	evaluate.evaluate_vcf(args, checkpoints, time_str)
	# finish timing
	helper.time_function("Total time to evaluate variants", checkpoints, time_str, final=True)

def savana_train(args):
	""" main function for savana train """
	outdir = helper.check_outdir(args.outdir)
	data_matrix = None
	if args.vcfs:
		# read in and create matrix from VCF files
		data_matrix = train.read_vcfs(args)
	elif args.load_matrix:
		# load data matrix from pickle file
		data_matrix = train.load_matrix(args)
	features, target = train.prepare_data(data_matrix, multiclass=args.multiclass)
	classifier = train.fit_classifier(features, target, outdir, split=0.2, downsample=args.downsample, hyperparameter=args.hyper, multiclass=args.multiclass)
	train.save_model(args, classifier, outdir)

def savana_main(args):
	""" default workflow for savana: savana_run, savana_classify, savana_evaluate """
	# call raw breakpoints
	savana_run(args)
	# get the input VCF for classification
	args.vcf=os.path.join(args.outdir,f'{args.sample}.sv_breakpoints.vcf')
	if not args.model and not args.params:
		# use ONT somatic only model by default
		from pathlib import Path
		args.model = os.path.join(Path(__file__).parent.parent.absolute(),'package_models/ONT-somatic-only.pkl')
		print(f'Using ONT somatic only model "{args.model}" to classify variants')
	# set the output VCF location
	args.output = os.path.join(args.outdir,f'{args.sample}.classified.sv_breakpoints.vcf')
	if not args.somatic_output:
		args.somatic_output = os.path.join(args.outdir,f'{args.sample}.somatic.vcf')
	savana_classify(args)
	if args.somatic:
		# evaluate against somatic/germline VCFs
		# set the input vcf as previous output
		args.input = args.output
		args.output = os.path.join(args.outdir,f'{args.sample}.evaluation.sv_breakpoints.vcf')
		args.stats = os.path.join(args.outdir,f'{args.sample}.evaluation.stats')
		savana_evaluate(args)


def main():
	""" main function for SAVANA - collects command line arguments and executes algorithm """

	# parse arguments - separated into subcommands
	global_parser = argparse.ArgumentParser(prog="savana", description="SAVANA - somatic SV caller")
	global_parser.add_argument('--version', action='version', version=f'SAVANA {helper.__version__}')
	subparsers = global_parser.add_subparsers(title="subcommands", help='SAVANA sub-commands', dest='command')

	# savana run
	run_parser = subparsers.add_parser("run", help="run SAVANA on tumour and normal long-read BAMs to detect SVs")
	run_parser.add_argument('--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
	run_parser.add_argument('--normal', nargs='?', type=str, required=True, help='Normal BAM file (must have index)')
	run_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
	run_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
	run_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")
	run_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
	run_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='MAPQ filter on reads which are considered (default=5)')
	run_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer when clustering adjacent potential breakpoints, excepting insertions (default=10)')
	run_parser.add_argument('--insertion_buffer', nargs='?', type=int, default=100, help='Buffer when clustering adjacent potential insertion breakpoints (default=100)')
	run_parser.add_argument('--depth', nargs='?', type=int, default=3, help='Minumum number of supporting reads from tumour OR normal to consider variant (default=3)')
	run_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
	run_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
	run_parser.add_argument('--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
	run_parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
	run_parser.set_defaults(func=savana_run)

	# savana classify
	classify_parser = subparsers.add_parser("classify", help="classify VCF using model")
	group = classify_parser.add_mutually_exclusive_group()
	group.add_argument('--model', nargs='?', type=str, required=False, help='Pickle file of machine-learning model')
	group.add_argument('--params', nargs='?', type=str, required=False, help='JSON file of custom filtering parameters')
	group.add_argument('--legacy', action='store_true', help='Legacy lenient/strict filtering')
	classify_parser.add_argument('--vcf', nargs='?', type=str, required=True, help='VCF file to classify')
	classify_parser.add_argument('--output', nargs='?', type=str, required=True, help='Output VCF with PASS columns and CLASS added to INFO')
	classify_parser.add_argument('--somatic_output', nargs='?', type=str, required=False, help='VCF with only PASS somatic variants')
	classify_parser.set_defaults(func=savana_classify)

	# savana evaluate
	evaluate_parser = subparsers.add_parser("evaluate", help="label SAVANA VCF with somatic/germline/missing given VCF(s) to compare against")
	evaluate_parser.add_argument('--input', nargs='?', type=str, required=True, help='VCF file to evaluate')
	evaluate_parser.add_argument('--somatic', nargs='?', type=str, required=True, help='Somatic VCF file to evaluate against')
	evaluate_parser.add_argument('--germline', nargs='?', type=str, required=False, help='Germline VCF file to evaluate against (optional)')
	evaluate_parser.add_argument('--overlap_buffer', nargs='?', type=int, default=100, help='Buffer for considering an overlap (default=100)')
	evaluate_parser.add_argument('--output', nargs='?', type=str, required=True, help='Output VCF with LABEL added to INFO')
	evaluate_parser.add_argument('--stats', nargs='?', type=str, required=False, help='Output file for statistics on comparison if desired (stdout otherwise)')
	group = evaluate_parser.add_mutually_exclusive_group()
	group.add_argument('--by_support', action='store_true', help='Comparison method: tie-break by read support')
	group.add_argument('--by_distance', action='store_true', default=True, help='Comparison method: tie-break by min. distance (default)')
	evaluate_parser.set_defaults(func=savana_evaluate)

	# savana train
	train_parser = subparsers.add_parser("train", help="train model on folder of input VCFs")
	group = train_parser.add_mutually_exclusive_group()
	group.add_argument('--vcfs', nargs='?', type=str, required=False, help='Folder of labelled VCF files to read in')
	group.add_argument('--load_matrix', nargs='?', type=str, required=False, help='Pre-loaded pickle file of VCFs')
	#TODO: implement the recurisve functionality lol whoops
	train_parser.add_argument('--downsample', nargs='?', type=float, default=0.5, help='Fraction to downsample majority class by')
	train_parser.add_argument('--recursive', action='store_true', help='Search recursively through input folder for input VCFs')
	train_parser.add_argument('--multiclass', action='store_true', help='Split the GERMLINE and SOMATIC labels so there are three classes')
	train_parser.add_argument('--hyper', action='store_true', help='Perform a randomised search on hyper parameters and use best')
	train_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
	train_parser.add_argument('--save_matrix', nargs='?', type=str, required=False, help='Output pickle file for data matrix of VCFs')
	train_parser.set_defaults(func=savana_train)

	try:
		global_parser.exit_on_error = False
		subparser = global_parser.parse_args().command
	except argparse.ArgumentError as e:
		# unable to parse args, set args to None
		subparser = None
	if not subparser:
		# arguments for default, main savana process: run, classify, evaluate
		global_parser.add_argument('--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
		global_parser.add_argument('--normal', nargs='?', type=str, required=True, help='Normal BAM file (must have index)')
		global_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
		global_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
		global_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")
		global_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
		global_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='MAPQ filter on reads which are considered (default=5)')
		global_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer when clustering adjacent potential breakpoints, excepting insertions (default=10)')
		global_parser.add_argument('--depth', nargs='?', type=int, default=3, help='Minumum number of supporting reads from tumour OR normal to consider variant (default=3)')
		global_parser.add_argument('--insertion_buffer', nargs='?', type=int, default=100, help='Buffer when clustering adjacent potential insertion breakpoints (default=100)')
		global_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
		global_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
		global_parser.add_argument('--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
		global_parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
		# classify args
		classify_group = global_parser.add_mutually_exclusive_group()
		classify_group.add_argument('--model', nargs='?', type=str, required=False, help='Pickle file of machine-learning model to classify with (default=ONT-somatic-only)')
		classify_group.add_argument('--params', nargs='?', type=str, required=False, help='JSON file of custom filtering parameters')
		classify_group.add_argument('--legacy', nargs='?', type=str, required=False, help='Legacy lenient/strict filtering')
		global_parser.add_argument('--somatic_output', nargs='?', type=str, required=False, help='VCF with only PASS somatic variants')
		# evaluate args
		global_parser.add_argument('--somatic', nargs='?', type=str, required=False, help='Somatic VCF file to evaluate against')
		global_parser.add_argument('--germline', nargs='?', type=str, required=False, help='Germline VCF file to evaluate against (optional)')
		global_parser.add_argument('--overlap_buffer', nargs='?', type=int, default=100, required=False, help='Buffer for considering an overlap (default=100)')
		evaluate_group = global_parser.add_mutually_exclusive_group()
		evaluate_group.add_argument('--by_support', action='store_true', help='Comparison method: tie-break by read support')
		evaluate_group.add_argument('--by_distance', action='store_true', default=True, help='Comparison method: tie-break by min. distance (default)')
		global_parser.set_defaults(func=savana_main)
		args = global_parser.parse_args()
	else:
		global_parser.exit_on_error = True
		args = global_parser.parse_args()

	print(logo)
	print(f'Version {helper.__version__} - beta')
	src_location = __file__
	print(f'Source: {src_location}\n')
	args.func(args)

if __name__ == "__main__":
	main()
