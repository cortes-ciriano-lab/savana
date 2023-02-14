"""
SAVANA strucural variant caller for long-read data - filtering sub-command
Created: 06/02/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import re
import sys
import os
import argparse
import csv

import savana.helper as helper
from savana.breakpoints import *
from savana.clusters import *

logo = """
███████  █████  ██    ██  █████  ███    ██  █████
██      ██   ██ ██    ██ ██   ██ ████   ██ ██   ██
███████ ███████ ██    ██ ███████ ██ ██  ██ ███████
     ██ ██   ██  ██  ██  ██   ██ ██  ██ ██ ██   ██
███████ ██   ██   ████   ██   ██ ██   ████ ██   ██
    _____ ____           
   / __(_) / /____  _____
  / /_/ / / __/ _ \/ ___/
 / __/ / / /_/  __/ /    
/_/ /_/_/\__/\___/_/     
"""

def apply_somatic_filters(validated_breakpoints):
	""" use heuristics to eliminate noise and germline variants """
	somatic_breakpoints_lenient = []
	for bp in validated_breakpoints:
		if bp.support['tumour'] == 0:
			continue
		elif bp.support['normal']/bp.support['tumour'] < 0.1:
			originating_cluster_stats = bp.originating_cluster.get_stats()
			if any(originating_cluster_stats.values()) == None:
					# unable to evaluate - skip for now
					continue	
			end_cluster_stats = bp.end_cluster.get_stats()
			if originating_cluster_stats['starts_std_dev'] < 150 and originating_cluster_stats['event_heuristic'] < 3:
				if bp.breakpoint_notation == "<INS>" and bp.support['tumour'] > 25:
					somatic_breakpoints_lenient.append(bp)
				elif bp.breakpoint_notation != "<INS>" and bp.support['tumour'] > 5:
					somatic_breakpoints_lenient.append(bp)

	somatic_breakpoints_strict = []
	for bp in validated_breakpoints:
		if bp.support['normal'] > 0:
			continue
		if bp.support['tumour'] > 7:
			originating_cluster_stats = bp.originating_cluster.get_stats()
			if any(originating_cluster_stats.values()) == None:
					# unable to evaluate - skip for now
					continue
			if any(end_cluster_stats.values()) == None:
					# unable to evaluate - skip for now
					continue
			if bp.support['tumour'] > 12 and originating_cluster_stats['uncertainty'] <= 15:
				if originating_cluster_stats['event_heuristic'] <= 0.025:
					somatic_breakpoints_strict.append(bp)
					continue
			elif bp.support['tumour'] > 12 and end_cluster_stats['uncertainty'] <= 30:
				somatic_breakpoints_strict.append(bp)
				continue
			elif end_cluster_stats['uncertainty'] <= 10:
				somatic_breakpoints_strict.append(bp)
				continue

	return somatic_breakpoints_lenient, somatic_breakpoints_strict

def filter_variants():
	somatic_breakpoints_lenient, somatic_breakpoints_strict = apply_somatic_filters(breakpoints)
	# output lenient vcf
	lenient_vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, breakpoints[0])
	for bp in somatic_breakpoints_lenient:
		lenient_vcf_string += bp.as_vcf(args.ref)
	with open(os.path.join(outdir, "somatic.sv_breakpoints.lenient.vcf"), 'w') as output:
		output.write(lenient_vcf_string)
	# output strict vcf
	strict_vcf_string = helper.generate_vcf_header(args.ref, args.ref_index, args.tumour, breakpoints[0])
	for bp in somatic_breakpoints_strict:
		strict_vcf_string += bp.as_vcf(args.ref)
	with open(os.path.join(outdir, "somatic.sv_breakpoints.strict.vcf"), 'w') as output:
		output.write(strict_vcf_string)

def main():
	""" main function for SAVANA validation - collects command line arguments and executes algorithm """
	parser = argparse.ArgumentParser(description="SAVANA - Filter")
	parser.add_argument('--input', nargs='?', type=str, required=False, help='VCF file to filter')
	parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
	parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
	parser.add_argument('--parameters', nargs='?', type=str, required=False, help='Parameters to use for filtering (optional)')
	parser.add_argument('--output', nargs='?', required=True, help='Output filtered VCF file name (default=filtered.vcf)')
	args = parser.parse_args()

	print(logo)
	print(f'Version {helper.__version__} - beta')
	src_location = __file__
	print(f'Source: {src_location}\n')

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
	
	# confirm output filtered VCF file does not exist
	if os.path.if_file(args.output):
		sys.exit(f'Output file "{args.output}" already exists. Please remove it or supply a different output filename.')

	filter_variants()

if __name__ == "__main__":
	main()
