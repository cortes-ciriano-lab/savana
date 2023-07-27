"""
SAVANA strucural variant caller for long-read data - evaluate sub-command
Created: 21/09/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import re
import os
import csv

import cyvcf2

import savana.helper as helper
from savana.breakpoints import *
from savana.clusters import *

def create_variant_dicts(vcf_file, label):
	""" given a vcf file, create a dict representation of relevant attributes for each variant """
	variant_dicts = []
	id_count = 0
	for variant in cyvcf2.VCF(vcf_file):
		variant_dict = {
			'label': label,
			'id': label+"_"+str(id_count),
			'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
			'start_loc': variant.start,
			'length': variant.INFO.get('SVLEN'),
			'type': variant.INFO.get('SVTYPE'),
			'within_buffer': [],
			'external_id': variant.ID
		}
		variant_dicts.append(variant_dict)
		id_count += 1

	return variant_dicts

def compute_statistics(args, compare_set, input_set, vcfs_string):
	""" given a compare set of variants, compute statistics on number of variants identifed """
	validation_str = []
	break_str = '----------'
	validation_str.append(f'Number of INPUT variants found in "{vcfs_string}": {len([x for x in input_set if x["within_buffer"]])}/{len(input_set)}')
	validation_str.append(f'Number of "SOMATIC" variants found in INPUT {len([x for x in compare_set if (x["label"] == "SOMATIC" and x["within_buffer"])])}/{len([x for x in compare_set if x["label"] == "SOMATIC"])}')
	if args.germline:
		validation_str.append(f'Number of "GERMLINE" variants found in INPUT {len([x for x in compare_set if (x["label"] == "GERMLINE" and x["within_buffer"])])}/{len([x for x in compare_set if x["label"] == "GERMLINE"])}')

	used_compare_set_ids = []
	used_input_set_ids = []

	# mark input set
	for variant in input_set:
		variant['validated'] = None
		within_buffer_sorted = sorted(variant['within_buffer'], key=lambda x: x[1])
		i = 0
		while not variant['validated'] and i < len(within_buffer_sorted):
			compare_variant, distance = within_buffer_sorted[i]
			if compare_variant['id'] not in used_compare_set_ids:
				variant['validated'] = (compare_variant, distance)
				used_compare_set_ids.append(compare_variant['id'])
			i+=1

	# mark compare sets
	for variant in compare_set:
		variant['validated'] = None
		within_buffer_sorted = sorted(variant['within_buffer'], key=lambda x: x[1])
		i = 0
		while not variant['validated'] and i < len(within_buffer_sorted):
			input_variant, distance = within_buffer_sorted[i]
			if input_variant['id'] not in used_input_set_ids:
				variant['validated'] = (input_variant, distance)
				used_input_set_ids.append(input_variant['id'])
			i+=1

	# SOMATIC CALCULATIONS
	# validated somatic variants (true positives)
	tp = [v for v in input_set if v['validated'] and v['validated'][0]['label'] == "SOMATIC"]
	# unvalidated variants (false positives)
	fp = [v for v in input_set if not v['validated']]
	# report missed somatic variants (false negative)
	fn = [v for v in compare_set if not v['validated'] and v['label'] == 'SOMATIC']
	validation_str.append(f'\nSTATISTICS FOR SOMATIC VARIANTS')
	validation_str.append('(Not allowing for two variants to be validated by the same event)')
	validation_str.append(break_str)
	validation_str.append(f'True Positives: {len(tp)}')
	validation_str.append(f'False Positives: {len(fp)}')
	validation_str.append(f'False Negatives: {len(fn)}')
	validation_str.append(break_str)

	try:
		precision = len(tp)/(len(tp)+len(fp))
		recall = len(tp)/(len(tp)+len(fn))
		f_measure = (2*precision*recall)/(precision+recall)
		validation_str.append(f'Precision: {round(precision, 3)}')
		validation_str.append(f'Recall: {round(recall, 3)}')
		validation_str.append(f'F-measure: {round(f_measure, 3)}')
		validation_str.append(break_str)
	except ZeroDivisionError as e:
		print("WARNING: Unable to calculate validation statistics due to divide by zero exception")
	except Exception as e:
		print(f'WARNING: Unable to calculate validation statistics due to "{str(e)}"')

	validation_str.append(f'\nMISSED VARIANTS')
	for v in fn:
		if v["length"]:
			validation_str.append(f'{v["external_id"]}\t{v["type"]}\t{int(v["length"])}')
		else:
			validation_str.append(f'{v["external_id"]}\t{v["type"]}')

	with open(args.stats, "w") as stats_file:
		for line in validation_str:
			stats_file.write(line+"\n")

	return

def evaluate_vcf(args, checkpoints, time_str):
	""" given the input, somatic, and germline VCFs, label the input VCF"""
	# create the comparison set from somatic & germline VCFs
	compare_set = create_variant_dicts(args.somatic, 'SOMATIC')
	vcfs_string=f'{args.somatic}'
	if args.germline:
		vcfs_string+=f' and {args.germline}'
		compare_set.extend(create_variant_dicts(args.germline, 'GERMLINE'))

	# read in input vcf, iterate through compare variants, store those that are within buffer
	input_vcf = cyvcf2.VCF(args.input)
	input_variants = []
	for variant in input_vcf:
		variant_chrom = variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM
		input_variants.append({
			'label': 'INPUT',
			'id': variant.ID,
			'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
			'start_loc': variant.start,
			'end_loc': variant.INFO.get('END'), #only for cuteSV
			'length': variant.INFO.get('SVLEN'),
			'type': variant.INFO.get('SVTYPE'),
			'within_buffer': []})
		if args.by_support:
			input_variants[-1]['tumour_support'] = int(variant.INFO.get('TUMOUR_SUPPORT'))
			input_variants[-1]['normal_support'] = int(variant.INFO.get('NORMAL_SUPPORT'))
		for compare_variant in compare_set:
			if compare_variant['start_chr'] == variant_chrom:
				distance = abs(compare_variant['start_loc'] - input_variants[-1]['start_loc'])
				if distance <= args.overlap_buffer:
					compare_variant['within_buffer'].append((input_variants[-1], distance))
					input_variants[-1]['within_buffer'].append((compare_variant, distance))
				if input_variants[-1]['end_loc']:
					distance = abs(compare_variant['start_loc'] - input_variants[-1]['end_loc'])
					if distance <= args.overlap_buffer:
						compare_variant['within_buffer'].append((input_variants[-1], distance))
						input_variants[-1]['within_buffer'].append((compare_variant, distance))

	# assign matched within buffer variants a label
	input_variant_labels = {}
	for variant in compare_set:
		variant['validated'] = False
		closest_value = None
		closest_variant = None
		for input_variant, distance in variant['within_buffer']:
			if input_variant['id'] in input_variant_labels:
				continue # only allow each variant to be labelled once
			if args.by_support:
				# (tie-break using support, don't allow label without support)
				closest_value = 0 if not closest_value else closest_value
				if variant['label'] == 'SOMATIC' and input_variant['tumour_support'] > closest_value:
					closest_variant = input_variant['id']
					closest_value = input_variant['tumour_support']
				elif variant['label'] == 'GERMLINE' and input_variant['normal_support'] > closest_value:
					closest_variant = input_variant['id']
					closest_value = input_variant['normal_support']
			elif args.by_distance:
				if not closest_value or distance < closest_value:
					closest_value = distance
					closest_variant = input_variant['id']
		if closest_variant:
			# match with the input variant with highest 'label' support
			input_variant_labels[closest_variant] = variant['label']

	# edit input_vcf to include LABEL in header
	input_vcf = cyvcf2.VCF(args.input)
	labels_list = ['SOMATIC', 'GERMLINE', 'NOT_IN_COMPARISON']
	desc_string = str(f'One of \'{"|".join(labels_list)}\' evaluating against {vcfs_string}')
	input_vcf.add_info_to_header({
		'ID': 'LABEL',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	for variant in input_vcf:
		variant.INFO['LABEL'] = input_variant_labels.get(variant.ID, 'NOT_IN_COMPARISON')
		out_vcf.write_record(variant)
	out_vcf.close()
	input_vcf.close()
	helper.time_function("Output labelled VCF", checkpoints, time_str)

	if args.stats:
		compute_statistics(args, compare_set, input_variants, vcfs_string)
		helper.time_function("Wrote statistics file", checkpoints, time_str)

	return

if __name__ == "__main__":
	print("Evaluate functions for SAVANA")
