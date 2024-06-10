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

def create_variant_dicts(vcf_file, label, qual_filter):
	""" given a vcf file, create a dict representation of relevant attributes for each variant """
	variant_dicts = []
	id_count = 0
	for variant in cyvcf2.VCF(vcf_file):
		if (qual_filter and helper.is_int(variant.QUAL) and variant.QUAL >= qual_filter) or not qual_filter or not helper.is_int(variant.QUAL):
			variant_dict = {
				'label': label,
				'id': label+"_"+str(id_count),
				'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
				'start_loc': variant.start,
				'length': variant.INFO.get('SVLEN'),
				'type': variant.INFO.get('SVTYPE'),
				'within_buffer': [],
				'external_id': variant.ID,
				'qual': round(variant.QUAL, 2) if variant.QUAL else variant.QUAL,
				'validated': None
			}
			variant_dicts.append(variant_dict)
			id_count += 1
			if variant.INFO.get("END"):
				# create another breakpoint object for alternate edge
				variant_dict = {
					'label': label,
					'id': label+"_"+str(id_count),
					'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
					'start_loc': variant.END,
					'length': variant.INFO.get('SVLEN'),
					'type': variant.INFO.get('SVTYPE'),
					'within_buffer': [],
					'external_id': variant.ID,
					'qual': round(variant.QUAL, 2) if variant.QUAL else variant.QUAL,
					'validated': None
				}
				variant_dicts.append(variant_dict)
				id_count += 1
			elif variant.INFO.get("END2"):
				# create another breakpoint object for alternate edge
				variant_dict = {
					'label': label,
					'id': label+"_"+str(id_count),
					'start_chr': variant.CHROM2[3:] if variant.CHROM2.startswith('chr') else variant.CHROM2,
					'start_loc': variant.END2,
					'length': variant.INFO.get('SVLEN'),
					'type': variant.INFO.get('SVTYPE'),
					'within_buffer': [],
					'external_id': variant.ID,
					'qual': round(variant.QUAL, 2) if variant.QUAL else variant.QUAL,
					'validated': None
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

	# match variants with each other
	used_compare_set_ids = []
	for in_variant in input_set:
		# sort compare set variants within the buffer
		within_buffer_sorted = sorted(in_variant['within_buffer'], key=lambda x: x[1])
		for within_buffer_compare_variant, distance in within_buffer_sorted:
			if within_buffer_compare_variant['id'] not in used_compare_set_ids:
				in_variant['validated'] = {
					'matched_id': within_buffer_compare_variant['id'],
					'distance': distance,
					'label': within_buffer_compare_variant['label']
				}
				# reciprocally label the compare variant
				for compare_variant in compare_set:
					if compare_variant['id'] == within_buffer_compare_variant['id']:
						compare_variant['validated'] = {
							'matched_id': in_variant['id'],
							'distance': distance,
							'label': in_variant['label']
						}
						break
				used_compare_set_ids.append(within_buffer_compare_variant['id'])
				break

	"""
	# mark input set
	for variant in input_set:
		variant['validated'] = None
		within_buffer_sorted = sorted(variant['within_buffer'], key=lambda x: x[1])
		i = 0
		while not variant['validated'] and i < len(within_buffer_sorted):
			compare_variant, distance = within_buffer_sorted[i]
			if compare_variant['id'] not in used_compare_set_ids:
				variant['validated'] = (compare_variant['id'], distance)
				# now go validate this variant

				used_compare_set_ids.append(compare_variant['id'])
			i+=1

	# mark compare sets
	for variant in compare_set:
		variant['validated'] = None
		within_buffer_sorted = sorted(variant['within_buffer'], key=lambda x: x[1])
		if variant['external_id'] == 'gridss45fb_2672h':
			print(within_buffer_sorted)
		i = 0
		while not variant['validated'] and i < len(within_buffer_sorted):
			input_variant, distance = within_buffer_sorted[i]
			if variant['external_id'] == 'gridss45fb_2672h':
				print(input_variant)
				print(distance)
				print(input_variant['id'] not in used_input_set_ids)
			if input_variant['id'] not in used_input_set_ids:
				variant['validated'] = (input_variant, distance)
				used_input_set_ids.append(input_variant['id'])
				if input_variant['id'] == "ID_118102_2":
					print("validated:")
					print(variant)
			i+=1
	"""

	# SOMATIC CALCULATIONS
	# validated somatic variants (true positives)
	tp = [v for v in input_set if v['validated'] and v['validated']['label'] == "SOMATIC"]
	# unvalidated variants (false positives)
	fp = [v for v in input_set if not v['validated']]
	# report missed somatic variants (false negative)
	fn = [v for v in compare_set if not v['validated'] and v['label'] == 'SOMATIC']
	validation_str.append('\nSTATISTICS FOR SOMATIC VARIANTS')
	if args.curate:
		validation_str.append('(Allowing for at most TWO variants to be validated by the same event)')
	else:
		validation_str.append('(Not allowing for more than ONE variant to be validated by the same event)')

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
	except ZeroDivisionError as _:
		print("WARNING: Unable to calculate validation statistics due to divide by zero exception")

	validation_str.append('\nMISSED VARIANTS')
	for v in fn:
		missed_variant_row = ''
		for key in ['external_id', 'type', 'length', 'qual']:
			missed_variant_row+= f'{v[key]}\t' if key in v and v[key] else ''
		validation_str.append(missed_variant_row)

	with open(args.stats, "w", encoding="utf-8") as stats_file:
		for line in validation_str:
			stats_file.write(line+"\n")

	return

def evaluate_vcf(args, checkpoints, time_str):
	""" given the input, somatic, and germline VCFs, label the input VCF"""
	# create the comparison set from somatic & germline VCFs
	compare_set = create_variant_dicts(args.somatic, 'SOMATIC', args.qual_filter)
	vcfs_string=f'{args.somatic}'
	if args.germline:
		vcfs_string+=f' and {args.germline}'
		compare_set.extend(create_variant_dicts(args.germline, 'GERMLINE', args.qual_filter))

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
		if variant.INFO.get("END"):
			# create another breakpoint object for alternate edge
			input_variants.append({
				'label': 'INPUT',
				'id': str(variant.ID)+"_2",
				'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
				'start_loc': variant.INFO.get('END'),
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
		elif variant.INFO.get("END2"):
			# create another breakpoint object for alternate edge
			input_variants.append({
				'label': 'INPUT',
				'id': str(variant.ID)+"_2",
				'start_chr': variant.CHROM2[3:] if variant.CHROM2.startswith('chr') else variant.CHROM2,
				'start_loc': variant.INFO.get('END2'),
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

	# label input variants with matched somatic/germline
	compare_variants_used = {}
	input_variant_labels = {}
	input_variants = sorted(input_variants, key=lambda x: min(d[1] for d in x['within_buffer']) if x['within_buffer'] else (args.overlap_buffer + 1))
	for variant in input_variants:
		variant['validated'] = False
		closest_variant = None
		closest_value = None
		for compare_variant, distance in variant['within_buffer']:
			# check how many times this variant has been used as a label
			used_count = len(compare_variants_used.get(compare_variant['id'],[]))
			if args.curate and used_count >= 2:
				continue
			elif not args.curate and used_count >= 1:
				# stricter when not curating a training set
				continue
			if args.by_support:
				# tie-break using support - default to zero to prevent no-support match
				closest_value = 0 if not closest_value else closest_value
				if compare_variant['label'] == 'SOMATIC' and variant['tumour_support'] > closest_value and variant['normal_support'] == 0:
					closest_variant = compare_variant
					closest_value = variant['tumour_support']
				elif compare_variant['label'] == 'GERMLINE' and variant['normal_support'] > closest_value:
					closest_variant = compare_variant
					closest_value = variant['normal_support']
			elif args.by_distance:
				# tie-break by closest variant
				if closest_value is None or distance < closest_value:
					closest_variant = compare_variant
					closest_value = distance
		if closest_variant:
			compare_variants_used.setdefault(closest_variant['id'], []).append(variant['id'])
			input_variant_labels[variant['id']] = (closest_variant['label'], closest_variant['external_id'], closest_value)

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
	input_vcf.add_info_to_header({
		'ID': 'LABEL_VARIANT_ID',
		'Type': 'String',
		'Number': '1',
		'Description': "ID of variant used to provide LABEL"
	})
	input_vcf.add_info_to_header({
		'ID': 'DISTANCE_TO_MATCH',
		'Type': 'Integer',
		'Number': '1',
		'Description': "Distance to matched variant in LABEL_VARIANT_ID"
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	for variant in input_vcf:
		label, match_id, distance = input_variant_labels.get(variant.ID, ('NOT_IN_COMPARISON', None, None))
		# if curating, don't allow mates to have incongruous labels
		if args.curate and not match_id:
			mate_id = variant.INFO.get('MATEID', None)
			label, match_id, distance = input_variant_labels.get(mate_id, ('NOT_IN_COMPARISON', None, None)) if mate_id else (label, match_id, distance)
			if match_id:
				match_id = "MATE_"+match_id
		variant.INFO['LABEL'] = label
		if match_id:
			variant.INFO['LABEL_VARIANT_ID'] = match_id
			variant.INFO['DISTANCE_TO_MATCH'] = distance
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
