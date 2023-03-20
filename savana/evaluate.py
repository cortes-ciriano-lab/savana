"""
SAVANA strucural variant caller for long-read data - validation sub-command
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

def create_variant_dicts(vcf_file, label, stats):
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
		}
		if stats:
			variant_dict['within_buffer'] = []
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

	f = open(args.stats, 'w')
	for line in validation_str:
		f.write(line)
	f.close()
	print(f'Evaluation statistics written to "{args.stats}"')

	return

def evaluate_vcf(args):
	""" given the input, somatic, and germline VCFs, label the input VCF"""
	# create the comparison set from somatic & germline VCFs
	compare_set = create_variant_dicts(args.somatic, 'SOMATIC', args.stats)
	vcfs_string=f'{args.somatic}'
	if args.germline:
		vcfs_string+=f' and {args.germline}'
		compare_set.extend(create_variant_dicts(args.germline, 'GERMLINE', args.stats))

	# read in input vcf, add to header, iterate through variants, add label field
	input_vcf = cyvcf2.VCF(args.input)	
	# add LABEL field to header (reference somatic/germline VCF files)
	labels_list = ['SOMATIC', 'GERMLINE', 'FOUND_IN_BOTH', 'NOT_IN_COMPARISON']
	desc_string = str(f'One of \'{"|".join(labels_list)}\' evaluating against {vcfs_string}')
	input_vcf.add_info_to_header({
		'ID': 'LABEL',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	# output vcf using modified input_vcf as template
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.stats:
		input_set = []
	for variant in input_vcf:
		labels = []
		variant_chrom = variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM
		if args.stats:
			input_set.append({
				'label': 'INPUT',
				'id': variant.ID,
				'start_chr': variant.CHROM[3:] if variant.CHROM.startswith('chr') else variant.CHROM,
				'start_loc': variant.start,
				'length': variant.INFO.get('SVLEN'),
				'type': variant.INFO.get('SVTYPE'),
				'within_buffer': []
			})
		for compare_variant in compare_set:
			if compare_variant['start_chr'] == variant_chrom:
				distance = abs(compare_variant['start_loc'] - variant.start)
				if distance <= args.buffer:
					labels.append(compare_variant['label'])
					if args.stats:
						compare_variant['within_buffer'].append((input_set[-1], distance))
						input_set[-1]['within_buffer'].append((compare_variant, distance))
		if not labels:
			variant.INFO['LABEL'] = labels_list[3]
		elif len(set(labels)) > 1:
			variant.INFO['LABEL'] = labels_list[2]
		else:
			variant.INFO['LABEL'] = labels[0]
		out_vcf.write_record(variant)
	out_vcf.close()
	input_vcf.close()

	if args.stats:
		compute_statistics(args, compare_set, input_set, vcfs_string)

	return

def validate_vcf(outdir, compare_vcf, validation_vcf):
	""" compare output vcf with 'truthset' validation vcf """
	validation_str = []
	validation_str.append(f'\nEvaluating compared to provided validation file: \'{validation_vcf}\'')
	truthset = []
	chrom_starts = None
	with open(validation_vcf) as f:
		reader = csv.reader(f, delimiter='\t')
		for line in list(reader):
			if line[0].startswith("#"):
				continue
			sv_type = re.search(r"SVTYPE=([A-Z]*)", line[7])
			try:
				sv_length_match = re.search(r"SVLEN=([0-9]*)", line[7])
				sv_length = sv_length_match.group(1)
			except Exception as e:
				sv_length = 0 # sometimes not available
			if not chrom_starts:
				if line[0].startswith("chr"):
					chrom_starts = True
				else:
					chrom_starts = False
			truthset.append({
				'start_chr': line[0],
				'start_loc': int(line[1]),
				'label': line[2],
				'alt': line[3],
				'info': line[7],
				'type': sv_type.group(1),
				'length': sv_length,
				'incorrectly_categorized': False,
				'seen': False
			})
	compareset = []
	with open(compare_vcf) as f:
		reader = csv.reader(f, delimiter='\t')
		for line in list(reader):
			if line[0].startswith("#"):
				continue
			#sv_type = re.search(r"SVTYPE=([A-Z]*)", line[7])
			sv_length = re.search(r"SVLEN=([0-9]*)", line[7])
			tumour_support = re.search(r"TUMOUR_SUPPORT=([0-9]*)", line[7])
			normal_support = re.search(r"NORMAL_SUPPORT=([0-9]*)", line[7])
			clusters = re.search(r"CLUSTER=([0-9,a-z]*)",line[7])
			bp_notation = re.search(r"BP_NOTATION=([\+,\-,\<INS\>]*)",line[7])
			# convert to same chrom format as truthset
			if chrom_starts and not line[0].startswith('chr'):
				converted_chrom = 'chr'+line[0]
			elif not chrom_starts and line[0].startswith('chr'):
				converted_chrom = line[0][3:]
			else:
				converted_chrom = line[0] # no need to convert
			# put into a dict
			compareset_entry = {
				'start_chr': converted_chrom,
				'start_loc': int(line[1]),
				'label': line[2],
				'alt': line[3],
				'info': line[7],
				'validated': False
			}
			# savana-specific variables
			if clusters:
				compareset_entry['clusters'] = clusters.group(1)
			if bp_notation:
				compareset_entry['bp_type'] = bp_notation.group(1)
			if sv_length:
				compareset_entry['length'] = sv_length.group(1)
				compareset_entry['tumour_support'] = int(tumour_support.group(1)) if tumour_support else 0
				compareset_entry['normal_support'] = int(normal_support.group(1)) if normal_support else 0
			compareset.append(compareset_entry)
	# sort by support so that those with less normal support and more tumour support are seen first
	somatic_compareset = list(sorted(compareset, key=lambda x: (x['normal_support'], (-1)*x['tumour_support'])))
	# get number of true positives, false positives, and false negatives
	tp, fp, fn = [], [], []
	already_validated = {}
	for bp in somatic_compareset:
		# remove chr prefix to compare
		for ts in truthset:
			# compare every somatic variant to each truthset variant
			if ts['start_chr'] == bp['start_chr']:
				if not ts['seen'] and not bp['validated'] and abs(ts['start_loc']-bp['start_loc']) <= 100:
					ts['seen'] = True
					ts['start_compare'] = abs(ts['start_loc']-bp['start_loc'])
					ts['clusters'] = bp['clusters'] if 'clusters' in bp else None
					bp['validated'] = True
					tp.append(bp)
					already_validated.pop(bp['label'], None)
				elif not bp['validated'] and abs(ts['start_loc']-bp['start_loc']) <= 100:
					already_validated[bp['label']] = [ts['label'], ts['clusters']]

	validation_str.append('\nDUPLICATE BREAKPOINTS')
	for bp_label, ts_info in already_validated.items():
		validation_str.append(f'{ts_info[0]} validated by {ts_info[1]} - {bp_label} marked as FP')

	fp = [bp for bp in somatic_compareset if not bp['validated']]
	fn = [ts for ts in truthset if not ts['seen']]

	validation_str.append('\nMISSED BREAKPOINTS')
	for bp in fn:
		validation_str.append(f'{bp["label"]} not validated')

	# calculate support by sv type
	sv_type_counts = {}
	for ts in truthset:
		if ts['type'] not in sv_type_counts:
			counts = {
				'seen': 0 if not ts['seen'] else 1,
				'missed': 0 if ts['seen'] else 1
			}
			sv_type_counts[ts['type']] = counts
		else:
			counts = {
				'seen': sv_type_counts[ts['type']]['seen'] if not ts['seen'] else sv_type_counts[ts['type']]['seen']+1,
				'missed': sv_type_counts[ts['type']]['missed'] if ts['seen'] else sv_type_counts[ts['type']]['missed']+1,
			}
			sv_type_counts[ts['type']] = counts

	validation_str.append('\nEVALUATION OF BREAKPOINTS')
	validation_str.append(f'True Positives: {len(tp)}')
	validation_str.append(f'False Positives: {len(fp)}')
	validation_str.append(f'False Negatives: {len(fn)}')
	try:
		precision = len(tp)/(len(tp)+len(fp))
		recall = len(tp)/(len(tp)+len(fn))
		f_measure = (2*precision*recall)/(precision+recall)
		validation_str.append(f'\nSTATISTICS')
		validation_str.append(f'Precision: {round(precision, 3)}')
		validation_str.append(f'Recall: {round(recall, 3)}')
		validation_str.append(f'F-measure: {round(f_measure, 3)}')
		#validation_str.append(f'{len(miscategorized)} True Variants Miscategorized as Germline')
	except ZeroDivisionError as e:
		print("WARNING: Unable to calculate validation statistics due to divide by zero exception")
	except Exception as e:
		print(f'WARNING: Unable to calculate validation statistics due to "{str(e)}"')

	validation_str.append('\nVALIDATION BY SV TYPE')
	for sv_type, values in sv_type_counts.items():
		pcnt = round(values['seen']/(values['missed']+values['seen'])*100, 2)
		validation_str.append(f'{sv_type}: identified {values["seen"]} of {values["seen"]+values["missed"]} ({pcnt}%)')

	# output validation statistics
	f = open(os.path.join(outdir, f'validation.stats'), "w+")
	for string in validation_str:
		f.write(string+"\n")
	f.close()

	return

if __name__ == "__main__":
	print("Evaluate functions for SAVANA")
