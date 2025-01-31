"""
SAVANA strucural variant caller for long-read data - train sub-command
Created: 19/04/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

from multiprocessing import Pool
from functools import partial
from collections import ChainMap
from statistics import mean
from scipy import stats

import pickle
import json
import pandas as pd
import numpy as np
import cyvcf2
import csv
import os
import re

import savana.train as train
from savana.breakpoints import *
from savana.clusters import *

def pool_predict(data_matrix, features_to_drop, model, threads, confidence=None):
	""" split the prediction across multiprocessing Pool """
	if confidence is not None:
		print(f' > Using Mondrian Conformal Prediction with Confidence of {confidence}')
		# untar nonconform scores if not already done
		models_dir = os.path.join(os.path.dirname(__file__),'models')
		for key in ['0','1']:
			txt_path = os.path.join(models_dir, f'{str(key)}_mondrian_dist.txt')
			tar_path = os.path.join(models_dir, f'{str(key)}_mondrian_dist.tar.gz')
			if os.path.isfile(txt_path):
				print(f'Using untarred {txt_path}')
			elif os.path.isfile(tar_path):
				import tarfile
				print(f'First time using MCP - will untar {tar_path}')
				tar = tarfile.open(tar_path, "r:gz")
				tar.extractall(models_dir)
				tar.close()
	data_matrices = np.array_split(data_matrix, threads)
	pool = Pool(threads)
	results = pool.map(partial(
		predict,
		features_to_drop=features_to_drop,
		model=model,
		confidence=confidence
	), data_matrices)
	pool.close()
	pool.join()

	return dict(ChainMap(*results))

def predict(data_matrix, features_to_drop, model, confidence):
	""" given features and a model, return the predictions """
	ids = data_matrix['ID']
	features = data_matrix.drop(features_to_drop, axis=1, errors='ignore')
	# reorder the columns based on model order
	feature_order = model.feature_names_in_
	features = features[feature_order]

	if confidence is not None:
		# mondrian
		predicted_probabilities = model.predict_proba(features)
		# get nonconform scores
		models_dir = os.path.join(os.path.dirname(__file__),'models')
		sig_level = 1.0 - confidence
		prediction_dict_mondrian = {}
		mondrian_dists = {}
		for key in [0, 1]:
			mondrian_dist_path = os.path.join(models_dir, f'{str(key)}_mondrian_dist.txt')
			mondrian_dists[key] = np.loadtxt(mondrian_dist_path, delimiter='\t', skiprows=1)
		for i, prob in enumerate(predicted_probabilities):
			for key in [0, 1]:
				#fraction = len(mondrian_dists[key][mondrian_dists[key] <= prob[key]])/len(mondrian_dists[key])
				p_value = stats.percentileofscore(mondrian_dists[key],  prob[key])
				if p_value > sig_level:
					existing_variant_predictions = prediction_dict_mondrian.get(ids.iloc[i], [])
					prediction_dict_mondrian[ids.iloc[i]] = existing_variant_predictions
					prediction_dict_mondrian[ids.iloc[i]] += [key]

		return prediction_dict_mondrian

	# perform prediction using model with standard prediction
	predicted = model.predict(features)
	prediction_dict = {}
	for i, value in enumerate(predicted):
		prediction_dict[ids.iloc[i]] = value

	return prediction_dict

def parse_cna(cna_file):
	""" using the cna path, store the copy number changepoints """
	cna_dict = {}
	with open(cna_file) as f:
		reader = csv.reader(f, delimiter="\t")
		next(reader, None) # skip header
		curr_chrom = None
		curr_copy_num = None
		curr_seg = None
		curr_minor_copy_num = None
		curr_seg_start = 0
		for row in reader:
			chrom, start, end, seg, _, _, _, copy_num, minor_copy_num, _, _ = row
			copy_num = round(float(copy_num))
			if not curr_copy_num or curr_chrom != chrom:
				# first segment of chromosome
				curr_chrom = chrom
				curr_copy_num = copy_num
				curr_seg = seg
				curr_seg_start = 0
			elif copy_num != curr_copy_num or minor_copy_num != curr_minor_copy_num:
				# change in total or minor copy number
				cna_dict.setdefault(chrom, []).append(
					[curr_seg_start, start, curr_seg]
				)
				# now reset
				curr_copy_num = copy_num
				curr_seg = seg
				curr_minor_copy_num = minor_copy_num
				curr_seg_start = start

	return cna_dict

def rescue_cna(args, checkpoints, time_str):
	""" update the filter column of variants which were rejected but have a supporting change in copy number """
	cna_dict = parse_cna(args.cna_rescue)
	classified_vcf = cyvcf2.VCF(args.output)
	rescue_dict = {} # store the id and distance of the rescued variants for each segment
	for variant in classified_vcf:
		if not args.tumour_only:
			# not tumour only can use NORMAL fields
			if variant.INFO['TUMOUR_READ_SUPPORT'] < args.min_support:
				continue
			elif variant.INFO['NORMAL_READ_SUPPORT'] != 0:
				continue
			elif int(variant.INFO['CLUSTERED_READS_NORMAL']) >= 3:
				continue
			elif variant.CHROM not in cna_dict:
				continue
			else:
				for seg_start, seg_end, seg_name in cna_dict[variant.CHROM]:
					# segment start
					distance_start = abs(variant.start - int(seg_start))
					distance_end = abs(variant.start - int(seg_end))
					if distance_start <= args.cna_rescue_distance:
						rescue_dict.setdefault(seg_name, []).append((variant.ID, distance_start))
					elif distance_end <= args.cna_rescue_distance:
						rescue_dict.setdefault(seg_name, []).append((variant.ID, distance_end))
		else:
			# tumour only - can't use NORMAL fields
			if variant.INFO['TUMOUR_READ_SUPPORT'] < args.min_support:
				continue
			elif variant.CHROM not in cna_dict:
				continue
			else:
				for seg_start, seg_end, seg_name in cna_dict[variant.CHROM]:
					# segment start
					distance_start = abs(variant.start - int(seg_start))
					distance_end = abs(variant.start - int(seg_end))
					if distance_start <= args.cna_rescue_distance:
						rescue_dict.setdefault(seg_name, []).append((variant.ID, distance_start))
					elif distance_end <= args.cna_rescue_distance:
						rescue_dict.setdefault(seg_name, []).append((variant.ID, distance_end))

	# sort by the distance
	rescue_dict = {k: sorted(v, key=lambda x: x[1]) for k, v in rescue_dict.items()}
	# flip so that the key is the variant id
	variant_rescue_dict = {v[0][0]: k for k,v in rescue_dict.items()}

	# open vcf(s) for reading/writing
	in_vcf = cyvcf2.VCF(args.output)
	desc_string = str('Copy number segment used to rescue the SV')
	in_vcf.add_info_to_header({
		'ID': 'CNA_SEG',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	cna_rescue_vcf = os.path.splitext(args.output)[0]+'.cna_rescue.vcf'
	cna_out_vcf = cyvcf2.Writer(cna_rescue_vcf, in_vcf)
	if args.somatic_output:
		cna_rescue_somatic_vcf = os.path.splitext(args.somatic_output)[0]+'.cna_rescue.vcf'
		cna_out_somatic_vcf = cyvcf2.Writer(cna_rescue_somatic_vcf, in_vcf)
	for variant in in_vcf:
		if not variant.FILTER and variant.ID in variant_rescue_dict:
			print(f'Would rescue {variant.ID} but it already passes!')
		elif variant.FILTER and variant.ID in variant_rescue_dict:
			variant.FILTER = 'PASS'
			variant.INFO['CLASS'] = 'RESCUED_CNA'
			seg = variant_rescue_dict[variant.ID]
			variant.INFO['CNA_SEG'] = seg
		elif variant.FILTER and variant.INFO.get('MATEID', None) in variant_rescue_dict:
			variant.FILTER = 'PASS'
			variant.INFO['CLASS'] = 'RESCUED_CNA'
			seg = variant_rescue_dict[variant.INFO.get('MATEID', None)]
			variant.INFO['CNA_SEG'] = seg
		cna_out_vcf.write_record(variant)
		if args.somatic_output and not variant.FILTER or variant.FILTER == 'PASS':
			# add to the somatic only VCF
			cna_out_somatic_vcf.write_record(variant)

	cna_out_vcf.close()
	if args.somatic_output:
		cna_out_somatic_vcf.close()

	helper.time_function("Rescued somatic SVs from CNAs", checkpoints, time_str)

	return

def legacy_pass_strict(variant, event_heuristic):
	""" apply legacy filter for strict thresholds """
	if variant['NORMAL_READ_SUPPORT'] > 0:
		return False
	if variant['TUMOUR_READ_SUPPORT'] > 7:
		origin_uncertainty = (variant['ORIGIN_STARTS_STD_DEV']+1) * (variant['ORIGIN_EVENT_SIZE_STD_DEV']+1)
		end_uncertainty = (variant['END_STARTS_STD_DEV']+1) * (variant['END_EVENT_SIZE_STD_DEV']+1)
		try:
			if variant['TUMOUR_READ_SUPPORT'] > 12 and origin_uncertainty <= 15:
				if event_heuristic <= 0.025:
					return True
			elif variant['TUMOUR_READ_SUPPORT'] > 12 and end_uncertainty <= 30:
				return True
			elif end_uncertainty <= 10:
				return True
		except Exception as _:
			# in case of None
			return False

def legacy_pass_lenient(variant, event_heuristic):
	""" apply legacy filter for lenient thresholds """
	if variant['TUMOUR_READ_SUPPORT'] == 0:
		return False
	elif variant['NORMAL_READ_SUPPORT']/variant['TUMOUR_READ_SUPPORT'] < 0.1:
		try:
			if variant['ORIGIN_STARTS_STD_DEV'] < 150 and event_heuristic < 3:
				if variant['BP_NOTATION'] == "<INS>" and variant['TUMOUR_READ_SUPPORT'] > 25:
					return True
				elif variant['BP_NOTATION'] != "<INS>" and variant['TUMOUR_READ_SUPPORT'] > 5:
					return True
		except Exception as _:
			# in case of None for any stat
			return False

def classify_legacy(args, checkpoints, time_str):
	""" classify using legacy lenient/strict filters """
	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix, args.tumour_only)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	strict_ids = {}
	lenient_ids = {}
	for _, variant in data_matrix.iterrows():
		if variant['ORIGIN_EVENT_SIZE_MEDIAN'] > 0:
			event_heuristic = variant['ORIGIN_EVENT_SIZE_STD_DEV']/variant['ORIGIN_EVENT_SIZE_MEDIAN']
		else:
			event_heuristic = None
		# apply filters
		if legacy_pass_strict(variant, event_heuristic):
			strict_ids[variant['ID']] = True
		if legacy_pass_lenient(variant, event_heuristic):
			lenient_ids[variant['ID']] = True

	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str(f'Variant class prediction from legacy strict/lenient filters {args.model}')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	# create filenames based on output
	if ".vcf" in args.output:
		strict_vcf_filename = args.output.replace(".vcf", ".strict.vcf")
		lenient_vcf_filename = args.output.replace(".vcf", ".lenient.vcf")
	else:
		strict_vcf_filename = args.output+".strict.vcf"
		lenient_vcf_filename = args.output+".lenient.vcf"
	strict_vcf = cyvcf2.Writer(strict_vcf_filename, input_vcf)
	lenient_vcf = cyvcf2.Writer(lenient_vcf_filename, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		passed_strict = strict_ids.get(variant_id, None)
		passed_lenient = lenient_ids.get(variant_id, None)
		if passed_strict:
			# update the INFO field if all sanity checks pass
			variant.INFO['CLASS'] = 'PASSED_SOMATIC_STRICT'
			# add to the somatic only VCF
			strict_vcf.write_record(variant)
		if passed_lenient:
			# update the INFO field if all sanity checks pass
			variant.INFO['CLASS'] = 'PASSED_SOMATIC_LENIENT'
			# add to the somatic only VCF
			lenient_vcf.write_record(variant)
	strict_vcf.close()
	lenient_vcf.close()
	input_vcf.close()

	helper.time_function("Output strict/lenient VCFs", checkpoints, time_str)

	return

def filter_with_comparator(variant_value, filter_value, comparator):
	""" given the variant's value, the filter, and the comparator (min/max), return whether it passes """
	# cast the filter value to the same type as the variant's value
	filter_value = type(variant_value)(filter_value)
	if comparator.upper() == "MAX":
		if variant_value > filter_value:
			return False
	elif comparator.upper() == "MIN":
		if variant_value < filter_value:
			return False
	else:
		raise ValueError(f'Comparator value "{comparator}" unrecognized. Must be one of "MIN" or "MAX"')

	return True

def classify_by_params(args, checkpoints, time_str):
	# read VCF into a dataframe, classify using a parameter JSON
	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix, args.tumour_only)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	# Read in the if/else statements and apply them
	filters = None
	with open(args.custom_params) as json_file:
		filters = json.load(json_file)

	# create dicts to store the VCF output and passing variants
	category_dicts = {}
	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str(f'Variant class as defined in params JSON {args.custom_params}')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	for category in filters.keys():
		filename = args.output.replace(".vcf", f'.{category}.vcf')
		category_dicts[category] = {
			'vcf': cyvcf2.Writer(filename, input_vcf),
			'passed_ids': {}
		}

	# use the matrix-representation of the variants to determine which ones pass the filter(s)
	for i, variant in data_matrix.iterrows():
		for category, filter_dict in filters.items():
			passed_filters = True
			for comp_field, filter_value in filter_dict.items():
				# separate MIN/MAX from field
				comparator, field = comp_field.split("_", 1)
				# if fails one filter, set to false
				passed_filters = False if not filter_with_comparator(variant[field], filter_value, comparator) else passed_filters
			if passed_filters:
				category_dicts[category]['passed_ids'][variant['ID']] = True

	# now, write out the vcf(s)
	for variant in input_vcf:
		variant_id = variant.ID
		for category, category_dict in category_dicts.items():
			if variant_id in category_dict['passed_ids']:
				variant.INFO['CLASS'] = 'PASSED_'+(category).upper()
				category_dict['vcf'].write_record(variant)

	for _, category_dict in category_dicts.items():
		category_dict['vcf'].close()

	return

def pacbio_pass_somatic(variant, min_support, min_af):
	""" custom manual filters for PacBio """
	if variant.INFO['NORMAL_READ_SUPPORT'] != 0:
		return False
	if variant.INFO['TUMOUR_READ_SUPPORT'] < min_support:
		return False
	if variant.INFO['ORIGIN_STARTS_STD_DEV'] > 50.0 or variant.INFO['END_STARTS_STD_DEV'] > 50.0:
		return False
	if variant.INFO['ORIGIN_MAPQ_MEAN'] < 40.0 or variant.INFO['END_MAPQ_MEAN'] < 40.0:
		return False
	if variant.INFO['ORIGIN_EVENT_SIZE_STD_DEV'] > 60.0 or variant.INFO['END_EVENT_SIZE_STD_DEV'] > 60.0:
		return False
	if variant.INFO['TUMOUR_AF'][0] < min_af:
		return False
	if variant.INFO['CLUSTERED_READS_NORMAL'] > 3:
		return False
	return True

def classify_pacbio(args, checkpoints, time_str):
	""" classify using PacBio filters (manually defined) """

	# increase min support and AF for PacBio (not below supplied mins)
	pacbio_min_support = max(7, args.min_support)
	pacbio_min_af = max(0.10, args.min_af)
	print(f'Using PacBio minimum support filters of {pacbio_min_support} reads and {pacbio_min_af} allele-fracion')

	if args.predict_germline:
		print('Germline PacBio filters not yet implemented.')

	input_vcf = cyvcf2.VCF(args.vcf)
	pass_dict = {}
	for variant in input_vcf:
		variant_id = variant.ID
		passed_somatic = pacbio_pass_somatic(variant, pacbio_min_support, pacbio_min_af)
		# update the pass dictionary
		pass_dict[variant_id] = passed_somatic
		if passed_somatic:
			# update mate as well
			pass_dict[variant.INFO.get('MATEID', None)] = passed_somatic

	# now that everything and mate has been passed, output files
	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str('Variant class prediction using PacBio filters')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	desc_string = str('Variant did not pass PacBio filters')
	input_vcf.add_filter_to_header({
		'ID': 'FAIL',
		'Description': desc_string
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.somatic_output:
		somatic_vcf = cyvcf2.Writer(args.somatic_output, input_vcf)
	if args.germline_output:
		print('Germline PacBio filters not yet implemented.')
		#germline_vcf = cyvcf2.Writer(args.germline_output, input_vcf)
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.somatic_output:
		somatic_vcf = cyvcf2.Writer(args.somatic_output, input_vcf)
	if args.germline_output:
		print('Germline PacBio filters not yet implemented.')
		#germline_vcf = cyvcf2.Writer(args.germline_output, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		if pass_dict[variant.ID]:
			# update the INFO field
			variant.INFO['CLASS'] = 'SOMATIC'
			if args.somatic_output:
				# add to somatic-only vcf
				somatic_vcf.write_record(variant)
		else:
			variant.FILTER = 'FAIL'
		out_vcf.write_record(variant)
	out_vcf.close()
	if args.somatic_output:
		somatic_vcf.close()
	input_vcf.close()

	helper.time_function("Output PacBio somatic VCF", checkpoints, time_str)

	return

def pass_rescue(variant):
	""" rescue variants predicted noise by model if they pass set of filters """
	# regardless of phase reject on these criteria
	if variant.INFO['NORMAL_READ_SUPPORT'] != 0:
		return False
	if variant.INFO['ORIGIN_MAPQ_MEAN'] < 15.0:
		return False
	if variant.INFO['TUMOUR_ALT_HP'][0] > 0 and variant.INFO['TUMOUR_ALT_HP'][1] > 0:
		# inconsistent phasing
		return False
	# is there high read support from a single haplotype?
	is_phased = True if (variant.INFO['TUMOUR_ALT_HP'][0] > 5 and variant.INFO['TUMOUR_ALT_HP'][1] == 0) else False
	is_phased = True if (variant.INFO['TUMOUR_ALT_HP'][0] == 0 and variant.INFO['TUMOUR_ALT_HP'][1] > 5) else is_phased
	if is_phased:
		# allow maximum one clustered normal read
		if variant.INFO['CLUSTERED_READS_NORMAL'] > 1:
			return False
		# must have more than 10% AF
		if variant.INFO['TUMOUR_AF'][0] <= 0.10:
			return False
	else:
		# if not high phased-reads support, require 7 reads of support
		if variant.INFO['TUMOUR_READ_SUPPORT'] < 7:
			return False
		# don't allow any clustered normal reads
		if variant.INFO['CLUSTERED_READS_NORMAL'] > 0:
			return False
		# require 20% AF
		if variant.INFO['TUMOUR_AF'][0] <= 0.20:
			return False

	return True

def classify_by_model(args, checkpoints, time_str):
	""" read VCF into a dataframe, classify using a model """

	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix, args.tumour_only)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	loaded_model = pickle.load(open(args.custom_model, "rb"))
	helper.time_function("Loaded classification model", checkpoints, time_str)

	prediction_dict = pool_predict(data_matrix, train.FEATURES_TO_DROP+["ID"], loaded_model, args.threads, args.confidence)

	helper.time_function("Performed prediction", checkpoints, time_str)
	class_dictionary = {
		-2: 'PREDICTED_BOTH', # mondrian
		-1: 'PREDICTED_NEITHER', # mondrian
		0: 'PREDICTED_NOISE',
		1: 'PREDICTED_SOMATIC',
		2: 'PREDICTED_GERMLINE',
		3: 'RESCUED_SOMATIC',
		4: 'REJECTED_LOW_SUPPORT',
		5: 'REJECTED_LOW_AF'
	}
	variant_classes = prediction_dict

	if args.confidence is not None:
		# convert the raw mondrian predictions into classes listed above
		for variant, raw_prediction in prediction_dict.items():
			# convert the raw predictions (inclu. mondrian)
			if type(raw_prediction) is list and len(raw_prediction) == 1:
				prediction_dict[variant] = raw_prediction[0]
			if type(raw_prediction) is list and len(raw_prediction) == 2:
				prediction_dict[variant] = -2
			if type(raw_prediction) is list and len(raw_prediction) == 0:
				prediction_dict[variant] = -1

	# assign the classes of variants in input vcf
	input_vcf = cyvcf2.VCF(args.vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		variant_mate_id = variant.INFO.get('MATEID', None)
		variant_class = prediction_dict.get(variant_id, None)
		# if no mate present, (as in ins, sbnds) set mate prediction to self prediction
		variant_mate_class = variant_classes.get(variant.INFO.get('MATEID', None), variant_class)
		if variant_class == 1 or variant_mate_class == 1:
			# check against hard filters that override prediction
			if not args.tumour_only:
				# filters can use NORMAL fields in hard filters
				if variant.INFO['TUMOUR_READ_SUPPORT'] < args.min_support:
					variant_classes[variant_id] = 4
					if variant_mate_id: # reject mate as well
						variant_classes[variant_mate_id] = 4
				elif not args.tumour_only and (variant.INFO['TUMOUR_READ_SUPPORT'] < variant.INFO['NORMAL_READ_SUPPORT']):
					variant_classes[variant_id] = 4
					if variant_mate_id: # reject mate as well
						variant_classes[variant_mate_id] = 4
				elif variant.INFO['TUMOUR_AF'][0] < args.min_af and variant.INFO['TUMOUR_AF'][1] < args.min_af:
					# determine if tumour-amplified
					avg_tumour_dp = mean([variant.INFO[dp] if isinstance(variant.INFO[dp], float) else variant.INFO[dp][0] for dp in ['TUMOUR_DP_BEFORE', 'TUMOUR_DP_AT', 'TUMOUR_DP_AFTER']])
					avg_normal_dp = mean([variant.INFO[dp] if isinstance(variant.INFO[dp], float) else variant.INFO[dp][0] for dp in ['NORMAL_DP_BEFORE', 'NORMAL_DP_AT', 'NORMAL_DP_AFTER']])
					if avg_tumour_dp <= avg_normal_dp*5:
						# no tumour amplification, low af holds
						variant_classes[variant_id] = 5
						if variant_mate_id: # reject mate as well
							variant_classes[variant_mate_id] = 5
				else:
					# both pass
					variant_classes[variant_id] = 1
					if variant_mate_id: # accept mate as well
						variant_classes[variant_mate_id] = 1
			else:
				# filters cannot use NORMAL fields in hard filters
				if variant.INFO['TUMOUR_READ_SUPPORT'] < args.min_support:
					variant_classes[variant_id] = 4
					if variant_mate_id: # reject mate as well
						variant_classes[variant_mate_id] = 4
				elif variant.INFO['TUMOUR_AF'][0] < args.min_af:
					variant_classes[variant_id] = 4
					if variant_mate_id: # reject mate as well
						variant_classes[variant_mate_id] = 4
				else:
					# both pass
					variant_classes[variant_id] = 1
					if variant_mate_id: # accept mate as well
						variant_classes[variant_mate_id] = 1
		elif variant_class == 2 and variant_mate_class == 2:
			# check against hard filters that override prediction
			if not args.tumour_only and (variant.INFO['NORMAL_READ_SUPPORT'] < args.min_support):
				variant_classes[variant_id] = 4
				if variant_mate_id: # reject mate as well
					variant_classes[variant_mate_id] = 4
			elif not args.tumour_only and (variant.INFO['NORMAL_AF'] < args.min_af):
				variant_classes[variant_id] = 5
				if variant_mate_id: # reject mate as well
					variant_classes[variant_mate_id] = 5
			else:
				# keep the true prediction
				variant_classes[variant_id] = variant_class
		elif variant_class == 0 or variant_mate_class == 0:
			# either edge was predicted as noise (and neither edge somatic or germline since failed earlier check)
			# update both edges to noise
			variant_classes[variant_id] = 0
			if variant_mate_id: # reject mate as well
				variant_classes[variant_mate_id] = 0
		elif (variant_class == -1 and variant_mate_class == -2) or (variant_class == -2 and variant_mate_class == -1):
			# tricky case - one edge predicted both and the other predicted neither.
			# update both edges to "BOTH"
			variant_classes[variant_id] = -2
			if variant_mate_id:
				variant_classes[variant_mate_id] = -2

	# all variants and mates have final classes, output files
	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str(f'Variant class prediction from model {args.custom_model} (after filtering for low AF or support)')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	desc_string = str(f'Variant filtered as likely noise by model: {args.custom_model}')
	input_vcf.add_filter_to_header({
		'ID': 'PREDICTED_NOISE',
		'Description': desc_string
	})
	if args.confidence is not None:
		desc_string = str('Variant predicted as neither noise nor somatic')
		input_vcf.add_filter_to_header({
			'ID': 'PREDICTED_NEITHER',
			'Description': desc_string
		})
		desc_string = str('Variant predicted as both noise and somatic')
		input_vcf.add_filter_to_header({
			'ID': 'PREDICTED_BOTH',
			'Description': desc_string
		})
	desc_string = str(f'Variant classified by model but removed due to low support (< {args.min_support} supporting reads)')
	input_vcf.add_filter_to_header({
		'ID': 'REJECTED_LOW_SUPPORT',
		'Description': desc_string
	})
	desc_string = str(f'Variant classified by model but removed due to AF < {args.min_af}')
	input_vcf.add_filter_to_header({
		'ID': 'REJECTED_LOW_AF',
		'Description': desc_string
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.somatic_output:
		somatic_vcf = cyvcf2.Writer(args.somatic_output, input_vcf)
	if args.germline_output:
		germline_vcf = cyvcf2.Writer(args.germline_output, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		variant_class = class_dictionary[variant_classes[variant_id]]
		variant.INFO['CLASS'] = variant_class
		# update the filter column here
		if variant_class not in ['PREDICTED_SOMATIC', 'PREDICTED_GERMLINE', 'RESCUED_SOMATIC']:
			variant.FILTER = variant_class
		# write germline and somatic files
		if args.somatic_output and variant_class in ['PREDICTED_SOMATIC', 'RESCUED_SOMATIC']:
			# add to the somatic only VCF
			somatic_vcf.write_record(variant)
		if args.germline_output and variant_class == 'PREDICTED_GERMLINE':
			germline_vcf.write_record(variant)
		out_vcf.write_record(variant)
	out_vcf.close()
	if args.somatic_output:
		somatic_vcf.close()
	if args.germline_output:
		germline_vcf.close()
	input_vcf.close()
	helper.time_function("Output classified VCF", checkpoints, time_str)
	if args.somatic_output:
		write_somatic_bedpe(args.somatic_output)
	helper.time_function("Output classified BEDPE", checkpoints, time_str)

	return

def write_somatic_bedpe(somatic_output):
	""" given a somatic vcf, wriite a somatic bedpe """

	# extract name for bedpe
	somatic_bedpe_file = f'{os.path.splitext(somatic_output)[0]}.bedpe'
	bedpe_string = ''
	for variant in cyvcf2.VCF(somatic_output):
		chrom = variant.CHROM
		start = variant.start + 1
		alt = variant.ALT[0]
		bp_id = variant.ID
		# print SV based on first edge, ignore second edge
		if bp_id.endswith('_1'):
			if alt == "<INS>" or '.' in alt:
				# ins or sbnd, no second edge
				alt_chr = chrom
				alt_pos = start + 1 # add one base for second pos
			else:
				# parse the alt column for chr/pos
				split_0, split_1 = alt.split(":")
				chr_match = re.search(r'[^[\]]*$', split_0)
				pos_match = re.search(r'^[^\[\]]*', split_1)
				alt_chr = chr_match.group(0)
				alt_pos = pos_match.group(0)

			# now get the length, support, and orientation from the INFO column
			support_string = ''
			for s in ['TUMOUR', 'NORMAL']:
				support = variant.INFO.get(s+"_READ_SUPPORT")
				support_string+=s+"_"+str(support)+"/" if support else ''
			support_string = support_string.rstrip("/")
			svlen = variant.INFO.get('SVLEN')
			orientation = variant.INFO.get('BP_NOTATION')
			sv_id = re.sub('_1$', '', bp_id) # remove pair identifier
			bedpe_string+=f'{chrom}\t{start}\t{start}\t{alt_chr}\t{alt_pos}\t{alt_pos}\t{sv_id}|{svlen}bp|{support_string}|{orientation}\n'

	with open(somatic_bedpe_file, 'w') as output:
		output.write(bedpe_string)


if __name__ == "__main__":
	print("Classification functions for SAVANA")
