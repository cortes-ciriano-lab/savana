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

import pickle
import json
import pandas as pd
import numpy as np
import cyvcf2

import savana.train as train
from savana.breakpoints import *
from savana.clusters import *

def pool_predict(data_matrix, features_to_drop, model, threads):
	""" split the prediction across multiprocessing Pool """
	data_matrices = np.array_split(data_matrix, threads)
	pool = Pool(threads)
	results = pool.map(partial(
		predict,
		features_to_drop=features_to_drop,
		model=model
	), data_matrices)
	pool.close()
	pool.join()

	return dict(ChainMap(*results))

def predict(data_matrix, features_to_drop, model):
	""" given features and a model, return the predictions """
	ids = data_matrix['ID']
	features = data_matrix.drop(features_to_drop, axis=1, errors='ignore')
	# reorder the columns based on model order
	feature_order = model.feature_names_in_
	features = features[feature_order]
	# perform prediction using model
	predicted = model.predict(features)
	prediction_dict = {}
	for i, value in enumerate(predicted):
		prediction_dict[ids.iloc[i]] = value

	return prediction_dict

def legacy_pass_strict(variant, event_heuristic):
	""" apply legacy filter for strict thresholds """
	if variant['NORMAL_SUPPORT'] > 0:
		return False
	if variant['TUMOUR_SUPPORT'] > 7:
		origin_uncertainty = (variant['ORIGIN_STARTS_STD_DEV']+1) * (variant['ORIGIN_EVENT_SIZE_STD_DEV']+1)
		end_uncertainty = (variant['END_STARTS_STD_DEV']+1) * (variant['END_EVENT_SIZE_STD_DEV']+1)
		try:
			if variant['TUMOUR_SUPPORT'] > 12 and origin_uncertainty <= 15:
				if event_heuristic <= 0.025:
					return True
			elif variant['TUMOUR_SUPPORT'] > 12 and end_uncertainty <= 30:
				return True
			elif end_uncertainty <= 10:
				return True
		except Exception as _:
			# in case of None
			return False

def legacy_pass_lenient(variant, event_heuristic):
	""" apply legacy filter for lenient thresholds """
	if variant['TUMOUR_SUPPORT'] == 0:
		return False
	elif variant['NORMAL_SUPPORT']/variant['TUMOUR_SUPPORT'] < 0.1:
		try:
			if variant['ORIGIN_STARTS_STD_DEV'] < 150 and event_heuristic < 3:
				if variant['BP_NOTATION'] == "<INS>" and variant['TUMOUR_SUPPORT'] > 25:
					return True
				elif variant['BP_NOTATION'] != "<INS>" and variant['TUMOUR_SUPPORT'] > 5:
					return True
		except Exception as _:
			# in case of None for any stat
			return False

def classify_legacy(args, checkpoints, time_str):
	""" classify using legacy lenient/strict filters """
	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix)
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
	data_matrix = train.format_data(data_matrix)
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
	if variant['NORMAL_SUPPORT'] != 0:
		return False
	if variant['TUMOUR_SUPPORT'] < min_support:
		return False
	if variant['ORIGIN_STARTS_STD_DEV'] > 50.0 or variant['END_STARTS_STD_DEV'] > 50.0:
		return False
	if variant['ORIGIN_MAPQ_MEAN'] < 40.0 or variant['END_MAPQ_MEAN'] < 40.0:
		return False
	if variant['ORIGIN_EVENT_SIZE_STD_DEV'] > 60.0 or variant['END_EVENT_SIZE_STD_DEV'] > 60.0:
		return False
	if variant.INFO['TUMOUR_AF'][0] < min_af:
		return False
	return True

def classify_pacbio(args, checkpoints, time_str):
	""" classify using PacBio filters (manually defined) """
	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	somatic_ids = {}
	for _, variant in data_matrix.iterrows():
		if pacbio_pass_somatic(variant, args.min_support, args.min_af):
			somatic_ids[variant['ID']] = True

	if args.predict_germline:
		print('Germline PacBio filters not yet implemented.')

	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str('Variant class prediction using PacBio filters')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	# create filenames based on output
	if ".vcf" in args.output:
		somatic_vcf_filename = args.output.replace(".vcf", ".somatic.vcf")
	else:
		somatic_vcf_filename = args.output+".somatic.vcf"
	somatic_vcf = cyvcf2.Writer(somatic_vcf_filename, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		passed_somatic = somatic_ids.get(variant_id, None)
		if passed_somatic:
			# update the INFO field if all sanity checks pass
			variant.INFO['CLASS'] = 'SOMATIC'
			# add to the somatic only VCF
			somatic_vcf.write_record(variant)
	somatic_vcf.close()
	input_vcf.close()

	helper.time_function("Output PacBio somatic VCF", checkpoints, time_str)

	return

def pass_rescue(variant):
	""" rescue variants predicted noise by model if they pass set of filters """
	if variant.INFO['NORMAL_SUPPORT'] != 0:
		return False
	if variant.INFO['TUMOUR_SUPPORT'] <= 6:
		return False
	if variant.INFO['TUMOUR_AF'][0] <= 0.10:
		return False
	if variant.INFO['CLUSTERED_READS_NORMAL'] > 1:
		return False
	if variant.INFO['ORIGIN_MAPQ_MEAN'] < 15.0:
		return False
	return True

def classify_by_model(args, checkpoints, time_str):
	""" read VCF into a dataframe, classify using a model """

	header, data = train.read_vcf(args.vcf)
	data_matrix = pd.DataFrame(data, columns=header)
	data_matrix = train.format_data(data_matrix)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	loaded_model = pickle.load(open(args.custom_model, "rb"))
	helper.time_function("Loaded classification model", checkpoints, time_str)

	prediction_dict = pool_predict(data_matrix, train.FEATURES_TO_DROP+["ID"], loaded_model, args.threads)

	helper.time_function("Performed prediction", checkpoints, time_str)
	# output vcf using modified input_vcf as template
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
		'ID': 'LIKELY_NOISE',
		'Description': desc_string
	})
	desc_string = str(f'Variant classified by model but removed due to low support (< {args.min_support} supporting reads)')
	input_vcf.add_filter_to_header({
		'ID': 'LOW_SUPPORT',
		'Description': desc_string
	})
	desc_string = str(f'Variant classified by model but removed due to AF < {args.min_af}')
	input_vcf.add_filter_to_header({
		'ID': 'LOW_AF',
		'Description': desc_string
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.somatic_output:
		somatic_vcf = cyvcf2.Writer(args.somatic_output, input_vcf)
	if args.germline_output:
		germline_vcf = cyvcf2.Writer(args.germline_output, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		variant_prediction = prediction_dict.get(variant_id, None)
		# if no mate present, (as in ins, sbnds) set mate prediction to self prediction
		variant_mate_prediction = prediction_dict.get(variant.INFO.get('MATEID', None), variant_prediction)
		if variant_prediction == 1 or variant_mate_prediction == 1:
			# AT LEAST ONE EDGE PREDICTED SOMATIC BY MODEL
			# perform sanity checks
			if variant.INFO['TUMOUR_SUPPORT'] < args.min_support:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
				variant.FILTER = 'LOW_SUPPORT'
				prediction_dict[variant_id] = 0 # update the dictionary as well for mate
			elif variant.INFO['TUMOUR_SUPPORT'] < variant.INFO['NORMAL_SUPPORT']:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
				variant.FILTER = 'LOW_SUPPORT'
				prediction_dict[variant_id] = 0 # update the dictionary as well for mate
			elif variant.INFO['TUMOUR_AF'][0] < args.min_af and variant.INFO['TUMOUR_AF'][1] < args.min_af:
				# determine if tumour-amplified
				avg_tumour_dp = mean([variant.INFO[dp] if isinstance(variant.INFO[dp], float) else variant.INFO[dp][0] for dp in ['TUMOUR_DP_BEFORE', 'TUMOUR_DP_AT', 'TUMOUR_DP_AFTER']])
				avg_normal_dp = mean([variant.INFO[dp] if isinstance(variant.INFO[dp], float) else variant.INFO[dp][0] for dp in ['NORMAL_DP_BEFORE', 'NORMAL_DP_AT', 'NORMAL_DP_AFTER']])
				if avg_tumour_dp <= avg_normal_dp*5:
					# no tumour amplification, low af holds
					variant.INFO['CLASS'] = 'PREDICTED_NOISE'
					variant.FILTER = 'LOW_AF'
					prediction_dict[variant_id] = 0 # update the dictionary as well for mate
				else:
					# tumour appears amplified, disregard that variant has low af
					variant.INFO['CLASS'] = 'PREDICTED_SOMATIC'
					if args.somatic_output:
						# add to the somatic only VCF
						somatic_vcf.write_record(variant)
			else:
				# update the INFO field if all sanity checks pass
				variant.INFO['CLASS'] = 'PREDICTED_SOMATIC'
				if args.somatic_output:
					# add to the somatic only VCF
					somatic_vcf.write_record(variant)
		elif variant_prediction == 2 and variant_mate_prediction == 2:
			# PREDICTED GERMLINE By MODEL
			# perform sanity checks
			if variant.INFO['NORMAL_SUPPORT'] < args.min_support:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
				variant.FILTER = 'LOW_SUPPORT'
				prediction_dict[variant_id] = 0 # update the dictionary as well for mate
			elif float(variant.INFO['NORMAL_AF']) < args.min_af:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
				variant.FILTER = 'LOW_AF'
				prediction_dict[variant_id] = 0 # update the dictionary as well for mate
			else:
				# update the INFO field if sanity check passes
				variant.INFO['CLASS'] = 'PREDICTED_GERMLINE'
				if args.germline_output:
					# add to the somatic only VCF
					germline_vcf.write_record(variant)
		else:
			# update the FILTER column
			variant.FILTER = 'LIKELY_NOISE'
			variant.INFO['CLASS'] = 'PREDICTED_NOISE'
		out_vcf.write_record(variant)
	out_vcf.close()
	if args.somatic_output:
		somatic_vcf.close()
	if args.germline_output:
		germline_vcf.close()
	input_vcf.close()
	helper.time_function("Output classified VCF", checkpoints, time_str)

	return

if __name__ == "__main__":
	print("Classification functions for SAVANA")
