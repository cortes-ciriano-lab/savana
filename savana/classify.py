"""
SAVANA strucural variant caller for long-read data - train sub-command
Created: 19/04/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import os

import pandas as pd
import numpy as np
import cyvcf2
import pickle
import json

from multiprocessing import Pool
from functools import partial
from collections import ChainMap

import savana.train as train
from savana.breakpoints import *
from savana.clusters import *

def get_info_fields(vcf):
	""" initialize the data header from INFO fields of sample VCF """
	fields = []
	for rec in vcf.header_iter():
		if rec['HeaderType'] == "INFO":
			fields.append(rec.info()['ID'])

	return fields

def pool_predict(data_matrix, features_to_drop, model, threads):
	""" """
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
	""" given feautres and a model, return the predictions """
	features = data_matrix.drop(features_to_drop, axis=1, errors='ignore')
	ids = data_matrix['ID']
	predicted = model.predict(features)
	prediction_dict = {}
	for i, value in enumerate(predicted):
		prediction_dict[ids.iloc[i]] = value

	return prediction_dict

def read_vcf(vcf_path):
	""" read VCF into a dataframe """
	data = []
	input_vcf = cyvcf2.VCF(vcf_path)
	header = ['ID'] + get_info_fields(input_vcf)
	for variant in input_vcf:
		row = [variant.ID]
		for field in header[1:]:
			row.append(variant.INFO.get(field))
		data.append(row)
	return pd.DataFrame(data, columns = header)

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
		except Exception as e:
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
		except Exception as e:
			# in case of None for any stat
			return False

def classify_legacy(args, checkpoints, time_str):
	""" classify using legacy lenient/strict filters """
	data_matrix = read_vcf(args.vcf)
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
	#read VCF into a dataframe, classify using a parameter JSON
	data_matrix = read_vcf(args.vcf)
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

def classify_by_model(args, checkpoints, time_str):
	""" read VCF into a dataframe, classify using a model """

	data_matrix = read_vcf(args.vcf)
	data_matrix = train.format_data(data_matrix)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	loaded_model = pickle.load(open(args.model, "rb"))
	helper.time_function("Loaded classification model", checkpoints, time_str)

	features_to_drop = [
		'ID', 'LABEL','MATEID','ORIGINATING_CLUSTER','END_CLUSTER',
		'TUMOUR_DP', 'NORMAL_DP', 'BP_NOTATION', 'SVTYPE', 'SVLEN',
		'ORIGIN_EVENT_SIZE_MEAN', 'ORIGIN_EVENT_SIZE_MEDIAN',
		'END_EVENT_SIZE_MEAN', 'END_EVENT_SIZE_MEDIAN'
		]
	prediction_dict = pool_predict(data_matrix, features_to_drop, loaded_model, 20)
	helper.time_function("Performed prediction", checkpoints, time_str)
	# output vcf using modified input_vcf as template
	input_vcf = cyvcf2.VCF(args.vcf)
	desc_string = str(f'Variant class prediction from model {args.model}')
	input_vcf.add_info_to_header({
		'ID': 'CLASS',
		'Type': 'String',
		'Number': '1',
		'Description': desc_string
	})
	desc_string = str(f'Variant filtered as likely noise by model: {args.model}')
	input_vcf.add_filter_to_header({
		'ID': 'LIKELY_NOISE',
		'Description': desc_string
	})
	out_vcf = cyvcf2.Writer(args.output, input_vcf)
	if args.somatic_output:
		somatic_vcf = cyvcf2.Writer(args.somatic_output, input_vcf)
	for variant in input_vcf:
		variant_id = variant.ID
		variant_prediction = prediction_dict.get(variant_id, None)
		if variant_prediction == 1:
			# PREDICTED SOMATIC BY MODEL
			# perform sanity checks
			if variant.INFO['TUMOUR_SUPPORT'] < 3:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
			elif variant.INFO['TUMOUR_SUPPORT'] < variant.INFO['NORMAL_SUPPORT']:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
			elif variant.INFO['NORMAL_SUPPORT']/variant.INFO['TUMOUR_SUPPORT'] > 0.1:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
			else:
				# update the INFO field if all sanity checks pass
				variant.INFO['CLASS'] = 'PREDICTED_SOMATIC'
				if args.somatic_output:
					# add to the somatic only VCF
					somatic_vcf.write_record(variant)
		elif variant_prediction == 2:
			# PREDICTED GERMLINE By MODEL
			# perform sanity checks
			if (variant.INFO['NORMAL_SUPPORT']+variant.INFO['TUMOUR_SUPPORT']) <= 3:
				variant.INFO['CLASS'] = 'PREDICTED_NOISE'
			else:
				# update the INFO field if sanity check passes
				variant.INFO['CLASS'] = 'PREDICTED_GERMLINE'
		else:
			# update the FILTER column
			variant.FILTER = 'LIKELY_NOISE'
			variant.INFO['CLASS'] = 'PREDICTED_NOISE'
		out_vcf.write_record(variant)
	out_vcf.close()
	if args.somatic_output:
		# add to the somatic only VCF
		somatic_vcf.close()
	input_vcf.close()
	helper.time_function("Output classified VCF", checkpoints, time_str)

	return

if __name__ == "__main__":
	print("Classification functions for SAVANA")
