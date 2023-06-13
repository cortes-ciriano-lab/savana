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

def classify_vcf(args, checkpoints, time_str):
	""" read VCF into a dataframe """
	data = []
	input_vcf = cyvcf2.VCF(args.vcf)
	header = ['ID'] + get_info_fields(input_vcf)
	for variant in input_vcf:
		row = [variant.ID]
		for field in header[1:]:
			row.append(variant.INFO.get(field))
		data.append(row)
	data_matrix = pd.DataFrame(data, columns = header)
	data_matrix = train.format_data(data_matrix)
	helper.time_function("Loaded raw breakpoints", checkpoints, time_str)

	loaded_model = pickle.load(open(args.model, "rb"))
	helper.time_function("Loaded classification model", checkpoints, time_str)

	features_to_drop = [
		'ID', 'LABEL','MATEID','ORIGINATING_CLUSTER','END_CLUSTER',
		'TUMOUR_DP', 'NORMAL_DP', 'BP_NOTATION', 'SVTYPE'
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
			# perform sanity checks on the model label
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
			# perform sanity checks on the model label
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
