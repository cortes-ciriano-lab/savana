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

from random import random

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import *
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint, uniform
from sklearn.utils.class_weight import *

# for plotting confusion matrix
import matplotlib.pyplot as plt

label_encoding = {
	'NOT_IN_COMPARISON': 0,
	'SOMATIC': 1,
	'FOUND_IN_BOTH': 2,
	'GERMLINE': 2
}
# inverted label_encoding
encoded_labels = {v: k for k, v in label_encoding.items()}

def format_data(data_matrix):
	""" parse columns, do conversions, one-hot-encoding """
	# split the DP tuples
	data_matrix[['TUMOUR_DP_BEFORE_0', 'TUMOUR_DP_BEFORE_1']] = data_matrix['TUMOUR_DP_BEFORE'].apply(pd.Series)
	data_matrix[['TUMOUR_DP_AT_0', 'TUMOUR_DP_AT_1']] = data_matrix['TUMOUR_DP_AT'].apply(pd.Series)
	data_matrix[['TUMOUR_DP_AFTER_0', 'TUMOUR_DP_AFTER_1']] = data_matrix['TUMOUR_DP_AFTER'].apply(pd.Series)
	data_matrix[['NORMAL_DP_BEFORE_0', 'NORMAL_DP_BEFORE_1']] = data_matrix['NORMAL_DP_BEFORE'].apply(pd.Series)
	data_matrix[['NORMAL_DP_AT_0', 'NORMAL_DP_AT_1']] = data_matrix['NORMAL_DP_AT'].apply(pd.Series)
	data_matrix[['NORMAL_DP_AFTER_0', 'NORMAL_DP_AFTER_1']] = data_matrix['NORMAL_DP_AFTER'].apply(pd.Series)
	# when nothing in second depth column (insertions), replace with value in first
	data_matrix['TUMOUR_DP_BEFORE_1'] = data_matrix['TUMOUR_DP_BEFORE_1'].fillna(data_matrix['TUMOUR_DP_BEFORE_0'])
	data_matrix['TUMOUR_DP_AT_1'] = data_matrix['TUMOUR_DP_AT_1'].fillna(data_matrix['TUMOUR_DP_AT_0'])
	data_matrix['TUMOUR_DP_AFTER_1'] = data_matrix['TUMOUR_DP_AFTER_1'].fillna(data_matrix['TUMOUR_DP_AFTER_0'])
	data_matrix['NORMAL_DP_BEFORE_1'] = data_matrix['NORMAL_DP_BEFORE_1'].fillna(data_matrix['NORMAL_DP_BEFORE_0'])
	data_matrix['NORMAL_DP_AT_1'] = data_matrix['NORMAL_DP_AT_1'].fillna(data_matrix['NORMAL_DP_AT_0'])
	data_matrix['NORMAL_DP_AFTER_1'] = data_matrix['NORMAL_DP_AFTER_1'].fillna(data_matrix['NORMAL_DP_AFTER_0'])
	data_matrix.replace([np.inf, -np.inf], -1, inplace=True)
	# create std_dev/mean_size ratio columns
	data_matrix['ORIGIN_STD_MEAN_RATIO'] = data_matrix['ORIGIN_STARTS_STD_DEV']/(data_matrix['ORIGIN_EVENT_SIZE_MEAN']+1.0)
	data_matrix['END_STD_MEAN_RATIO'] = data_matrix['END_STARTS_STD_DEV']/(data_matrix['END_EVENT_SIZE_MEAN']+1.0)
	data_matrix['ORIGIN_STD_MEAN_RATIO'] = data_matrix['ORIGIN_STD_MEAN_RATIO'].fillna(0)
	data_matrix['END_STD_MEAN_RATIO'] = data_matrix['END_STD_MEAN_RATIO'].fillna(0)
	# convert the SVTYPE to 0/1/2
	data_matrix['SVTYPE'] = data_matrix['SVTYPE'].map({'BND':0,'INS':1, 'SBND': 2})
	# convert the SOURCE to 0/1/2
	data_matrix['SOURCE'] = data_matrix['SOURCE'].map({'CIGAR':0,'SUPPLEMENTARY':1, 'CIGAR/SUPPLEMENTARY': 2})
	# one-hot-encoding of BP_NOTATION
	sv_type_one_hot = pd.get_dummies(data_matrix['BP_NOTATION'])
	# check to make sure all bp types are present
	for bp_type in ["++","+-","-+","--"]:
		if bp_type not in sv_type_one_hot:
			sv_type_one_hot[bp_type] = False
	data_matrix.drop('BP_NOTATION', axis=1)
	data_matrix = data_matrix.join(sv_type_one_hot)

	return data_matrix

def prepare_data(data_matrix, germline_class):
	""" add predictor and split into train/test """
	# reformat/parse columns
	data_matrix = format_data(data_matrix)
	# not enough ONT evidence to give these a 'SOMATIC' label
	# relabelling to avoid confusing the model
	condition = (data_matrix['LABEL'] == 'SOMATIC') & (data_matrix['TUMOUR_SUPPORT'] < 3)
	data_matrix.loc[condition, 'LABEL'] = 'NOT_IN_COMPARISON'
	condition = (data_matrix['LABEL'] == 'SOMATIC') & (data_matrix['TUMOUR_SUPPORT'] < data_matrix['NORMAL_SUPPORT'])
	data_matrix.loc[condition, 'LABEL'] = 'NOT_IN_COMPARISON'
	condition = (data_matrix['LABEL'] == 'SOMATIC') & (data_matrix['NORMAL_SUPPORT']/data_matrix['TUMOUR_SUPPORT'] > 0.1)
	data_matrix.loc[condition, 'LABEL'] = 'NOT_IN_COMPARISON'
	# similar threshold for 'GERMLINE'
	condition = (data_matrix['LABEL'] == 'GERMLINE') & (data_matrix['NORMAL_SUPPORT']+data_matrix['TUMOUR_SUPPORT'] < 3)
	data_matrix.loc[condition, 'LABEL'] = 'NOT_IN_COMPARISON'
	if germline_class:
		# encode the labels to 0/1/2 for FALSE/SOMATIC/GERMLINE
		data_matrix['LABEL'] = data_matrix['LABEL'].map(
			label_encoding
		)
	else:
		# encode the labels to 0/1 for FALSE/FOUND
		label_encoding_reduced = label_encoding
		label_encoding_reduced['GERMLINE'] = 0
		data_matrix['LABEL'] = data_matrix['LABEL'].map(
			label_encoding_reduced
		)

	# drop irrelevant/redundant columns (some have been encoded in a different format)
	features = data_matrix.drop([
		'LABEL','MATEID','ORIGINATING_CLUSTER','END_CLUSTER',
		'TUMOUR_DP_BEFORE', 'TUMOUR_DP_AT', 'TUMOUR_DP_AFTER',
		'NORMAL_DP_BEFORE', 'NORMAL_DP_AT', 'NORMAL_DP_AFTER',
		'BP_NOTATION', '<INS>', 'LABEL_VARIANT_ID',
		'ORIGIN_EVENT_SIZE_MEAN', 'ORIGIN_EVENT_SIZE_MEDIAN',
		'END_EVENT_SIZE_MEAN', 'END_EVENT_SIZE_MEDIAN', 'CLASS',
		'REPEAT', 'BLACKLIST', 'INS_PON', 'MICROSATELLITE'
		], axis=1, errors='ignore')
	target = data_matrix['LABEL']

	return features, target

def format_value_counts(value_counts):
	""" print this out prettier """
	values = value_counts.keys().tolist()
	counts = value_counts.tolist()
	for i in range(0,len(values)):
		print(f'{round(counts[i], 4)} - {encoded_labels[values[i]]}')


def fit_classifier(X, y, outdir, split, downsample, hyperparameter, germline_class):
	""" given the features (X) and target (y), split into test/train and fit the model """
	# split into train/test
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split)

	# randomly downsample the 'NOT_IN_COMPARISON'/0 by 'downsample' %
	train_matrix = pd.concat([X_train, y_train], axis=1)
	print(f'\nCOUNTS PRE-DOWNSAMPLING:')
	format_value_counts(train_matrix['LABEL'].value_counts(normalize=False))
	print(f'PROPORTIONS PRE-DOWNSAMPLING:')
	format_value_counts(train_matrix['LABEL'].value_counts(normalize=True))
	# perform downsampling
	print(f'\n>DOWNSAMPLING MAJORITY CLASS BY {downsample}')
	mask = train_matrix['LABEL'] == 0
	train_matrix = pd.concat([train_matrix[mask].sample(frac=(1.0-downsample)), train_matrix[~mask]])
	print(f'\nCOUNTS POST-DOWNSAMPLING:')
	format_value_counts(train_matrix['LABEL'].value_counts(normalize=False))
	print(f'PROPORTIONS POST-DOWNSAMPLING:')
	format_value_counts(train_matrix['LABEL'].value_counts(normalize=True))

	# put back into X_train and y_train
	X_train = train_matrix.drop('LABEL', axis=1)
	y_train = train_matrix['LABEL']

	# fit the random forest
	random_forest = None
	if hyperparameter:
		# hyper-parameter testing
		param_dist = {'n_estimators': randint(400,500), 'max_depth': randint(10,20)}
		rf = RandomForestClassifier(n_jobs=16)
		rand_search = RandomizedSearchCV(rf,
					param_distributions=param_dist,
					n_iter=5,
					cv=5,
					scoring='f1_macro')
		rand_search.fit(X_train, y_train)
		print('Using Best Hyperparameters:')
		print(', '.join(f'{key}: {value}' for key, value in rand_search.best_params_.items()))
		# store the best model
		random_forest = rand_search.best_estimator_
	else:
		# use default pre-deteremined hyperparameters
		random_forest = RandomForestClassifier(max_depth=20, n_estimators=400, n_jobs=16)
		random_forest.fit(X_train, y_train)
	y_pred = random_forest.predict(X_test)

	# FIRST STATS METHOD
	average_method = 'weighted'
	stats = {
		'Precision': precision_score(y_test, y_pred, average=average_method),
		'Recall': recall_score(y_test, y_pred, average=average_method),
		'F-score': f1_score(y_test, y_pred, average=average_method)
	}
	print('\n- Weighted Stats -')
	for stat, value in stats.items():
		print(f'{stat}: {round(value, 3)}')
	# SECOND STATS METHOD
	average_method = 'macro'
	stats = {
		'Precision': precision_score(y_test, y_pred, average=average_method),
		'Recall': recall_score(y_test, y_pred, average=average_method),
		'F-score': f1_score(y_test, y_pred, average=average_method)
	}
	print('\n- Macro Stats -')
	for stat, value in stats.items():
		print(f'{stat}: {round(value, 3)}')

	# print out the feature importances
	print('\nSorted Feature Importances:')
	feature_importances = pd.Series(random_forest.feature_importances_, index=X_train.columns).sort_values(ascending=False)
	print(feature_importances.to_string())

	print('\nConfusion Matrix::')
	labels = np.array(['FALSE', 'SOMATIC', 'GERMLINE']) if germline_class else np.array(['FALSE', 'TRUE'])
	cm=confusion_matrix(y_test, y_pred)
	disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
	disp.plot()
	plt.savefig(os.path.join(outdir, 'confusion_matrix.png'))
	plt.show()

	report = classification_report(y_test, y_pred)
	print(report)
	with open(os.path.join(outdir, 'model_stats.txt'), 'w') as f:
		f.write(report)

	# export y_test and pred
	test_matrix = X_test
	test_matrix['TRUE_LABEL'] = y_test
	test_matrix['PREDICTED_LABEL'] = y_pred
	# only interested in the FP/FN
	test_matrix = test_matrix[test_matrix['TRUE_LABEL'] != test_matrix['PREDICTED_LABEL']]
	test_matrix.to_csv(os.path.join(outdir, 'test_set_incorrect.tsv'), sep="\t")
	# only interested in the TP
	test_matrix = test_matrix[test_matrix['TRUE_LABEL'] == test_matrix['PREDICTED_LABEL']]
	test_matrix.to_csv(os.path.join(outdir, 'test_set_correct.tsv'), sep="\t")

	return random_forest

def save_model(args, model, outdir):
	""" save model and its info to pickle file """
	model_path = os.path.join(outdir, 'random_forest_model.pkl')
	print(f'Saving classifier to {model_path}')
	model.set_params(n_jobs=1) # reset jobs before saving
	pickle.dump(model, open(model_path, "wb"))
	# also save arguments used to create model
	cmd_string = 'savana train'
	for arg, value in vars(args).items():
		if value and arg != "func":
			cmd_string+=f' --{arg} {value}'
	cmd_string+='"'
	with open(os.path.join(outdir, 'model_arguments.txt'), 'w') as f:
		f.write(cmd_string)

def load_matrix(args):
	"""  read in dataframe from pickle file """
	print(f'Loading data matrix from {args.load_matrix}')
	df = pd.read_pickle(args.load_matrix)
	return df

def get_info_fields(vcf):
	""" initialize the data header from INFO fields of sample VCF """
	fields = []
	for rec in vcf.header_iter():
		if rec['HeaderType'] == "INFO":
			fields.append(rec.info()['ID'])

	return fields

def read_vcfs(args):
	""" given the folder of labelled input VCFs, return an output dataframe """
	header = []
	data = []
	for root, _, file_names in os.walk(args.vcfs):
		for file in file_names:
			f = os.path.join(root, file)
			if os.path.isfile(f) and file.endswith('.vcf'):
				print(f'Loading {f} into matrix')
				input_vcf = cyvcf2.VCF(f)
				if not header:
					header = get_info_fields(input_vcf)
				for variant in input_vcf:
					row = []
					for field in header:
						row.append(variant.INFO.get(field))
					data.append(row)
	df = pd.DataFrame(data, columns = header)
	if args.save_matrix:
		print(f'Saving data matrix to pickle file {args.save_matrix}')
		df.to_pickle(args.save_matrix)
	return df

if __name__ == "__main__":
	print("Training functions for SAVANA")
