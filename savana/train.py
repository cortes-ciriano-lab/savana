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
	for bp_type in ["++","+-","-+","--","<INS>", "<SBND>"]:
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

def cross_conformal_classifier(X, y, outdir, split, downsample):
	""" """
	# perform 80/20 split
	X_train_with_ids, X_reserved_test_with_ids, y_train, y_reserved_test = train_test_split(X, y, test_size=split)

	# remove IDs
	X_train = X_train_with_ids.iloc[: , 2:]
	X_reserved_test = X_reserved_test_with_ids.copy(deep=True) # prevent ids from being dropped
	X_reserved_test = X_reserved_test.iloc[: , 2:]

	#TODO try downsampling?

	print(f'\nDistributions of Training Data:')
	format_value_counts(y_train.value_counts(normalize=False))
	format_value_counts(y_train.value_counts(normalize=True))
	print(f'\nDistributions of Testing Data:')
	format_value_counts(y_reserved_test.value_counts(normalize=False))
	format_value_counts(y_reserved_test.value_counts(normalize=True))

	# split the data into 5 folds for cross-conformal prediction
	kf = StratifiedKFold(n_splits=5)
	class_nconform_scores = {}
	fold = 1
	print('')
	for train_index, test_index in kf.split(X_train, y_train):
		print(f'Training on Fold {fold}')
		X_k_train, X_k_test = X_train.iloc[train_index], X_train.iloc[test_index]
		y_k_train, y_k_test = y_train.iloc[train_index], y_train.iloc[test_index]
		# train on k-fold split
		random_forest = RandomForestClassifier(max_depth=20, n_estimators=400, n_jobs=16)
		random_forest.fit(X_k_train, y_k_train)
		# test on k-fold split
		probabilities = random_forest.predict_proba(X_k_test)
		fold_nconform_scores = {}
		for i, prob in enumerate(probabilities):
			true_class = y_k_test.iloc[i]
			nonconform_score = prob[true_class]
			if true_class not in fold_nconform_scores:
				# initialize
				fold_nconform_scores[true_class] = np.array([])
			fold_nconform_scores[true_class] = np.append(fold_nconform_scores[true_class], nonconform_score)
		# once all the non-conform scores collected, append to full nonconfom scores list
		for key in fold_nconform_scores.keys():
			if key not in class_nconform_scores:
				# initialize
				#class_nconform_scores[key] = np.array([])
				# alt - for averaging the p-values
				pass
				#class_nconform_scores.setdefault(key, {})[fold] = np.array([])
				"""
				if fold not in class_nconform_scores[key]:
					class_nconform_scores[key][fold] = np.array([])
				"""
			class_nconform_scores.setdefault(key, {})[fold] = np.array([])
			#class_nconform_scores[key] = np.append(class_nconform_scores[key], fold_nconform_scores[key])
			class_nconform_scores[key][fold] = np.append(class_nconform_scores[key][fold], fold_nconform_scores[key])
		fold += 1

	# sort class non-conformity scores
	for key in class_nconform_scores.keys():
		#class_nconform_scores[key] = np.sort(class_nconform_scores[key])
		for fold in class_nconform_scores[key].keys():
			class_nconform_scores[key][fold] = np.sort(class_nconform_scores[key][fold])

	# now fit random forest on full training set
	random_forest = RandomForestClassifier(max_depth=20, n_estimators=400, n_jobs=16)
	random_forest.fit(X_train, y_train)
	# test on reserved test set
	test_probabilities = random_forest.predict_proba(X_reserved_test)
	confidence = 0.80
	epsilon = 1.0 - confidence
	# for testing, write out the Mondrian arrays of non-conformity scores
	# requires concatenating and sorting all the folds
	sorted_class_scores = {}
	for key, fold_nconform_scores in class_nconform_scores.items():
		sorted_class_scores[key] = np.array([])
		for fold, scores in fold_nconform_scores.items():
			sorted_class_scores[key] = np.append(sorted_class_scores[key], scores)
		# sort once all folds appended
		sorted_class_scores[key] = np.sort(sorted_class_scores[key])
		# then write out to file
		np.savetxt(f'{outdir}/{key}_mondrian_dist.txt', sorted_class_scores[key], delimiter='\t', header='probability', comments='')

	#for key, scores in class_nconform_scores.items():
	#	np.savetxt(f'{outdir}/{key}.txt', scores, delimiter='\t', header='probability', comments='')
	print(f'\nUsing Epsilon of {round(epsilon, 3)}')
	# initialize blank arrays for 0, 1 and 'truth'
	p_class_assignments = {k: np.array([]) for k in class_nconform_scores.keys()}
	p_class_assignments['truth'] = np.array([])
	for i, prob in enumerate(test_probabilities):
		p_class_assignments['truth'] = np.append(p_class_assignments['truth'], y_reserved_test.iloc[i])
		for key, class_scores in class_nconform_scores.items():
			#fraction = len(scores[scores <= prob[key]])/len(scores)
			fractions = []
			for fold, fold_scores in class_scores.items():
				fractions.append(len(fold_scores[fold_scores <= prob[key]])/len(fold_scores))
			fraction = np.mean(fractions)
			if fraction >= epsilon:
				p_class_assignments[key] = np.append(p_class_assignments[key], 1)
			else:
				p_class_assignments[key] = np.append(p_class_assignments[key], 0)

	count_neither = 0
	count_both = 0
	tp_somatic, tn_somatic, fp_somatic, fn_somatic = 0, 0, 0, 0
	tp_noise, tn_noise, fp_noise, fn_noise = 0, 0, 0, 0
	for i in range(len(test_probabilities)):
		if p_class_assignments[0][i] == 0 and p_class_assignments[1][i] == 0:
			count_neither += 1
		elif p_class_assignments[0][i] == 1 and p_class_assignments[1][i] == 1:
			count_both += 1
		# SOMATIC CONFUSION MATRIX VALUES
		if p_class_assignments[1][i] == 1 and p_class_assignments['truth'][i] == 1:
			# somatic = true and class = somatic
			tp_somatic += 1
		if p_class_assignments[1][i] == 0 and p_class_assignments['truth'][i] == 0:
			# somatic = false and class = noise
			tn_somatic += 1
		if p_class_assignments[1][i] == 1 and p_class_assignments['truth'][i] == 0:
			# somatic = true but class = noise
			fp_somatic += 1
		if p_class_assignments[1][i] == 0 and p_class_assignments['truth'][i] == 1:
			# somatic = false but class = somatic
			fn_somatic += 1
		# NOISE CONFUSION MATRIX VALUES
		if p_class_assignments[0][i] == 1 and p_class_assignments['truth'][i] == 0:
			# noise = true and class = noise
			tp_noise += 1
		if p_class_assignments[0][i] == 0 and p_class_assignments['truth'][i] == 1:
			# noise = false and class = somatic
			tn_noise += 1
		if p_class_assignments[0][i] == 1 and p_class_assignments['truth'][i] == 1:
			# noise = true but class = somatic
			fp_noise += 1
		if p_class_assignments[0][i] == 0 and p_class_assignments['truth'][i] == 0:
			# noise = false but class = noise
			fn_noise += 1

	print('\nResults using Cross-conformal prediction')
	print(f'Neither: {count_neither} | Both: {count_both} | Total: {len(test_probabilities)}')
	# SOMATIC SCORES
	precision = tp_somatic/(tp_somatic+fp_somatic)
	recall = tp_somatic/(tp_somatic+fn_somatic)
	f_measure = (2*precision*recall)/(precision+recall)
	print(f'Somatic Scores:')
	print(f'TP Somatic: {tp_somatic} | FN Somatic: {fn_somatic}')
	print(f'FP Somatic: {fp_somatic} | TN Somatic: {tn_somatic}')
	print(f'Precision: {round(precision,2)} | Recall: {round(recall,2)} | F-measure: {round(f_measure,2)}')
	# NOISE SCORES
	precision = tp_noise/(tp_noise+fp_noise)
	recall = tp_noise/(tp_noise+fn_noise)
	f_measure = (2*precision*recall)/(precision+recall)
	print(f'Noise Scores:')
	print(f'TP Noise: {tp_noise} | FN Noise: {fn_noise}')
	print(f'FP Noise: {fp_noise} | TN Noise: {tn_noise}')
	print(f'Precision: {round(precision,2)} | Recall: {round(recall,2)} | F-measure: {round(f_measure,2)}')

	groups = {}
	for group in ['both','null','noise', 'somatic']:
		groups[group] = [0]*2
	for i in range(len(test_probabilities)):
		truth = int(p_class_assignments['truth'][i])
		if p_class_assignments[0][i] == 1 and p_class_assignments[1][i] == 1:
			groups['both'][truth]+=1
		if p_class_assignments[0][i] == 0 and p_class_assignments[1][i] == 0:
			groups['null'][truth]+=1
		if p_class_assignments[0][i] == 1 and p_class_assignments[1][i] == 0:
			groups['noise'][truth]+=1
		if p_class_assignments[0][i] == 0 and p_class_assignments[1][i] == 1:
			groups['somatic'][truth]+=1

	with open(f'{outdir}/class_counts.txt', 'w') as fp:
		fp.write('\t'.join(['group','count_0', 'count_1']))
		fp.write('\n')
		for key, counts in groups.items():
			fp.write('\t'.join([str(i) for i in [key,counts[0],counts[1]]]))
			fp.write('\n')

	# USING REGULAR NON-CONFORMAL VOTE
	print('\nResults using standard prediction')
	voting_class_assignments = {k: np.array([]) for k in class_nconform_scores.keys()}
	voting_class_assignments['truth'] = np.array([])
	for i, prob in enumerate(test_probabilities):
		voting_class_assignments['truth'] = np.append(voting_class_assignments['truth'], y_reserved_test.iloc[i])
		if prob[0] > prob[1]:
			voting_class_assignments[0] = np.append(voting_class_assignments[0], 1)
			voting_class_assignments[1] = np.append(voting_class_assignments[1], 0)
		else:
			voting_class_assignments[1] = np.append(voting_class_assignments[1], 1)
			voting_class_assignments[0] = np.append(voting_class_assignments[0], 0)

	count_neither = 0
	count_both = 0
	tp_somatic, tn_somatic, fp_somatic, fn_somatic = 0, 0, 0, 0
	tp_noise, tn_noise, fp_noise, fn_noise = 0, 0, 0, 0
	for i in range(len(test_probabilities)):
		if voting_class_assignments[0][i] == 0 and voting_class_assignments[1][i] == 0:
			count_neither += 1
		elif voting_class_assignments[0][i] == 1 and voting_class_assignments[1][i] == 1:
			count_both += 1
		# SOMATIC CONFUSION MATRIX VALUES
		if voting_class_assignments[1][i] == 1 and voting_class_assignments['truth'][i] == 1:
			# somatic = true and class = somatic
			tp_somatic += 1
		if voting_class_assignments[1][i] == 0 and voting_class_assignments['truth'][i] == 0:
			# somatic = false and class = noise
			tn_somatic += 1
		if voting_class_assignments[1][i] == 1 and voting_class_assignments['truth'][i] == 0:
			# somatic = true but class = noise
			fp_somatic += 1
		if voting_class_assignments[1][i] == 0 and voting_class_assignments['truth'][i] == 1:
			# somatic = fasle but class = somatic
			fn_somatic += 1
		# NOISE CONFUSION MATRIX VALUES
		if voting_class_assignments[0][i] == 1 and voting_class_assignments['truth'][i] == 0:
			# noise = true and class = noise
			tp_noise += 1
		if voting_class_assignments[0][i] == 0 and voting_class_assignments['truth'][i] == 1:
			# noise = false and class = somatic
			tn_noise += 1
		if voting_class_assignments[0][i] == 1 and voting_class_assignments['truth'][i] == 1:
			# noise = true but class = somatic
			fp_noise += 1
		if voting_class_assignments[0][i] == 0 and voting_class_assignments['truth'][i] == 0:
			# noise = false but class = noise
			fn_noise += 1

	print(f'Neither: {count_neither} | Both: {count_both} | Total: {len(test_probabilities)}')
	# SOMATIC SCORES
	precision = tp_somatic/(tp_somatic+fp_somatic)
	recall = tp_somatic/(tp_somatic+fn_somatic)
	f_measure = (2*precision*recall)/(precision+recall)
	print(f'Somatic Scores:')
	print(f'TP Somatic: {tp_somatic} | FN Somatic: {fn_somatic}')
	print(f'FP Somatic: {fp_somatic} | TN Somatic: {tn_somatic}')
	print(f'Precision: {round(precision,2)} | Recall: {round(recall,2)} | F-measure: {round(f_measure,2)}')

	# NOISE SCORES
	precision = tp_noise/(tp_noise+fp_noise)
	recall = tp_noise/(tp_noise+fn_noise)
	f_measure = (2*precision*recall)/(precision+recall)
	print(f'Noise Scores:')
	print(f'TP Noise: {tp_noise} | FN Noise: {fn_noise}')
	print(f'FP Noise: {fp_noise} | TN Noise: {tn_noise}')
	print(f'Precision: {round(precision,2)} | Recall: {round(recall,2)} | F-measure: {round(f_measure,2)}')

	with open(f'{outdir}/class_counts_voted.txt', 'w') as fp:
		fp.write('\t'.join(['group','count_0', 'count_1']))
		fp.write('\n')
		for key, counts in groups.items():
			fp.write('\t'.join([str(i) for i in [key,counts[0],counts[1]]]))
			fp.write('\n')

	# EVALUATE
	print('\n- Precision-Recall Curve -')
	# now drop for the pr-curve # TODO remove this, it's unneeded i think
	#X_reserved_test_no_ids = X_reserved_test.drop(['SAMPLE','ID'], axis=1, errors='ignore')
	probs = random_forest.predict_proba(X_reserved_test)
	precision, recall, _ = precision_recall_curve(y_reserved_test, probs[:, 1])
	print(f'AUC: {round(auc(recall, precision),2)}')

	# print out the feature importances
	print('\nSorted Feature Importances:')
	feature_importances = pd.Series(random_forest.feature_importances_, index=X_train.drop(['SAMPLE','ID'], axis=1, errors='ignore').columns).sort_values(ascending=False)
	print(feature_importances.to_string())

	# export the data and predictions
	test_matrix = X_reserved_test_with_ids.copy(deep=True)
	y_pred = random_forest.predict(X_reserved_test)
	test_matrix['TRUE_LABEL'] = y_reserved_test.copy()
	test_matrix['STANDARD_PREDICTED_LABEL'] = y_pred
	test_matrix['CONFORMAL_PREDICT_NOISE'] = p_class_assignments[0]
	test_matrix['CONFORMAL_PREDICT_SOMATIC'] = p_class_assignments[1]
	try:
		test_matrix['PREDICT_SCORE'] = probs[:, 1]
	except Exception as e:
		print("Couldn't add prediction score")

	# just print all
	test_matrix.to_csv(os.path.join(outdir, f'test_set_all_random_forest.tsv'), sep="\t", index=False)

	# this step takes forever - do in R instead
	#test_cutoffs(random_forest, X_reserved_test_no_ids, y_reserved_test, outdir, 'random-forest')

	return random_forest

def test_cutoffs(model, X_test, y_test, outdir, model_type):
	""" """
	# TESTING WITH DIFFERENT CUTOFFS
	print(f'\nResults using different cutoffs for model: {model_type}')
	y_pred_proba = model.predict_proba(X_test)
	cutoffs = [0+(0.01*i) for i in range(0,101)]
	cutoff_results = {c: {'somatic_precision': None, 'somatic_recall': None, 'noise_precision': None, 'noise_recall': None} for c in cutoffs}
	for cutoff in cutoffs:
		cutoff_class_assignments = {k: np.array([]) for k in [0,1]}
		cutoff_class_assignments['truth'] = np.array([])
		for i, prob in enumerate(y_pred_proba):
			cutoff_class_assignments['truth'] = np.append(cutoff_class_assignments['truth'], y_test.iloc[i])
			if prob[0] > cutoff:
				cutoff_class_assignments[0] = np.append(cutoff_class_assignments[0], 1)
				cutoff_class_assignments[1] = np.append(cutoff_class_assignments[1], 0)
			else:
				cutoff_class_assignments[1] = np.append(cutoff_class_assignments[1], 1)
				cutoff_class_assignments[0] = np.append(cutoff_class_assignments[0], 0)

		tp_somatic, tn_somatic, fp_somatic, fn_somatic = 0, 0, 0, 0
		tp_noise, tn_noise, fp_noise, fn_noise = 0, 0, 0, 0
		for i in range(len(y_pred_proba)):
			# SOMATIC CONFUSION MATRIX VALUES
			if cutoff_class_assignments[1][i] == 1 and cutoff_class_assignments['truth'][i] == 1:
				# somatic = true and class = somatic
				tp_somatic += 1
			if cutoff_class_assignments[1][i] == 0 and cutoff_class_assignments['truth'][i] == 0:
				# somatic = false and class = noise
				tn_somatic += 1
			if cutoff_class_assignments[1][i] == 1 and cutoff_class_assignments['truth'][i] == 0:
				# somatic = true but class = noise
				fp_somatic += 1
			if cutoff_class_assignments[1][i] == 0 and cutoff_class_assignments['truth'][i] == 1:
				# somatic = fasle but class = somatic
				fn_somatic += 1
			# NOISE CONFUSION MATRIX VALUES
			if cutoff_class_assignments[0][i] == 1 and cutoff_class_assignments['truth'][i] == 0:
				# noise = true and class = noise
				tp_noise += 1
			if cutoff_class_assignments[0][i] == 0 and cutoff_class_assignments['truth'][i] == 1:
				# noise = false and class = somatic
				tn_noise += 1
			if cutoff_class_assignments[0][i] == 1 and cutoff_class_assignments['truth'][i] == 1:
				# noise = true but class = somatic
				fp_noise += 1
			if cutoff_class_assignments[0][i] == 0 and cutoff_class_assignments['truth'][i] == 0:
				# noise = false but class = noise
				fn_noise += 1

		# SOMATIC SCORES
		try:
			cutoff_results[cutoff]['somatic_precision'] = tp_somatic/(tp_somatic+fp_somatic)
		except ZeroDivisionError as e:
			cutoff_results[cutoff]['somatic_precision'] = 0.0
		try:
			cutoff_results[cutoff]['somatic_recall'] = tp_somatic/(tp_somatic+fn_somatic)
		except ZeroDivisionError as e:
			cutoff_results[cutoff]['somatic_recall'] = 0.0
		#f_measure = (2*precision*recall)/(precision+recall)
		try:
			cutoff_results[cutoff]['noise_precision'] = tp_noise/(tp_noise+fp_noise)
		except ZeroDivisionError as e:
			cutoff_results[cutoff]['noise_precision'] = 0.0
		try:
			cutoff_results[cutoff]['noise_recall'] = tp_noise/(tp_noise+fn_noise)
		except ZeroDivisionError as e:
			cutoff_results[cutoff]['noise_recall'] = 0.0
		#f_measure = (2*precision*recall)/(precision+recall)

	with open(f'{outdir}/cutoff_results_{model_type}.txt', 'w') as fp:
		fp.write('\t'.join(['cutoff','somatic_precision', 'somatic_recall', 'noise_precision', 'noise_recall']))
		fp.write('\n')
		for cutoff, results in cutoff_results.items():
			line = [
				str(cutoff),
				str(results['somatic_precision']),
				str(results['somatic_recall']),
				str(results['noise_precision']),
				str(results['noise_recall'])
			]
			fp.write('\t'.join(line))
			fp.write('\n')

def random_forest_training(hyperparameter, X_train, y_train):
	""" train and return a random forest classifier """
	if hyperparameter:
		# hyper-parameter testing
		param_dist = {'n_estimators': randint(100,1000), 'max_depth': randint(10,25)}
		rf = RandomForestClassifier(n_jobs=16)
		rand_search = RandomizedSearchCV(
			rf,
			param_distributions=param_dist,
			n_iter=5,
			cv=3,
			scoring='average_precision',
			verbose=3
		)
		rand_search.fit(X_train, y_train)
		print('Using Best Hyperparameters:')
		print(', '.join(f'{key}: {round(value,2)}' for key, value in rand_search.best_params_.items()))
		# store the best model
		random_forest = rand_search.best_estimator_
	else:
		# use default pre-deteremined hyperparameters
		random_forest = RandomForestClassifier(max_depth=20, n_estimators=400, n_jobs=16)
		random_forest.fit(X_train, y_train)

	return random_forest

def gradient_boost_training(hyperparameter, X_train, y_train):
	""" train and return a gradient boosting classifier """
	if hyperparameter:
		param_dist = {'n_estimators': randint(100,1000), 'learning_rate': uniform(0.01, 1.0), 'max_depth': randint(10,25)}
		gb = GradientBoostingClassifier()
		rand_search = RandomizedSearchCV(
			gb,
			param_distributions=param_dist,
			n_iter=5,
			cv=3,
			scoring='average_precision',
			verbose=3
		)
		rand_search.fit(X_train, y_train)
		print('Using Best Hyperparameters:')
		print(', '.join(f'{key}: {round(value,2)}' for key, value in rand_search.best_params_.items()))
		# store the best model
		gboostm = rand_search.best_estimator_
	else:
		gboostm = GradientBoostingClassifier()
		gboostm.fit(X_train, y_train)

	return gboostm

def evaluate_model(model, X_train, X_test, y_test, y_pred, outdir, model_type, germline_class):
	""" calculate performance stats of model """
	print(f'\nStatistics for model: {model_type}')
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

	print('\n- Precision-Recall Curve -')
	probs = model.predict_proba(X_test.drop(['SAMPLE','ID'], axis=1))
	precision, recall, _ = precision_recall_curve(y_test, probs[:, 1])
	print(f'F1: {round(f1_score(y_test, y_pred),2)} | AUC: {round(auc(recall, precision),2)}')

	# print out the feature importances
	print('\nSorted Feature Importances:')
	feature_importances = pd.Series(model.feature_importances_, index=X_train.drop(['SAMPLE','ID'], axis=1).columns).sort_values(ascending=False)
	print(feature_importances.to_string())

	print('\nConfusion Matrix:')
	labels = np.array(['FALSE', 'SOMATIC', 'GERMLINE']) if germline_class else np.array(['FALSE', 'TRUE'])
	cm=confusion_matrix(y_test, y_pred)
	disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
	disp.plot()
	plt.savefig(os.path.join(outdir, f'{model_type}_confusion_matrix.png'))
	plt.show()

	report = classification_report(y_test, y_pred)
	print(report)
	with open(os.path.join(outdir, f'{model_type}_model_stats.txt'), 'w') as f:
		f.write(report)

	# export y_test and pred
	test_matrix = X_test.copy(deep=True)
	test_matrix['TRUE_LABEL'] = y_test.copy()
	test_matrix['PREDICTED_LABEL'] = y_pred.copy()
	# only interested in the FP/FN
	test_matrix_missed = test_matrix[test_matrix['TRUE_LABEL'] != test_matrix['PREDICTED_LABEL']]
	test_matrix_missed.to_csv(os.path.join(outdir, f'test_set_incorrect_{model_type}.tsv'), sep="\t", index=False)
	# only interested in the TP/TN
	test_matrix_found = test_matrix[test_matrix['TRUE_LABEL'] == test_matrix['PREDICTED_LABEL']]
	test_matrix_found.to_csv(os.path.join(outdir, f'test_set_correct_{model_type}.tsv'), sep="\t", index=False)


def fit_classifier(X, y, outdir, split, downsample, hyperparameter, germline_class):
	""" given the features (X) and target (y), split into test/train and fit the model """
	# split into train/test
	X_train_with_ids, X_test_with_ids, y_train, y_test = train_test_split(X, y, test_size=split)

	# TODO: is downsampling needed if stratify arg is used for train/test?
	#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.6, stratify=y)

	# in train, randomly downsample the 'NOT_IN_COMPARISON'/0 by 'downsample' %
	train_matrix = pd.concat([X_train_with_ids, y_train], axis=1)
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
	X_train_with_ids = train_matrix.drop('LABEL', axis=1)
	y_train = train_matrix['LABEL']
	# remove IDs from both downsampled train and untouched test
	X_train = X_train_with_ids.iloc[: , 2:]
	X_test = X_test_with_ids.iloc[: , 2:]

	# train Random Forest
	random_forest_model = random_forest_training(hyperparameter, X_train, y_train)
	# use and evaluate
	y_pred_rf = random_forest_model.predict(X_test)
	evaluate_model(random_forest_model, X_train_with_ids, X_test_with_ids, y_test, y_pred_rf, outdir, "random_forest", germline_class)
	if False:
		# option to test randomforestclassifier at different cutoffs
		test_cutoffs(random_forest_model, X_test, y_test, outdir, "random_forest")
	# try with custom cutoff
	#custom_cutoff = 0.64
	custom_cutoff = None
	if custom_cutoff:
		y_pred_rf = (random_forest_model.predict_proba(X_test)[:,1] >= custom_cutoff).astype(bool)
		evaluate_model(random_forest_model, X_train, X_test, y_test, y_pred_rf, outdir, f'random_forest_cutoff-{custom_cutoff}', germline_class)

	"""
	gboost_model = gradient_boost_training(hyperparameter, X_train, y_train)
	y_pred_gb = gboost_model.predict(X_test)
	evaluate_model(gboost_model, X_train, X_test, y_test, y_pred_gb, outdir, "gradient_boosting", germline_class)
	"""

	return random_forest_model

def save_model(args, model, outdir):
	""" save model and its info to pickle file """
	model_path = os.path.join(outdir, 'random_forest_model.pkl')
	print(f'\nSaving classifier to {model_path}')
	try:
		model.set_params(n_jobs=1) # reset jobs before saving
	except AttributeError as e:
		# model doesn't have parallelisation option
		pass
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
				samples = input_vcf.samples
				if len(samples) != 1:
					raise ValueError(f'One sample per VCF permitted: {len(samples)} samples found in {f} ({",".join(samples)})')
				if not header:
					header = ['SAMPLE', 'ID']
					header.extend(get_info_fields(input_vcf))
				for variant in input_vcf:
					row = [samples[0], variant.ID]
					for field in header[2:]:
						row.append(variant.INFO.get(field))
					data.append(row)
	df = pd.DataFrame(data, columns = header)
	if args.save_matrix:
		print(f'Saving data matrix to pickle file {args.save_matrix}')
		df.to_pickle(args.save_matrix)
	return df

if __name__ == "__main__":
	print("Training functions for SAVANA")
