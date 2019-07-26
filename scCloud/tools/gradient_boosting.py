import time
import numpy as np
import pandas as pd
from collections import defaultdict

from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
# from xgboost import XGBClassifier
from lightgbm import LGBMClassifier


def find_markers(data, label_attr, n_jobs = 1, min_gain = 1.0, random_state = 0, remove_ribo = False):
	if remove_ribo:
		data = data[:,np.vectorize(lambda x: not x.startswith('RPL') and not x.startswith('RPS'))(data.var_names)]

	X_train, X_test, y_train, y_test = train_test_split(data.X, data.obs[label_attr], test_size = 0.1, random_state = random_state, stratify = data.obs[label_attr])

	# start = time.time()
	# xgb = XGBClassifier(n_jobs = n_jobs, n_gpus = 0)
	# xgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], eval_metric = 'merror')
	# # print(xgb.evals_result())
	# end = time.time()
	# print("XGBoost used {:.2f}s to train.".format(end - start))

	start = time.time()
	lgb = LGBMClassifier(n_jobs = n_jobs, metric = 'multi_error', importance_type = 'gain')
	lgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], early_stopping_rounds = 1)
	end = time.time()
	print("LightGBM used {:.2f}s to train.".format(end - start))

	ntot = (lgb.feature_importances_ >= min_gain).sum()
	ords = np.argsort(lgb.feature_importances_)[::-1][:ntot]

	ncat = data.obs[label_attr].cat.categories.size
	log_exprs = ['mean_log_expression_{}'.format(i + 1) for i in range(ncat)]
	titles = [('down', 'down_gain'), ('weak', 'weak_gain'), ('strong', 'strong_gain')]
	markers = defaultdict(lambda: defaultdict(list))

	kmeans = KMeans(n_clusters = 3, random_state = random_state)
	for gene_id in ords:
		gene_symbol = data.var_names[gene_id]
		mydat = data.var.loc[gene_symbol, log_exprs].values.reshape(-1, 1)
		kmeans.fit(mydat)
		kmeans_label_mode = pd.Series(kmeans.labels_).mode()[0]
		for i, kmeans_label in enumerate(np.argsort(kmeans.cluster_centers_[:,0])):
			if kmeans_label != kmeans_label_mode:
				for clust_label in (kmeans.labels_ == kmeans_label).nonzero()[0]:
					markers[clust_label][titles[i][0]].append(gene_symbol)
					markers[clust_label][titles[i][1]].append('{:.2f}'.format(lgb.feature_importances_[gene_id]))

	end = time.time()
	print("find_markers took {:.2f}s to finish.".format(end - start))

	return markers



