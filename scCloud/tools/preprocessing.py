import time
import re
import numpy as np
import pandas as pd

from collections import Counter

from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def update_var_names(data):
    start = time.time()
    if 'genome' in data.uns:
        prefix = re.compile('^(' + '|'.join(data.uns['genome'].split(',')) + ')_+')
        if prefix.match(data.var_names[0]):
            data.var['gene_ids'] = [prefix.sub('', x) for x in data.var['gene_ids']]
            data.var_names = pd.Index([prefix.sub('', x) for x in data.var_names])

    gsyms = data.var_names.values

    dup_ids = Counter()
    for i in range(gsyms.size):
        idn = dup_ids[gsyms[i]]
        dup_ids[gsyms[i]] += 1
        if idn > 0:
            gsyms[i] = gsyms[i] + ".{}".format(idn)

    data.var_names = pd.Index(gsyms)

    end = time.time()
    print("update_var_names is finished. Time spent = {:.2f}s.".format(end - start))


def qc_metrics(data, mito_prefix='MT-', percent_cells=0.0005):
    """
    Sets n_genes, n_counts, percent_mito on adata.obs and n_cells, percent_cells, and robust on data.var

    :param data:
        Annotated data matrix
    :param mito_prefix: str
        String that mitochrondrial genes start with
    :param percent_cells: float
        Cutoff for a feature to be `robust`
    """
    data.obs['n_genes'] = data.X.getnnz(axis=1)
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    mito_prefixes = mito_prefix.split(',')

    def startswith(name):
        for prefix in mito_prefixes:
            if name.startswith(prefix):
                return True
        return False

    mito_genes = data.var_names.map(startswith).values.nonzero()[0]
    data.obs['percent_mito'] = data.X[:, mito_genes].sum(axis=1).A1 / np.maximum(data.obs['n_counts'].values, 1.0)
    data.var['n_cells'] = data.X.getnnz(axis=0)
    data.var['percent_cells'] = data.var['n_cells'] / data.shape[0]
    data.var['robust'] = data.var['percent_cells'] >= percent_cells


def filter_cells_cite_seq(data, max_cells):
    assert issparse(data.X)
    data.obs['n_counts'] = data.X.sum(axis=1).A1
    obs_index = np.zeros(data.shape[0], dtype=bool)
    obs_index[np.argsort(data.obs['n_counts'].values)[::-1][:max_cells]] = True
    data._inplace_subset_obs(obs_index)
    data.var['robust'] = True
    print("After filteration, {nc} cells are kept, with the minimum nUMI = {numi}.".format(nc=max_cells,
        numi=data.obs['n_counts'].min()))


def log_norm(data, norm_count):
    """ Normalization and then take log """
    start = time.time()

    assert issparse(data.X)
    mat = data.X[:, data.var['robust'].values]
    scale = norm_count / mat.sum(axis=1).A1
    data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
    data.X = data.X.log1p()

    end = time.time()
    print("Normalization is finished. Time spent = {:.2f}s.".format(end - start))


def pca(data, standardize=True, max_value=10, nPC=50, random_state=0, features=None):
    start = time.time()
    orig_data = data
    if features is not None:
        data = data[:, data.var[features]].copy()

    if issparse(data.X):
        data.X = data.X.toarray()

    if standardize:
        scaler = StandardScaler(copy=False)
        scaler.fit_transform(data.X)

    if max_value is not None:
        data.X[data.X > max_value] = max_value
        data.X[data.X < -max_value] = -max_value

    pca = PCA(n_components=nPC, random_state=random_state)
    X_pca = pca.fit_transform(data.X)
    orig_data.obsm['X_pca'] = X_pca
    # orig_data.varm['PCs'] = pca.components_.T
    orig_data.uns['pca'] = {}
    orig_data.uns['pca']['components'] = pca.components_
    orig_data.uns['pca']['variance'] = pca.explained_variance_
    orig_data.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
    end = time.time()
    print("PCA is done. Time spent = {:.2f}s.".format(end - start))


# def run_rpca(data, scale = False, max_value = 10.0, nPC = 50, random_state = 0):
# 	""" smooth outliers, then no center/scale data """
# 	start = time.time()

# 	# Smooth out outliers
# 	means, variances = mean_variance_axis(data.X, axis = 0)
# 	stds = np.sqrt(variances * (data.X.shape[0] / (data.X.shape[0] - 1))) # make it unbiased
# 	assert (stds == 0.0).sum() == 0

# 	data_new = (data.X.data - means[data.X.indices]) / stds[data.X.indices]
# 	outliers = data_new > max_value
# 	data.X.data[outliers] = max_value * stds[data.X.indices[outliers]] + means[data.X.indices[outliers]]

# 	if scale:
# 		data.X.data /= stds[data.X.indices]

# 	U, S, VT = randomized_svd(data.X, n_components = nPC, random_state = random_state)
# 	data.obsm['X_rpca'] = U * S

# 	end = time.time()
# 	print("RPCA is done. Time spent = {:.2f}s.".format(end - start))
