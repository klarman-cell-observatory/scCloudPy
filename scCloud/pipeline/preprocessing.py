import time

import numpy as np
import pandas as pd


def filter_data(data, output_filt = None, plot_filt = None, plot_filt_figsize = None, mito_prefix = 'MT-', min_genes = 500, max_genes = 6000, min_umis = 100, max_umis = 600000, percent_mito = 0.1, percent_cells = 0.0005, min_genes_on_raw = 100):
	start = time.time()

	data.obs['n_genes'] = data.X.getnnz(axis = 1)
	if data.obs['n_genes'].min() == 0: # if raw.h5
		data._inplace_subset_obs(data.obs['n_genes'].values >= min_genes_on_raw)
	data.obs['n_counts'] = data.X.sum(axis = 1).A1

	mito_prefixes = mito_prefix.split(',')

	def startswith(name):
		for prefix in mito_prefixes:
			if name.startswith(prefix):
				return True
		return False

	mito_genes = data.var_names.map(startswith).values.nonzero()[0]
	data.obs['percent_mito'] = data.X[:, mito_genes].sum(axis=1).A1 / np.maximum(data.obs['n_counts'].values, 1.0)

	if output_filt is not None:
		writer = pd.ExcelWriter(output_filt + '.filt.xlsx', engine='xlsxwriter')
		gb1 = data.obs.groupby('Channel')
		df_before = gb1.median()
		df_before = df_before.assign(total = gb1.size())
		df_before.rename(columns = {'n_genes' : 'median_n_genes_before', 'n_counts' : 'median_n_umis_before', 'percent_mito' : 'median_percent_mito_before'}, inplace = True)

	if plot_filt is not None:
		df_plot_before = data.obs[['Channel', 'n_genes', 'n_counts', 'percent_mito']].copy()
		df_plot_before.reset_index(drop = True, inplace = True)
		df_plot_before['status'] = 'original'

	# Filter cells
	obs_index = np.logical_and.reduce((data.obs['n_genes'] >= min_genes,
									   data.obs['n_genes'] < max_genes,
									   data.obs['n_counts'] >= min_umis,
									   data.obs['n_counts'] < max_umis,
									   data.obs['percent_mito'] < percent_mito))
	data._inplace_subset_obs(obs_index)

	if output_filt is not None:
		gb2 = data.obs.groupby('Channel')
		df_after = gb2.median()
		df_after = df_after.assign(kept = gb2.size())
		df_after.rename(columns = {'n_genes' : 'median_n_genes', 'n_counts' : 'median_n_umis', 'percent_mito' : 'median_percent_mito'}, inplace = True)
		df = pd.concat((df_before, df_after), axis = 1, sort = False)
		df.fillna(0, inplace = True)
		df['kept'] = df['kept'].astype(int)
		df['filt'] = df['total'] - df['kept']
		df = df[['kept', 'median_n_genes', 'median_n_umis', 'median_percent_mito', 'filt', 'total', 'median_n_genes_before', 'median_n_umis_before', 'median_percent_mito_before']]
		df.sort_values('kept', inplace = True)
		df.to_excel(writer, sheet_name = "Cell filtration stats")

	if plot_filt is not None:
		df_plot_after = data.obs[['Channel', 'n_genes', 'n_counts', 'percent_mito']].copy()
		df_plot_after.reset_index(drop = True, inplace = True)
		df_plot_after['status'] = 'filtered'
		df_plot = pd.concat((df_plot_before, df_plot_after), axis = 0)
		from scCloud.plotting import plot_qc_violin
		figsize = None
		if plot_filt_figsize is not None:
			width, height = plot_filt_figsize.split(',')
			figsize = (int(width), int(height))
		plot_qc_violin(df_plot, 'count', plot_filt + '.filt.UMI.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		plot_qc_violin(df_plot, 'gene', plot_filt + '.filt.gene.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		plot_qc_violin(df_plot, 'mito', plot_filt + '.filt.mito.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		print("Filtration plots are generated.")

	# Filter genes
	data.var['n_cells'] = data.X.getnnz(axis = 0)
	data.var['percent_cells'] = data.var['n_cells'] / data.shape[0]
	data.var['robust'] = data.var['percent_cells'] >= percent_cells

	data.var['highly_variable_genes'] = data.var['robust'] # default all robust genes are "highly" variable
	data.var['hvg_rank'] = -1 # default all ranks are -1

	if output_filt is not None:
		idx = data.var['robust'] == False
		df = pd.DataFrame({'n_cells': data.var.loc[idx, 'n_cells'], 'percent_cells': data.var.loc[idx, 'percent_cells']})
		df.index.name = 'gene'
		df.sort_values('n_cells', ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Gene filtration stats")
		writer.save()
		print("Filtration results are written.")

	var_index = (data.var['n_cells'] > 0).values
	data._inplace_subset_var(var_index)
	print("After filteration, {nc} cells and {ng} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(nc = data.shape[0], ng = data.shape[1], nrb = data.var['robust'].sum()))

	end = time.time()
	print("filter_data is finished. Time spent = {:.2f}s.".format(end - start))
