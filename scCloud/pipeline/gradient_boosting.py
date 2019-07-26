import numpy as np
import pandas as pd

from scCloud.tools import read_input, find_markers


def run_find_markers(input_h5ad_file, output_file, label_attr = 'louvain_labels', n_jobs = 1, min_gain = 1.0, random_state = 0, remove_ribo = False):
	data = read_input(input_h5ad_file, mode = 'a')
	markers = find_markers(data, label_attr, n_jobs = n_jobs, min_gain = min_gain, random_state = random_state, remove_ribo = remove_ribo)

	nclust = len(markers)
	keywords = [('strong', 'strong_gain'), ('weak', 'weak_gain'), ('down', 'down_gain')]

	writer = pd.ExcelWriter(output_file, engine='xlsxwriter')

	for i in range(nclust):
		sizes = []
		for keyword in keywords:
			sizes.append(len(markers[i][keyword[0]]))

		arr = np.zeros((max(sizes), 8), dtype = object)
		arr[:] = ''

		for j in range(3):
			arr[0:sizes[j], j * 3] = markers[i][keywords[j][0]]
			arr[0:sizes[j], j * 3 + 1] = markers[i][keywords[j][1]]

		df = pd.DataFrame(data = arr, columns = ['strongly up-regulated', 'gain', '', 'weakly up-regulated', 'gain', '', 'down-regulated', 'gain'])
		df.to_excel(writer, sheet_name = "{}".format(i + 1), index = False)

	writer.save()
