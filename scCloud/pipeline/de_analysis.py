import time

from scCloud.tools import read_input
from scCloud.tools.de_analysis import __diff_exp, diff_exp_to_excel


def run_de_analysis(input_file, output_excel_file, labels, n_jobs, alpha, run_fisher, run_mwu, run_roc, subset_string,
                    temp_folder):
    start = time.time()
    if subset_string is None:
        data = read_input(input_file, mode='r+')
        X = data.X[:]
        output_file = input_file
    else:
        attr, value = subset_string.split(':')
        data = read_input(input_file, mode='a')
        data = data[data.obs[attr] == value].copy()
        X = data.X
        import os
        output_file = os.path.splitext(output_excel_file)[0] + '.h5ad'

    end = time.time()
    print("{0} is loaded. Time spent = {1:.2f}s.".format(input_file, end - start))
    __diff_exp(data, X, labels=labels, n_jobs=n_jobs, run_fisher=run_fisher, run_mwu=run_mwu, run_roc=run_roc,
        temp_folder=temp_folder)
    data.write(output_file)

    print("Differential expression results are written back to h5ad file.")

    diff_exp_to_excel(data.var, output_excel_file, alpha=alpha)
