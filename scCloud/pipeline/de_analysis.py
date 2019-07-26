import time

import numpy as np
import pandas as pd
import xlsxwriter
from natsort import natsorted

from scCloud.tools import read_input
from scCloud.tools.de_analysis import test2fields, format_short_output_cols, __markers


def save_markers_to_excel(output_file, df, alpha=0.05):
    clusts = natsorted([x[10:] for x in df.columns if x.startswith("WAD_score_")])
    tests = [x for x in ['t', 'fisher', 'mwu'] if "{0}_qval_{1}".format(x, clusts[0]) in df.columns]
    has_roc = "auc_{0}".format(clusts[0]) in df.columns

    cols = ["percentage", "percentage_other", "percentage_fold_change", "mean_log_expression", "log_fold_change",
            "WAD_score"]
    if has_roc:
        cols.extend(["auc", "predpower"])
    if has_roc:
        cols.extend(["tpr_at_fpr01", "tpr_at_fpr025", "tpr_at_fpr03", "tpr_at_fpr05"])
    cols_short_format = cols.copy()
    for test in tests:
        cols.extend(test2fields[test])

    workbook = xlsxwriter.Workbook(output_file, {'nan_inf_to_errors': True})
    workbook.formats[0].set_font_size(9)
    for clust_id in clusts:
        idx = df["{0}_qval_{1}".format(tests[0], clust_id)] <= alpha
        for test in tests[1:]:
            idx = idx & (df["{0}_qval_{1}".format(test, clust_id)] <= alpha)

        idx_up = idx & (df["WAD_score_{0}".format(clust_id)] > 0.0)
        idx_down = idx & (df["WAD_score_{0}".format(clust_id)] < 0.0)
        assert idx_up.sum() + idx_down.sum() == idx.sum()

        col_names = ["{0}_{1}".format(x, clust_id) for x in cols]
        df_up = pd.DataFrame(df.loc[idx_up.values, col_names])
        df_up.rename(columns=lambda x: '_'.join(x.split('_')[:-1]), inplace=True)
        df_up.sort_values(by="auc" if has_roc else "WAD_score", ascending=False, inplace=True)
        # format output as excel table
        df_up = format_short_output_cols(df_up, cols_short_format)
        worksheet = workbook.add_worksheet(name="{0} up".format(clust_id))
        df_up.reset_index(inplace=True)
        df_up.rename(index=str, columns={"index": "gene"}, inplace=True)
        if len(df_up.index) > 0:
            worksheet.add_table(0, 0, len(df_up.index), len(df_up.columns) - 1,
                {'data': np.array(df_up), 'style': 'Table Style Light 1',
                 'first_column': True, 'header_row': True,
                 'columns': [{'header': x} for x in df_up.columns.values]})
        else:
            worksheet.write_row(0, 0, df_up.columns.values)

        df_down = pd.DataFrame(df.loc[idx_down.values, col_names])
        df_down.rename(columns=lambda x: '_'.join(x.split('_')[:-1]), inplace=True)
        df_down.sort_values(by="auc" if has_roc else "WAD_score", ascending=True, inplace=True)
        # format output as excel table
        worksheet = workbook.add_worksheet(name="{0} dn".format(clust_id))
        df_down = format_short_output_cols(df_down, cols_short_format)
        df_down.reset_index(inplace=True)
        df_down.rename(index=str, columns={"index": "gene"}, inplace=True)
        if len(df_up.index) > 0:
            worksheet.add_table(0, 0, len(df_down.index), len(df_down.columns) - 1,
                {'data': np.array(df_down), 'style': 'Table Style Light 1',
                 'first_column': True, 'header_row': True,
                 'columns': [{'header': x} for x in df_down.columns.values]})
        else:
            worksheet.write_row(0, 0, df_down.columns.values)
    workbook.close()

    print("Excel spreadsheet is written.")


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
    __markers(data, X, labels=labels, n_jobs=n_jobs, run_fisher=run_fisher, run_mwu=run_mwu, run_roc=run_roc,
        temp_folder=temp_folder)
    data.write(output_file)

    print("Differential expression results are written back to h5ad file.")

    save_markers_to_excel(output_excel_file, data.var, alpha=alpha)
