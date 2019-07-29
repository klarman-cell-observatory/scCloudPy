import pandas as pd
import scipy
import numpy as np
import scCloud as sc


def assert_excel_equal(test_case, path, test_path):
    test_output = pd.read_excel(test_path, sheet_name=None)
    output = pd.read_excel(path, sheet_name=None)
    test_case.assertListEqual(list(test_output.keys()), list(output.keys()))
    for key in output:
        pd.testing.assert_frame_equal(test_output[key], output[key])


def assert_dict_of_arrays_equal(test_case, dict1, dict2):
    test_case.assertListEqual(list(dict1.keys()), list(dict2.keys()))

    for key in dict1.keys():
        val1 = dict1[key]
        val2 = dict2[key]
        if scipy.sparse.issparse(val1):
            val1 = val1.toarray()
        if scipy.sparse.issparse(val2):
            val2 = val2.toarray()
        if isinstance(val1, np.ndarray):
            np.testing.assert_array_almost_equal(val1, val2, err_msg=str(key) + ' not equal')
        else:
            np.testing.assert_array_equal(val1, val2, err_msg=str(key) + ' not equal')


def assert_adata_files_equal(test_case, path, test_path):
    test_data = sc.tools.read_input(test_path, mode='a')
    data = sc.tools.read_input(path, mode='a')
    assert_adata_equal(test_case, data, test_data)


def assert_adata_equal(test_case, data1, data2):
    if scipy.sparse.issparse(data1.X):
        data1.X = data1.X.toarray()
    if scipy.sparse.issparse(data2.X):
        data2.X = data2.X.toarray()
    np.testing.assert_array_almost_equal(data1.X, data2.X)
    pd.testing.assert_frame_equal(data1.obs, data2.obs)
    pd.testing.assert_frame_equal(data1.var, data2.var)
    assert_dict_of_arrays_equal(test_case, data1.uns, data2.uns)
    assert_dict_of_arrays_equal(test_case, data1.obsm, data2.obsm)
