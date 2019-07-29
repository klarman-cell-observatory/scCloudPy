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

    for key in dict1:
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


def assert_adata_equal(test_case, path, test_path):
    test_data = sc.tools.read_input(test_path, mode='a')
    data = sc.tools.read_input(path, mode='a')

    if scipy.sparse.issparse(data.X):
        data.X = data.X.toarray()
    if scipy.sparse.issparse(test_data.X):
        test_data.X = test_data.X.toarray()
    np.testing.assert_array_almost_equal(data.X, test_data.X)
    pd.testing.assert_frame_equal(test_data.obs, data.obs)
    pd.testing.assert_frame_equal(test_data.var, data.var)
    assert_dict_of_arrays_equal(test_case, test_data.uns, data.uns)
    assert_dict_of_arrays_equal(test_case, test_data.obsm, data.obsm)
