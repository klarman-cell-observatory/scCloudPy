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


def assert_adata_equal(test_case, path, test_path):
    test_data = sc.tools.read_input(test_path)
    data = sc.tools.read_input(path)
    X = data.X[()]
    test_X = test_data.X[()]
    if scipy.sparse.issparse(X):
        X = X.toarray()
        test_X = test_X.toarray()
    np.testing.assert_array_almost_equal(X, test_X)
    pd.testing.assert_frame_equal(test_data.obs, data.obs)
    pd.testing.assert_frame_equal(test_data.var, data.var)
    test_case.assertListEqual(list(test_data.uns.keys()), list(data.uns.keys()))
    test_case.assertListEqual(list(test_data.obsm_keys()), list(data.obsm_keys()))
    for key in data.uns_keys():
        test_val = test_data.uns[key]
        val = data.uns[key]
        if scipy.sparse.issparse(val):
            val = val.toarray()
            test_val = test_val.toarray()
        np.testing.assert_array_almost_equal(test_val, val)
    for key in data.obsm_keys():
        test_val = test_data.obsm[key]
        val = data.obsm[key]
        if scipy.sparse.issparse(val):
            val = val.toarray()
            test_val = test_val.toarray()
        np.testing.assert_array_almost_equal(test_val, val)
