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


def assert_dict_of_arrays_equal(test_case, dict1, dict2, blacklist):
    # TODO handle nested
    if blacklist is None:
        blacklist = set()
    list1_keys = set(dict1.keys())
    list2_keys = set(dict2.keys())
    test_case.assertSetEqual(list1_keys.difference(blacklist), list2_keys.difference(blacklist))
    for key in dict1.keys():
        if key not in blacklist:
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


def assert_adata_files_equal(test_case, path, test_path, obs_blacklist=None, var_blacklist=None, uns_blacklist=None,
                             obsm_blacklist=None):
    test_data = sc.tools.read_input(test_path, mode='a')
    data = sc.tools.read_input(path, mode='a')
    assert_adata_equal(test_case, data, test_data, obs_blacklist, var_blacklist, uns_blacklist, obsm_blacklist)


def assert_adata_equal(test_case, data1, data2, obs_blacklist=None, var_blacklist=None, uns_blacklist=None,
                       obsm_blacklist=None):
    if scipy.sparse.issparse(data1.X):
        data1.X = data1.X.toarray()
    if scipy.sparse.issparse(data2.X):
        data2.X = data2.X.toarray()
    np.testing.assert_array_almost_equal(data1.X, data2.X)
    if obs_blacklist is not None:
        obs1 = data1.obs.drop(obs_blacklist, axis=1)
        obs2 = data2.obs.drop(obs_blacklist, axis=1)
    else:
        obs1 = data1.obs
        obs2 = data2.obs

    pd.testing.assert_frame_equal(obs1, obs2)
    if var_blacklist is not None:
        var1 = data1.var.drop(var_blacklist, axis=1)
        var2 = data2.var.drop(var_blacklist, axis=1)
    else:
        var1 = data1.var
        var2 = data2.var
    pd.testing.assert_frame_equal(var1, var2)
    assert_dict_of_arrays_equal(test_case, data1.uns, data2.uns, uns_blacklist)
    assert_dict_of_arrays_equal(test_case, data1.obsm, data2.obsm, obsm_blacklist)
