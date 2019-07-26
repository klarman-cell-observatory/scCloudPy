import unittest
import subprocess
import os
import scCloud as sc
import pandas as pd
import numpy as np
import scipy.sparse


class TestClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.remove('test_cluster.h5ad')

    def test_cluster(self):
        subprocess.call(
            ['scCloud', 'cluster', os.path.join('data', '3k_pbmc'),
             'test_cluster', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne'])
        test_data = sc.tools.read_input('test_cluster.h5ad')
        data = sc.tools.read_input(os.path.join('output', 'test_cluster.h5ad'))
        self.assertEqual((test_data.X[()] != data.X[()]).sum(), 0)
        pd.testing.assert_frame_equal(test_data.obs, data.obs)
        pd.testing.assert_frame_equal(test_data.var, data.var)
        self.assertListEqual(list(test_data.uns.keys()), list(data.uns.keys()))
        self.assertListEqual(list(test_data.obsm_keys()), list(data.obsm_keys()))
        for key in data.uns_keys():
            test_val = test_data.uns[key]
            val = data.uns[key]
            if scipy.sparse.issparse(val):
                val = val.toarray()
                test_val = test_val.toarray()
            np.testing.assert_array_equal(test_val, val)
        for key in data.obsm_keys():
            test_val = test_data.obsm[key]
            val = data.obsm[key]
            if scipy.sparse.issparse(val):
                val = val.toarray()
                test_val = test_val.toarray()
            np.testing.assert_array_equal(test_val, val)


# '--run-approximated-louvain',
# '--run-louvain',
# '--run-fle',
# '--run-net-fle'
# '--plot-hvg',


if __name__ == '__main__':
    unittest.main()
