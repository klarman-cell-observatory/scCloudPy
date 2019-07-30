import unittest
import os
import scCloud as sc
import scCloud.commands
from .test_util import assert_adata_files_equal


class TestClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.remove('test_cluster.h5ad')

    def test_cluster(self):
        cmd = scCloud.commands.cluster(
            ['cluster', os.path.join('tests', 'data', '3k_pbmc'),
             'test_cluster', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne'])
        cmd.execute()
        assert_adata_files_equal(self, os.path.join('tests', 'output', 'test_cluster.h5ad'), 'test_cluster.h5ad',
            obs_blacklist=set(['approx_leiden_labels']), var_blacklist=set(['hvg_rank']))
        # approx_leiden_labels for 3rd cell differs between mac and linux
        test_data = sc.tools.read_input('test_cluster.h5ad', mode='a')
        data = sc.tools.read_input(os.path.join('tests', 'output', 'test_cluster.h5ad'), mode='a')
        n = (data.obs['approx_leiden_labels'] != test_data.obs['approx_leiden_labels']).sum()
        self.assertLessEqual(n, 1, '{} approx_leiden_labels differ'.format(n))
        n = (data.var['hvg_rank'] != test_data.var['hvg_rank']).sum()
        self.assertLessEqual(n, 1, '{} hvg_rank differ'.format(n))
        # '--run-approximated-louvain',
        # '--run-louvain',
        # '--run-fle',
        # '--run-net-fle'
        # '--plot-hvg',


if __name__ == '__main__':
    unittest.main()
