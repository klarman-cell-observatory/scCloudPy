import unittest
import os
import scCloud.commands
from .test_util import assert_files_equal, assert_adata_files_equal


class TestClusterPipeline(unittest.TestCase):
    # def tearDown(self):
    #     os.remove('test_cluster.h5ad')
    #     os.remove('test_cluster.hvg.pdf')

    def test_cluster(self):
        cmd = scCloud.commands.cluster(
            ['cluster', os.path.join('tests', 'data', '3k_pbmc'),
             'test_cluster', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne', '--plot-hvg'])

        cmd.execute()
        # if is_running_in_docker:
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_cluster.h5ad'), 'test_cluster.h5ad')
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_cluster.hvg.pdf'), 'test_cluster.hvg.pdf')
        # else:
        # TODO, test umap reproducibility
        assert_adata_files_equal(self, os.path.join('tests', 'output', 'test_cluster.h5ad'), 'test_cluster.h5ad',
            obsm_blacklist=set(['X_net_umap_pred', 'X_net_umap', 'X_umap']))
        # '--run-louvain'
        # '--run-fle',
        # '--run-net-fle'


if __name__ == '__main__':
    unittest.main()
