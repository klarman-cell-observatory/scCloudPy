import unittest
import os
import scCloud.commands
from .test_util import assert_adata_files_equal


class TestClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.path.exists('test_cluster.h5ad') and os.remove('test_cluster.h5ad')
        os.path.exists('test_cluster.hvg.pdf') and os.remove('test_cluster.hvg.pdf')

    def test_cluster(self):
        cmd = scCloud.commands.cluster(
            ['cluster', os.path.join('tests', 'scCloud-test-data', 'input', '3k_pbmc'),
             'test_cluster', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne', '--run-fle',
             '--run-net-fle', '--plot-hvf', '--run-louvain', '--run-approximated-louvain'])
        cmd.execute()

        # TODO diff pdfs
        assert_adata_files_equal(self, os.path.join('tests', 'scCloud-test-data', 'output', 'test_cluster.h5ad'),
            'test_cluster.h5ad')


if __name__ == '__main__':
    unittest.main()