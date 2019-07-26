import unittest
import subprocess
import os
import shutil
import pandas as pd


class TestDePipeline(unittest.TestCase):

    def tearDown(self):
        os.remove('out_marker.de.xlsx')

    def test_de_analysis(self):
        # de_analysis modifies h5ad file
        shutil.copy(os.path.join('output', 'test_cluster.h5ad'), 'test.h5ad')
        subprocess.run(
            ['scCloud', 'de_analysis', 'test.h5ad', 'out_marker.de.xlsx', '--fisher',
             '--mwu', '--roc', '--labels', 'leiden_labels'])
        test_output = pd.read_excel('out_marker.de.xlsx', sheet_name=None)
        output = pd.read_excel(os.path.join('output', 'out_marker.de.xlsx'), sheet_name=None)
        self.assertListEqual(list(test_output.keys()), list(output.keys()))
        for key in output:
            pd.testing.assert_frame_equal(test_output[key], output[key])


if __name__ == '__main__':
    unittest.main()
