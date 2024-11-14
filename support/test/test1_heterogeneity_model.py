#!/usr/bin/env python

import unittest
import os
import IMP
from IMP.spatiotemporal import prepare_protein_library
import IMP.test


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
csv_path = os.path.join(TOPDIR, 'modeling', 'Input_Information', 'gen_FCS')
topology_path = os.path.join(TOPDIR, 'modeling', 'Snapshots', 'Snapshots_Modeling')

times = ["0min", "1min", "2min"]
exp_comp = {
    'A': f'{csv_path}/exp_compA.csv',
    'B': f'{csv_path}/exp_compB.csv',
    'C': f'{csv_path}/exp_compC.csv'
}
expected_subcomplexes = ['A', 'B', 'C']
nmodels = 3
template_topology = f'{topology_path}/spatiotemporal_topology.txt'
template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}

class TestPrepareProteinLibraryContentRealData(unittest.TestCase):

    def test_3_0min_topol_content(self):
        """
        This test tests prepare_protein_library by ensuring that 3_0min is correctly configured.
        """
        # create output directory
        with IMP.test.temporary_directory() as tmpdir:
            os.chdir(tmpdir)
            # Calling function here directly
            prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                            template_topology=template_topology,
                                                            template_dict=template_dict)

            # Testing if 3_0min_topol.txt is created in testing directory
            generated_file_path = os.path.join(tmpdir, "3_0min_topol.txt")
            self.assertTrue(os.path.exists(generated_file_path), "3_0min_topol.txt should exist")

            # Read the content of the generated file
            with open(generated_file_path, 'r') as f:
                generated_content = f.read()

            # Read in old topology file
            old_topol_file = os.path.join(topology_path, "3_0min_topol.txt")
            with open(old_topol_file, 'r') as f:
                old_content = f.read()

            # Compare both contents
            self.assertEqual(generated_content, old_content,
                             "The content of 3_0min_topol.txt does not match the previous topology file.")

if __name__ == '__main__':
    unittest.main()

