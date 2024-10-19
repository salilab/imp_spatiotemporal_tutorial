#!/usr/bin/env python

import unittest
import os
import IMP
from IMP.spatiotemporal import prepare_protein_library


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
csv_path = os.path.join(TOPDIR, 'modeling', 'Input_Information', 'gen_FCS')
topology_path = os.path.join(TOPDIR, 'modeling', 'Snapshots', 'Snapshots_Modeling')

times = ["0min", "1min", "2min"]
exp_comp = {
    'A': f'{csv_path}/exp_compA.csv',
    'B': f'{csv_path}/exp_compB.csv',
    'C': f'{csv_path}/exp_compC.csv'
}
print(exp_comp)
expected_subcomplexes = ['A', 'B', 'C']
nmodels = 3
template_topology = f'{topology_path}/spatiotemporal_topology.txt'
print(template_topology)
template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}

class TestPrepareProteinLibraryContentRealData(unittest.TestCase):

    def test_3_0min_topol_content(self):
        """
        This test is focused only on 3_0min snapshot and tests two things:
        -if 3_0min_topol.txt topology is generated
        -if the content of 3_0min_topol.txt matches the expected content
        """
        # Calling function here directly
        prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                        template_topology=template_topology,
                                                        template_dict=template_dict)

        # Testing if 3_0min_topol.txt is created in testing directory
        generated_file_path = "3_0min_topol.txt"
        self.assertTrue(os.path.exists(generated_file_path), "3_0min_topol.txt should exist")

        # Comparing if content in 3_0min_topol.txt matches the expected content
        expected_content = ("|molecule_name | color | fasta_fn | fasta_id | pdb_fn | chain | residue_range | pdb_offset | bead_size | em_residues_per_gaussian | rigid_body | super_rigid_body | chain_of_super_rigid_bodies |\n\n"
        "|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|-1,18|2|1|10|1|1||\n"
        "|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|19,END|2|1|10|2|1||\n"
        "|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|1,83|0|1|10|3|2||\n"
        "|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|84,END|0|1|10|4|2||")

        # Read the content of the generated file
        with open(generated_file_path, 'r') as f:
            generated_content = f.read()

        # Compare both contents
        self.assertEqual(generated_content, expected_content,
                         "The content of 3_0min_topol.txt does not match the expected content.")
        # if this assert.Equal doesnt work with expected_content variable, maybe we can try that test is comparing
        # content of previous 3_0min_topol.txt (instead of expected_content variable) and newly created 3_0min_topol.txt

if __name__ == '__main__':
    unittest.main()

