#!/usr/bin/env python

import unittest
import os
import IMP.spatiotemporal as spatiotemporal

# General paths
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
data_path = os.path.join(TOPDIR, 'modeling', 'Trajectories', 'Trajectories_Modeling', 'data')
old_pdf_path = os.path.join(TOPDIR, 'modeling', 'Trajectories', 'Trajectories_Modeling', 'output')

# Path where the test output of spatiotemporal.create_DAG should be saved
output = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'output'))

# Paths for expected and generated files
expected_pdf_path = os.path.join(old_pdf_path, "pdf.txt")
generated_pdf_path = os.path.join(output, "pdf.txt")

# parameters
state_dict = {'0min': 3, '1min': 3, '2min': 1}
npaths = 3
expected_subcomplexes = ['A', 'B', 'C']
exp_comp = {
    'A': f'{data_path}/exp_compA.csv',
    'B': f'{data_path}/exp_compB.csv',
    'C': f'{data_path}/exp_compC.csv'
}

class TestSpatiotemporalDAG(unittest.TestCase):

    def test_create_dag_and_check_pdf(self):
        """Test if spatiotemporal.create_DAG creates pdf.txt and if the content matches the expected file."""

        # Call the function to generate output
        nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(
            state_dict,
            out_pdf=True,
            npaths=npaths,
            input_dir=data_path,
            scorestr='_scores.log',
            output_dir=output,
            spatio_temporal_rule=True,
            expected_subcomplexes=expected_subcomplexes,
            score_comp=True,
            exp_comp_map=exp_comp,
            draw_dag=False # there is no need for heatmap
        )

        # Check if the pdf.txt file was created
        self.assertTrue(os.path.exists(generated_pdf_path), "pdf.txt should exist in the output directory")

        # Check if the content of the generated pdf.txt matches the expected content
        with open(generated_pdf_path, 'r') as generated_file:
            generated_content = generated_file.read()

        with open(expected_pdf_path, 'r') as expected_file:
            expected_content = expected_file.read()

        # Compare the file contents
        self.assertEqual(generated_content, expected_content,
                         "The content of pdf.txt does not match the expected content")

    # Do we need clear out also generated files? I think GitHub Actions delete them after running test

if __name__ == '__main__':
    unittest.main()








