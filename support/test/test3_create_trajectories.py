#!/usr/bin/env python

import unittest
import os
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import IMP.test

# General paths
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
data_path = os.path.join(TOPDIR, 'modeling', 'Trajectories', 'Trajectories_Modeling', 'data')
old_pdf_path = os.path.join(TOPDIR, 'modeling', 'Trajectories', 'Trajectories_Modeling', 'output')



# Paths for expected and generated files
expected_pdf_path = os.path.join(old_pdf_path, "labeled_pdf.txt")

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

        # create output directory
        with IMP.test.temporary_directory() as tmpdir:
            output = os.path.join(tmpdir, 'output/')

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
            # Check if the labeled_pdf.txt file was created
            generated_pdf_path = os.path.join(output, "labeled_pdf.txt")
            self.assertTrue(os.path.exists(generated_pdf_path), "labeled_pdf.txt should exist in the output directory")

            # Use temporal_precision to check between old and new model
            analysis.temporal_precision(expected_pdf_path,generated_pdf_path,output_fn=output+'trj_test.txt')
            f=open(output+'trj_test.txt','r')
            # First line is description
            f.readline()
            # 2nd line is the temporal precision
            trj_comp=float(f.readline())
            f.close()
            self.assertAlmostEqual(trj_comp, 1.0, delta=1e-5)
            # Test the precision of the trajecotry
            analysis.precision(generated_pdf_path,output_fn=output+'precision.txt')
            f = open(output + 'precision.txt', 'r')
            # First line is description
            f.readline()
            # 2nd line is the precision
            precision = float(f.readline())
            f.close()
            self.assertAlmostEqual(precision, 1.0, delta=1e-5)

if __name__ == '__main__':
    unittest.main()








