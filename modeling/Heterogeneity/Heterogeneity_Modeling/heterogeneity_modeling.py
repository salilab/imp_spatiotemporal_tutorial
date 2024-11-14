"""
This code models the heterogeneity of the system. This is mostly performed by an internal function of IMP:
-prepare_protein_library (IMP.spatiotemporal.prepare_protein_library.prepare_protein_library): based on
time - dependent stoichiometry data this function calculates configurations of possible snapshot models for
each time together with corresponding pmi topology file.
For example: for our system for first two time points we calculated three most possible configurations
(check state_dict - important dictionary that describe our snapshot model).
For more information regarding prepare_protein_library check:
https://integrativemodeling.org/nightly/doc/ref/namespaceIMP_1_1spatiotemporal_1_1prepare__protein__library.html
"""

from IMP.spatiotemporal import prepare_protein_library


if __name__ == "__main__":
    # 1a - parameters for prepare_protein_library:
    times = ["0min", "1min", "2min"]
    exp_comp = {'A': '../../Input_Information/gen_FCS/exp_compA.csv',
                'B': '../../Input_Information/gen_FCS/exp_compB.csv',
                'C': '../../Input_Information/gen_FCS/exp_compC.csv'}
    expected_subcomplexes = ['A', 'B', 'C']
    template_topology = 'spatiotemporal_topology.txt'
    template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}
    nmodels = 3

    # 1b - calling prepare_protein_library
    prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                    template_topology=template_topology, template_dict=template_dict)
