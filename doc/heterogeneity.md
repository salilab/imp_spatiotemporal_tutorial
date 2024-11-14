Modeling of heterogeneity {#heterogeneity}
====================================

# Heterogeneity modeling step 1: gather information {#heterogeneity1}

# Heterogeneity modeling step 2: representation, scoring, and search process {#heterogeneity2}

We first must select which snapshots to model. Here, we choose only to model snapshots at 0 minutes, 1 minute, and 2 minutes because ET and SAXS data are only available at those time points. We know this complex has three protein chains (A, B, and C), and we choose to model these chains based on their protein copy number data. We then use `prepare_protein_library`, [documented here](https://integrativemodeling.org/nightly/doc/ref/namespaceIMP_1_1spatiotemporal_1_1prepare__protein__library.html), to calculate the protein copy numbers for each snapshot model and to use the topology file of the full complex (`spatiotemporal_topology.txt`) to generate a topology file for each of these snapshot models. Here, we choose to model 3 protein copy numbers at each time point, and restrict the final time point to have the same protein copy numbers as the PDB structure. 

\code{.py}
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
IMP.spatiotemporal.prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                template_topology=template_topology, template_dict=template_dict)
\endcode

From the output of `prepare_protein_library`, we see that there are 3 snapshot models at each time point (it is possible to have more snapshot models than copy numbers if multiple copies of the protein exist in the complex). We then wrote `generate_all_snapshots`, which creates a directory for each snapshot, copies the necessary files into that directory, and submits a job script to run sampling. The job script will likely need to be customized for the user's computer or cluster.


# Heterogeneity modeling step 3: assessment {#heterogeneity_assess}

