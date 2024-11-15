Modeling of heterogeneity {#heterogeneity}
====================================

Here, we describe the first modeling problem in our composite workflow, how to build models of heterogeneity modeling using IMP. In this tutorial, heterogeneity modeling only includes protein copy number; however, in general, other types of information, such as the coarse location in the final state, could also be included in heterogeneity models.

# Heterogeneity modeling step 1: gathering of information {#heterogeneity1}

We begin heterogeneity modeling with the first step of integrative modeling, gathering information. Heterogeneity modeling will rely on copy number information about the complex. In this case, we utilize the X-ray crystal structure of the fully assembled Bmi1/Ring1b-UbcH5c complex from the protein data bank (PDB), and synthetically generated protein copy numbers during the assembly process.

\image html Input_heterogeneity.png width=600px

The PDB structure of the complex informs the final state of our model and constrains the maximum copy number for each protein, while the protein copy number data informs the stoichiometry of the complex during assembly.

# Heterogeneity modeling step 2: representation, scoring function, and search process {#heterogeneity2}

Next, we represent, score and search for heterogeneity models models. These operations are performed by the `heterogeneity_modeling.py` in the `Heterogeneity/Heterogeneity_Modeling` folder. A single heterogeneity model is a set of protein copy numbers, scored according to its fit to experimental copy number data at that time point. As ET and SAXS data, are only available at 0 minutes, 1 minute, and 2 minutes, we choose to create heterogeneity models at these three time points. We then use `prepare_protein_library`, [documented here](https://integrativemodeling.org/nightly/doc/ref/namespaceIMP_1_1spatiotemporal_1_1prepare__protein__library.html), to calculate the protein copy numbers for each snapshot model and to use the topology file of the full complex (`spatiotemporal_topology.txt`) to generate a topology file for each of these snapshot models. The choices made in this topology file are important for the representation, scoring function, and search process for snapshot models, and are [discussed later.] (@ref snapshot_representation) For heterogeneity modeling, we choose to model 3 protein copy numbers at each time point, and restrict the final time point to have the same protein copy numbers as the PDB structure. 

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

From the output of `prepare_protein_library`, we see that there are 3 heterogeneity models at each time point (it is possible to have more snapshot models than copy numbers if multiple copies of the protein exist in the complex). For each heterogeneity model, we see 2 files:
- *.config, a file with a list of proteins represented in the heterogeneity model
- *_topol.txt, a topology file for snapshot modeling corresponding to this heterogeneity model.

# Heterogeneity modeling step 3: assessment {#heterogeneity_assess}

\image html Heterogeneity_Assessment.png width=600px

