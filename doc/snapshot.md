Modeling of snapshots {#snapshots}
====================================

Here, we describe the second modeling problem in our composite workflow, how to build models of static snapshot models using IMP. We note that this process is similar to previous tutorials of [actin](https://integrativemodeling.org/tutorials/actin/) and [RNA PolII](https://integrativemodeling.org/tutorials/rnapolii_stalk/).

# Snapshot modeling step 1: gathering of information {#snapshots1}

We begin snapshot modeling with the first step of integrative modeling, gathering information. Snapshot modeling utilizes structural information about the complex. In this case, we utilize heterogeneity models, the X-ray crystal structure of the fully assembled Bmi1/Ring1b-UbcH5c complex from the protein data bank (PDB), synthetically generated electron tomography (ET) density maps during the assembly process, and physical theories.

\image html Input_snapshot.png width=600px

The heterogeneity models inform protein copy numbers for the snapshot models. The PDB structure of the complex informs the structure of the individual proteins. The time-dependent ET data informs the size and shape of the assembling complex. Physical theories inform connectivity and excluded volume.

# Snapshot modeling step 2: representation, scoring function, and search process {#snapshots2}

Next, we represent, score and search for snapshot models. To do so, navigate to the `Snapshots/Snapshots_Modeling/` folder. Here, you will find two python scripts. The first, `static_snapshot.py`, uses IMP to represent, score, and search for models of a single static snapshot. The second, `start_sim.py`, automates the creation of a snapshot model for each heterogeneity model.

## Modeling one snapshot

Here, we will describe the process of modeling a single snapshot model, as performed by running `static_snapshot.py`.

### Representing the model {#snapshot_representation}

We begin by representing the data and the model. In general, the *representation* of a system is defined by all the variables that need to be determined.

For our model of a protein complex, we use a combination of two representations. The first is a series of *spherical beads*, which can correspond to portions of the biomolecules of interest, such as atoms or groups of atoms. The second is a series of *3D Gaussians*, which help calculate the overlap between our model and the density from ET data.

Beads and Gaussians in our model belong to either a *rigid body* or *flexible string*. The positions of all beads and Gaussians in a single rigid body are constrained during sampling and do not move relative to each other. Meanwhile, flexible beads can move freely during sampling, but are restrained by sequence connectivity.

To begin, we built a topology file with the representation for the model of the complete system, `spatiotemporal_topology.txt`, located in `Heterogeneity/Heterogeneity_Modeling/`. This complete topology was used as a template to build topologies of each snapshot. Based on our observation of the structure of the complex, we chose to represent each protein with at least 2 separate rigid bodies, and left the first 28 residues of protein C as flexible beads. Rigid bodies were described with 1 bead for every residue, and 10 residues per Gaussian. Flexible beads were described with 1 bead for every residue and 1 residue per Gaussian. A more complete description of the options available in topology files is available in the the [TopologyReader](@ref IMP::pmi::topology::TopologyReader) documentation.

\code{.txt}
|molecule_name | color | fasta_fn | fasta_id | pdb_fn | chain | residue_range | pdb_offset | bead_size | em_residues_per_gaussian | rigid_body | super_rigid_body | chain_of_super_rigid_bodies | 

|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|-1,18|2|1|10|1|1||
|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|19,147|2|1|10|2|1||
|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|3,83|-2|1|10|3|2||
|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|84,101|-2|1|10|4|2||
|E3-ubi-RING2|green|3rpg.fasta.txt|E3-ubi-RING2|BEADS|C|16,44|-15|1|1|5|3||
|E3-ubi-RING2|green|3rpg.fasta.txt|E3-ubi-RING2|3rpg.pdb|C|45,116|-15|1|10|6|3||
\endcode

Next, we must prepare `static_snapshot.py` to read in this topology file. We begin by defining the input variables, `state` and `time`, which define which topology to use, as well as the paths to other pieces of input information.

\code{.py}
# Running parameters to access correct path of ET_data for EM restraint
# and topology file for certain {state}_{time}_topol.txt
state = sys.argv[1]
time = sys.argv[2]

# Topology file
topology_file = f"../{state}_{time}_topol.txt"
# Paths to input data for topology file
pdb_dir = "../../../../Input_Information/PDB"
fasta_dir = "../../../../Input_Information/FASTA"
# Path where forward gmms are created with BuildSystem (based ont topology file)
# If gmms exist, they will be used from this folder
forward_gmm_dir = "../forward_densities/"
# Path to experimental gmms
exp_gmm_dir= '../../../../Input_Information/ET_data/add_noise'
\endcode

Next, we build the system, using the topology tile, described above.
\code{.py}
# Create a system from a topology file. Resolution is set on 1.
bs = IMP.pmi.macros.BuildSystem(mdl, resolutions= 1, name= f'Static_snapshots_{state}_{time}')
bs.add_state(t)
\endcode

Then, we prepare for later sampling steps by setting which Monte Carlo moves will be performed. Rotation (`rot`) and translation (`trans`) parameters are separately set for super rigid bodies (`srb`), rigid bodies (`rb`), and beads (`bead`).
\code{.py}
#  Macro execution: It gives hierarchy and degrees of freedom (dof).
# In dof we define how much can each (super) rigid body translate and rotate between two adjacent Monte Carlo steps
root_hier, dof = bs.execute_macro(max_rb_trans=1.0,
                                  max_rb_rot=0.5, max_bead_trans=2.0,
                                  max_srb_trans=1.0, max_srb_rot=0.5)
\endcode

### Scoring the model {#snapshot_scoring}

After building the model representation, we choose a scoring function to score the model based on input information. This scoring function is represented as a series of restraints that serve as priors.

#### Connectivity

We begin with a connectivity restraint, which restrains beads adjacent in sequence to be close in 3D space.

\code{.py}
# Adding Restraints
# Empty list where the data from restraints should be collected
output_objects=[]

# Two common restraints: ConnectivityRestraint and ExcludedVolumeSphere
# ConnectivityRestraint is added for each "molecule" separately
for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)
\endcode

#### Excluded volume

Next is an excluded volume restraint, which restrains beads to minimize their spatial overlap.

\code{.py}
# Add excluded volume
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                     included_objects=[root_hier],
                                     resolution=1000)
output_objects.append(evr)
evr.add_to_model()
\endcode

#### Electron tomography

Finally, we restrain our models based on their fit to ET density maps. Both the experimental map and the forward protein density are represented as Gaussian mixture models (GMMs) to speed up scoring. The score is based on the log of the correlation coefficient between the experimental density and the forward protein density.

\code{.py}
# Applying time-dependent EM restraint. Point to correct gmm / mrc file at each time point
# Path to corresponding .gmm file (and .mrc file)
em_map = exp_gmm_dir + f"/{time}_noisy.gmm"

# Create artificial densities from hierarchy
densities = IMP.atom.Selection(root_hier,
                 representation_type=IMP.atom.DENSITIES).get_selected_particles()

# Create EM restraint based on these densities
emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,
        target_fn=em_map,
        slope=0.000001,
        scale_target_to_mass=True,
        weight=1000)
output_objects.append(emr)
emr.add_to_model()
\endcode

### Searching for good scoring models {#snapshot_searching}

After building a scoring function that scores alternative models based on their fit to the input information, we aim to search for good scoring models. For complicated systems, stochastic sampling techniques such as Monte Carlo (MC) sampling are often the most efficient way to compute good scoring models. Here, we generate a random initial configuration and then perform temperature replica exchange MC sampling with 16 temperatures from different initial configurations. By performing multiple runs of replica exchange MC from different initial configurations, we can later ensure that our sampling is sufficiently converged.

\code{.py}
# Generate random configuration
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50)

# Perform replica exchange sampling
rex=IMP.pmi.macros.ReplicaExchange(mdl,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory='output', # name 'output' is the best for imp sampcon select_good
        output_objects=output_objects,
        monte_carlo_steps=200, # Number of MC steps between writing frames.
        number_of_best_scoring_models=0,
        number_of_frames=500) # number of frames to be saved
# In our case, for each snapshot we generated 25000 frames altogether (50*500)
rex.execute_macro()
\endcode

After performing sampling, a variety of outputs will be created. These outputs include `.rmf` files, which contain multi-resolution models output by IMP, and `.out` files which contains a variety of information about the run such as the value of the restraints and the MC acceptance rate.

## Generalizing modeling to all snapshots {#snapshot_combine}

Next, we will describe the process of modeling a multiple static snapshots, as performed by running `start_sim.py`.

From heterogeneity modeling, we see that there are 3 heterogeneity models at each time point (it is possible to have more snapshot models than copy numbers if multiple copies of the protein exist in the complex), each of which has a corresponding topology file in `Heterogeneity/Heterogeneity_Modeling/`. We wrote a function, `generate_all_snapshots`, which creates a directory for each snapshot, copies the python script and topology file into that directory, and submits a job script to run sampling. The job script will likely need to be customized for the user's computer or cluster.

\code{.py}
# 1a - parameters for generate_all_snapshots
# state_dict - universal parameter
state_dict = {'0min': 3, '1min': 3, '2min': 1}

main_dir = os.getcwd()
topol_dir = os.path.join(os.getcwd(), '../../Heterogeneity/Heterogeneity_Modeling')
items_to_copy = ['static_snapshot.py']  # additionally we need to copy only specific topology file
# jobs script will likely depend on the user's cluster / configuration
job_template = ("#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -N Tutorial\n#$ -pe smp 16\n"
                "#$ -l h_rt=48:00:00\n\nmodule load Sali\nmodule load imp\nmodule load mpi/openmpi-x86_64\n\n"
                "mpirun -np $NSLOTS python3 static_snapshot.py {state} {time}")
number_of_runs = 50

# 1b - calling generate_all_snapshots
generate_all_snapshots(state_dict, main_dir, topol_dir, items_to_copy, job_template, number_of_runs)

\endcode

We note that sometimes errors such as the one below can arise during sampling. These errors are caused by issues generating forward GMM files, which is done stochastically. If such issues arrise, remove all files in the `forward_densities` folder for that snapshot and resubmit the corresponding jobs.

\code{.py}
  File "/imp/main/20240607-af6f9d6a95/lib/release8/IMP/isd/gmm_tools.py", line 35, in decorate_gmm_from_text
    weight = float(fields[2])
IndexError: list index out of range
\endcode

# Snapshot modeling step 3: assessment {#snapshot_assess}

Now, we have a variety of alternative snapshot models. In general, we would like to assess these models in at least 4 ways: estimate the sampling precision, compare the model to data used to construct it, validate the model against data not used to construct it, and quantify the precision of the model. Here, we will focus specifically on estimating the sampling precision of the model, while quantitative comparisons between the model and experimental data will be reserved for the final step, when we assess [trajectories](https://integrativemodeling.org/tutorials/spatiotemporal/trajectory_assess.html). To assess these snapshot models, we navigate to the `Snapshots/Snapshots_Assessment` folder and run `snapshot_assessment.py`. This script performs the following analysis.

## Filtering good scoring models {#snapshot_filter}

Initially, we want to filter the various alternative models to select those that meet certain parameter thresholds. In this case, we filter the structural models in each snapshot by the median cross correlation with EM data. We note that this filtering criteria is subjective, and developing a Bayesian method to objectively weigh different restraints for filtering remains an interesting future development in integrative modeling.

The current filtering procedure involves three steps. In the first step, we look through the `stat.*.out` files to write out the cross correlation with EM data for each model, which, in this case, is labeled column `3`, `GaussianEMRestraint_None_CCC`. In other applications, the column that corresponds to each type of experimental data may change, depending on the scoring terms for each model. For each snapshot, a new file is written with this data (`{state}_{time}_stat.txt`).

\code{.py}
# state_dict - universal parameter
state_dict = {'0min': 3, '1min': 3, '2min': 1}
# current directory
main_dir = os.getcwd()

# 1 calling extracting_stat_files function and related parameters
keys_to_extract = [3]
runs_nr = 50
replica_nr = 16
replica_output_name = 'output'
decimals_nr = 16

extracting_stat_files(state_dict, runs_nr, replica_nr, replica_output_name, keys_to_extract, decimals_nr)
print("extracting_stat_files is DONE")
print("")
print("")
\endcode

In the second step, we want to determine the median value of EM cross correlation for each snapshot. We wrote `general_rule_calculation` to look through the `general_rule_column` for each `{state}_{time}_stat.txt` file and determine both the median value and the number of structures generated.

\code{.py}
# 2 calling general_rule_calculation and related parameters
general_rule_column = '3'

general_rule_calculation(state_dict, general_rule_column)

print("general_rule_calculation is DONE")
print("")
print("")
\endcode

In the third step, we use the `imp_sampcon select_good` tool to filter each snapshot, according to the median value determined in the previous step. For each snapshot, this function produces a file, `good_scoring_models/model_ids_scores.txt`, which contains the run, replicaID, scores, and sampleID for each model that passes filtering. It also saves RMF files with each model from two independent groups of sampling runs from each snapshot to `good_scoring_models/sample_A` and `good_scoring_models/sample_B`, writes the scores for the two independent groups of sampling runs to `good_scoring_models/scoresA.txt` and `good_scoring_models/scoresB.txt`, and writes `good_scoring_models/model_sample_ids.txt` to connect each model to its division of sampling run. More information on `imp_sampcon` is available in the analysis portion of the [actin tutorial](https://integrativemodeling.org/tutorials/actin/analysis.html).

\code{.py}
# 3 calling general_rule_filter_independent_samples
general_rule_filter_independent_samples(state_dict, main_dir)
print("general_rule_filter_independent_samples is DONE")
print("")
print("")
\endcode

## Plotting data, clustering models, and determining sampling precision {#snapshot_sampling_precision}

Next, scores can be plotted for analysis. Here, we wrote the `create_histograms` function to run `imp_sampcon plot_score` so that it plots distributions for various scores of interest. Each of these plots are saved to `histograms{state}_{time}/{score}.png`, where score is an object listed in the `score_list`. These plots are useful for debugging the modeling protocol, and should appear roughly Gaussian.

\code{.py}
# 4 calling create_histograms and related parameters
score_list = [
    'Total_Score',
    'ConnectivityRestraint_Score',
    'ExcludedVolumeSphere_Score',
    'GaussianEMRestraint_None',
    'GaussianEMRestraint_None_CCC'
] # list of histograms we want to create in each histograms{state}_{time} directory

create_histograms(state_dict, main_dir, score_list)
print("create_histograms is DONE")
print("")
print("")
\endcode

We then check the number of models in each sampling run though our function, `count_rows_and_generate_report`, which writes the `independent_samples_stat.txt` file. Empirically, we have found that ensuring the overall number of models in each independent sample after filtering is roughly equal serves a good first check on sampling convergence.

\code{.py}
# 5 calling count_rows_and_generate_report
count_rows_and_generate_report(state_dict)
print("count_rows_and_generate_report is DONE")
print("")
print("")
\endcode

Next, we write the density range dictionaries, which are output as `{state}_{time}_density_ranges.txt`. These dictionaries label each protein in each snapshot model, which will be passed into `imp_sampcon` to calculate the localization density of each protein.

\code{.py}
# 6 calling create_density_dictionary:
create_density_dictionary_files(state_dict, main_dir)
print("create_density_dictionary is DONE")
print("")
print("")
\endcode

Next, we run `imp_sampcon exhaust` on each snapshot. This code performs checks on the exhaustiveness of the sampling. Specifically it analyzes the convergence of the model score, whether the two model sets were drawn from the same distribution, and whether each structural cluster includes models from each sample proportionally to its size. The output for each snapshot is written out to the `exhaust_{state}_{time}` folder.

\code{.py}
# 7 calling exhaust
exhaust(state_dict, main_dir)
print("exhaust is DONE")
print("")
print("")
\endcode

Plots for determining the sampling precision are shown below for a single snapshot, 1_2min. (a) Tests the convergence of the lowest scoring model (`snapshot_{state}_{time}.Top_Score_Conv.pdf`). Error bars represent standard deviations of the best scores, estimated by selecting different subsets of models 10 times. The light-blue line indicates a lower bound reference on the total score. (b) Tests that the scores of two independently sampled models come from the same distribution (`snapshot_{state}_{time}.Score_Dist.pdf`). The difference between the two distributions, as measured by the KS test statistic (D) and KS test p-value (p) indicates that the difference is both statistically insignificant (p>0.05) and small in magnitude (D<0.3). (c) Determines the structural precision of a snapshot model (`snapshot_{state}_{time}.ChiSquare.pdf`). RMSD clustering is performed at 1 Å intervals until the clustered population (% clustered) is greater than 80%, and either the χ<sup>2</sup> p-value is greater than 0.05 or Cramer’s V is less than 0.1. The sampling precision is indicated by the dashed black line. (d) Populations from sample 1 and sample 2 are shown for each cluster (`snapshot_{state}_{time}.Cluster_Population.pdf`).

\image html Snapshot_Exhaust.png width=1200px

Further structural analysis can be calculated by using the `cluster.*` files. The `cluster.*.{sample}.txt` files contain the model number for the models in that cluster, where `{sample}` indicates which round of sampling the models came from. The `cluster.*` folder contains an RMF for centroid model of that cluster, along with the localization densities for each protein. The localization densities of each protein from each independent sampling can be compared to ensure independent samplings produce the same results.

Ideally, each of these plots should be checked for each snapshot. As a way to summarize the output of these checks, we can gather the results of the KS test and the sampling precision test for all snapshots. This is done by running `extract_exhaust_data` and `save_exhaust_data_as_png`, which write `KS_sampling_precision_output.txt` and `KS_sampling_precision_output.png`, respectively.

\code{.py}
# 8 calling extract_exhaust_data
extract_exhaust_data(state_dict)
print("extract_exhaust_data is DONE")
print("")
print("")

# 9 calling save_exhaust_data_as_png
save_exhaust_data_as_png()
print("save_exhaust_data_as_png is DONE")
print("")
print("")
\endcode

These codes write a table that include the KS two sample test statistic (D), the KS test p-value, and the sampling precision for each snapshot model, which is replotted below.

\image html Snapshot_sampling.png width=600px

## Visualizing models {#snapshot_visualization}

The resulting RMF files and localization densities from this analysis can be viewed in [UCSF Chimera](https://www.rbvi.ucsf.edu/chimera/) (version>=1.13) or [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/).

Here, we plotted each centroid model (A - blue, B - orange, and C - purple) from the most populated cluster for each snapshot model and compared that model to the experimental EM profile (gray).

\image html static_snapshots_noCC.png width=600px

Finally, now that snapshot models have been assessed, we can perform [Modeling of trajectories.] (@ref trajectories)
