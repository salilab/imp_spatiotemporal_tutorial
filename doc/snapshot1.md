Snapshot modeling steps 2-4: representation, scoring, and search process {#snapshot1}
====================================

Here, we describe how to build models of static snapshots using IMP. We note that this process is similar to previous tutorials of [actin](https://integrativemodeling.org/tutorials/actin/) and [RNA PolII](https://integrativemodeling.org/tutorials/rnapolii_stalk/).

Navigate to the `Snapshots/Snapshots_Modeling/` folder. Here, you will find two python scripts. The first, `static_snapshot.py`, uses IMP to represent, score, and search for models of a single static snapshot. The second, `start_sim.py`, automates the creation of multiple static snapshots at each time point, which score well according to protein copy number data.

# Modeling one snapshot

Here, we will describe the process of modeling a single snapshot model, as performed by running `static_snapshot.py`.

## Representing the model

The second step in integrative modeling is representing the data and the model. In general, the *representation* of a system is defined by all the variables that need to be determined.

For our model of a protein complex, we use a combination of two representations. The first is a series of *spherical beads*, which can correspond to portions of the biomolecules of interest, such as atoms or groups of atoms. The second is a series of *3D Gaussians*, which help calculate the overlap between our model and the density from ET data.

Beads and Gaussians in our model belong to either a *rigid body* or *flexible string*. The positions of all beads and Gaussians in a single rigid body are constrained during sampling and do not move relative to each other. Meanwhile, flexible beads can move freely during sampling, but are restrained by sequence connectivity.

To begin, we built a topology file with the representation for the model of the complete system, `spatiotemporal_topology.txt`. Later, we will use this complete topology as a template to build topologies of each snapshot. Based on our observation of the structure of the complex, we chose to represent each protein with at least 2 separate rigid bodies, and left the first 28 residues of protein C as flexible beads. Rigid bodies were described with 1 bead for every residue, and 10 residues per Gaussian. Flexible beads were described with 1 bead for every residue and 1 residue per Gaussian. A more complete description of the options available in topology files is available in the the [TopologyReader](@ref IMP::pmi::topology::TopologyReader) documentation.

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

## Scoring the model

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

## Searching for good scoring models

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

# Generalizing modeling to all snapshots

Next, we will describe the process of modeling a multiple static snapshots, as performed by running `start_sim.py`.

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

\code{.py}
# 2a - parameters for generate_all_snapshots
# state_dict - universal parameter
state_dict = {'0min': 3, '1min': 3, '2min': 1}

main_dir = os.getcwd()
items_to_copy = ['static_snapshot.py']  # additionally we need to copy only specific topology file
# jobs script will likely depend on the user's cluster / configuration
job_template = ("#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -N Tutorial\n#$ -pe smp 16\n"
                "#$ -l h_rt=48:00:00\n\nmodule load Sali\nmodule load imp\nmodule load mpi/openmpi-x86_64\n\n"
                "mpirun -np $NSLOTS python3 static_snapshot.py {state} {time}")
number_of_runs = 50

# 2b - calling generate_all_snapshots
generate_all_snapshots(state_dict, main_dir, items_to_copy, job_template, number_of_runs)
\endcode

Finally, we note that sometimes errors such as the one below can arise during sampling. These errors are caused by issues generating forward GMM files, which is done stochastically. If such issues arrise, remove all files in the `forward_densities` folder for that snapshot and resubmit the corresponding jobs.

\code{.py}
  File "/imp/main/20240607-af6f9d6a95/lib/release8/IMP/isd/gmm_tools.py", line 35, in decorate_gmm_from_text
    weight = float(fields[2])
IndexError: list index out of range
\endcode

Finally, once sampling is complete, we aim to [assess the snapshot models.] (@ref snapshot_assess)