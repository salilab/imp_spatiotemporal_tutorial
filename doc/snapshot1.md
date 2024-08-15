Snapshot modeling steps 2-4: representation, scoring, and search process {#snapshot1}
====================================

Here, we describe how to build models of static snapshots using IMP. We note that this process is similar to previous tutorials (https://integrativemodeling.org/tutorials/actin/ and https://integrativemodeling.org/tutorials/rnapolii_stalk/).

Navigate to the `Snapshots/Snapshots_Modeling/` folder. Here, you will find two python scripts. The first, `static_snapshot.py`, uses IMP to represent, score, and search for models of a single static snapshot. The second, `start_sim.py`, automates the creation of multiple static snapshots at each time point, which score well according to protein copy number data.

# Modeling one snapshot

Here, we will describe the process of modeling a single snapshot model, as performed by running `static_snapshot.py`.

## Representing the model

The second step in integrative modeling is representing the data and the model. In general, the *representation* of a system is defined by all the variables that need to be determined.

For our model of a protein complex, we use a combination of two representations. The first is a series of *spherical beads*, which can correspond to portions of the biomolecules of interest, such as atoms or groups of atoms. The second is a series of *3D Gaussians*, which help calculate the overlap between our model and the density from ET data.

Beads and Gaussians in our model belong to either a *rigid body* or *flexible string*. The positions of all beads and Gaussians in a single rigid body are constrained during sampling and do not move relative to each other. Meanwhile, flexible beads can move freely during sampling, but are restrained by sequence connectivity.

To begin, we built a topology file for the complete system, `spatiotemporal_topology.txt`. Later, we will use this complete topology as a template to build topologies of each snapshot. Based on our observation of the structure of the complex, we chose to represent each protein with at least 2 separate rigid bodies, and left the first 28 residues of protein C as flexible beads. Rigid bodies were described with 1 bead for every residue, and 10 residues per Gaussian. Flexible beads were described with 1 bead for every residue and 1 residue per Gaussian. A more complete description of the options available in topology files is available in the the [TopologyReader](@ref IMP::pmi::topology::TopologyReader) documentation.

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

Sampling

Describe sampling output

# Generalizing modeling to all snapshots

Next, we will describe the process of modeling a multiple static snapshots, as performed by running `start_sim.py`.

Selection of snapshots to model

Run loop above for all snapshots