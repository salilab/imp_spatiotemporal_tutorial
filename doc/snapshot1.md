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

To begin, we built a topology file for the complete system, `spatiotemporal_topology.txt`. Based on our observation of the structure of the complex, we chose to represent each protein with at least 2 separate rigid bodies, and left the first 28 residues of protein C as flexible beads. Rigid bodies were described with 1 bead for every residue, and 10 residues per Gaussian. Flexible beads were described with 1 bead for every residue and 1 residue per Gaussian. A more complete description of the options available in topology files is available in the the [TopologyReader](@ref IMP::pmi::topology::TopologyReader) documentation.

\code{.txt}
|molecule_name | color | fasta_fn | fasta_id | pdb_fn | chain | residue_range | pdb_offset | bead_size | em_residues_per_gaussian | rigid_body | super_rigid_body | chain_of_super_rigid_bodies | 

|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|-1,18|2|1|10|1|1||
|Ubi-E2-D3|blue|3rpg.fasta.txt|Ubi-E2-D3|3rpg.pdb|A|19,147|2|1|10|2|1||
|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|3,83|-2|1|10|3|2||
|BMI-1|red|3rpg.fasta.txt|BMI-1|3rpg.pdb|B|84,101|-2|1|10|4|2||
|E3-ubi-RING2|green|3rpg.fasta.txt|E3-ubi-RING2|BEADS|C|16,44|-15|1|1|5|3||
|E3-ubi-RING2|green|3rpg.fasta.txt|E3-ubi-RING2|3rpg.pdb|C|45,116|-15|1|10|6|3||
\endcode

## Scoring the model

Restraints:
Connectivity
Excluded volume
Cross-linking
EM restraint (Representation of EM data)

## Searching for good scoring models

Sampling

Describe sampling output

# Generalizing modeling to all snapshots

Next, we will describe the process of modeling a multiple static snapshots, as performed by running `start_sim.py`.

Selection of snapshots to model

Run loop above for all snapshots