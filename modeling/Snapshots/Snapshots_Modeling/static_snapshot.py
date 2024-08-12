'''
This code creates snapshot models (.rmf frames) using pmi. More about model representation with pmi and running
ReplicaExchange can be found here:
https://integrativemodeling.org/tutorials/actin/pmidesign.html
'''
import IMP
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.dof
import ihm.cross_linkers
import IMP.atom
import IMP.saxs
import os
import sys

# Running parameters to access correct path of ET_data for EM restraint and topology file for certain {state}_{time}_topol.txt
state = sys.argv[1]
time = sys.argv[2]

# Topology file
topology_file = f"../{state}_{time}_topol.txt"
# Paths to input data for topology file
pdb_dir = "../../Input_Information/PDB"
fasta_dir = "../../Input_Information/FASTA"
# Path where output .gmm created with BuildSystem (based ont topology file) and execute_macro are saved for each snapshot{state}_{time}
gmm_dir = "../ET_data/"

# Setting model
mdl = IMP.Model()

# Read the topology file
t = IMP.pmi.topology.TopologyReader(topology_file, pdb_dir=pdb_dir, fasta_dir=fasta_dir, gmm_dir=gmm_dir)


# Create a system from a topology file. Resolution is set on 1.
bs = IMP.pmi.macros.BuildSystem(mdl, resolutions= 1, name= f'Static_snapshots_{state}_{time}')
bs.add_state(t)

#  Macro execution: It gives hierarchy and degrees of freedom (dof).
# In dof we define how much can each (super) rigid body translate and rotate between two adjacent Monte Carlo steps
root_hier, dof = bs.execute_macro(max_rb_trans=1.0,
                                  max_rb_rot=0.5, max_bead_trans=2.0,
                                  max_srb_trans=1.0, max_srb_rot=0.5)

# Adding Restraints
# Empty list where the data from restraints should be collected
output_objects=[]

# Two common restraints: ConnectivityRestraint and ExcludedVolumeSphere
# ConnectivityRestraint is added for each "molecule" separately
for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)


evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                     included_objects=[root_hier],
                                     resolution=1000)
output_objects.append(evr)

# Applying time-dependent EM restraint
# Path to corresponding .gmm file (and .mrc file)
em_map = f"../ET_data/experimental/{time}_noisy.gmm" #Path to .gmm file

# Create artificial densities from hierarchy
densities = IMP.atom.Selection(root_hier,
                 representation_type=IMP.atom.DENSITIES).get_selected_particles()

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,
        target_fn=em_map, #
        slope=0.000001,
        scale_target_to_mass=True,
        weight=1000)
output_objects.append(emr)


# ReplicaExchange and sampling
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50)  #Is this optimal?

evr.add_to_model()
emr.add_to_model()



rex=IMP.pmi.macros.ReplicaExchange(mdl,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory='output', # name 'output' is the best for imp sampcon select_good
        output_objects=output_objects,
        monte_carlo_steps=200, # Number of MC steps between writing frames.
        number_of_best_scoring_models=0,
        number_of_frames=500)
# It is optimal for pmi analysis that up to cca 30000 frames are generated otherwise analysis may take to long
# In our case, for each snapshot we generated 25000 frames altogether (50*500)
rex.execute_macro()

