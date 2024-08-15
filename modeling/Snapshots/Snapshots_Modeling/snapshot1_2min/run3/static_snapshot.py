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

# Running parameters will be mainly used to access correct path
state = sys.argv[1]
time = sys.argv[2]


## Topology Filespatiotemporal_topology.txt
topology_file = f"../{state}_{time}_topol.txt"
pdb_dir = "../../PDB/"
fasta_dir = "../../FASTA/"
gmm_dir = "../ET_data/" #Outputs will be saved here

#Model representation

mdl = IMP.Model()

## Read the topology file
t = IMP.pmi.topology.TopologyReader(topology_file, pdb_dir=pdb_dir, fasta_dir=fasta_dir, gmm_dir=gmm_dir)


## Create a BuildSystem from a topology file. Resolution is set on 1. Then for each timepoint one state can be created. (so far- we have only one state for final model)
bs = IMP.pmi.macros.BuildSystem(mdl, resolutions= 1, name= f'Static_snapshots_{state}_{time}')
bs.add_state(t)

## Macro execution. It gives us hierarchy and degrees of freedom (dof).
## In dof we define how much can each (super) rigid body translate and rotate during one step of sampling

root_hier, dof = bs.execute_macro(max_rb_trans=1.0,
                                  max_rb_rot=0.5, max_bead_trans=2.0,
                                  max_srb_trans=1.0, max_srb_rot=0.5)
#Lets see how it works with parameters

#Also we created list of molecules:
molecules = t.get_components()
print(molecules)


#Adding Restraints

##Two common restraints

output_objects=[] #Here all the data from restraints should be collected

#Completely new approach of adding restraints:
for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print(output_objects)


evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                     included_objects=[root_hier],
                                     resolution=1000)  #based on what do you set this value - lowest possible resolution -> very expensive restraint
output_objects.append(evr) #Also this part was added

##"Scoring restraint" is time depandent ET
## First we need to create densities from model representation (is this even necessary??)


em_map = f"../ET_data/experimental/{time}_noisy.gmm" #Path to gmm file
# SECOND CHANGE! -> ALSO CHANGE EACH GMM FILE TO NOISY INSTEAD OF FITTED

densities = IMP.atom.Selection(root_hier,
                 representation_type=IMP.atom.DENSITIES).get_selected_particles()

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,        # We compare artificially created densities from model representation to EM map approximated as GMM
        target_fn=em_map, #
        slope=0.000001, # Do I leave it as a default? 0.00000001, I remove 2 zeros..Is it this slope too steep or too gentle?
        scale_target_to_mass=True, # Normalizes the mass of the model wrs: EM map
        weight=1000)       #I guess it is ok like that! -> is it even reevant to
output_objects.append(emr)

print(output_objects)

#Sampling
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50)  #Is this optimal?

dof.optimize_flexible_beads(500)  #Also this part was added

evr.add_to_model()
emr.add_to_model()



rex=IMP.pmi.macros.ReplicaExchange(mdl,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory='output',                     # We set the output directory for this run.
        output_objects=output_objects,                          # Write these items to the stat file
        monte_carlo_steps=200,                              # Number of MC steps between writing frames. Is that ok?
        number_of_best_scoring_models=0,                  # set >0 to store best scoring PDB files
        number_of_frames=500)                             # We also lower this
rex.execute_macro()

