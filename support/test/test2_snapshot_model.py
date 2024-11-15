#!/usr/bin/env python
import unittest
import os
import subprocess
import IMP
import pickle
import IMP.test
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

def run_sim(state, time, nframes, shuffle=True):
    """
    This function creates snapshot models (.rmf frames) using pmi. More about model representation with pmi and running
    ReplicaExchange can be found here:
    https://integrativemodeling.org/tutorials/actin/pmidesign.html
    """
    # Topology file
    topology_file = os.path.join(TOPDIR, 'modeling', 'Heterogeneity', 'Heterogeneity_Modeling', f"{state}_{time}_topol.txt")
    # Paths to input data for topology file
    pdb_dir = os.path.join(TOPDIR, 'modeling', 'Input_Information', 'PDB')
    fasta_dir = os.path.join(TOPDIR, 'modeling', 'Input_Information', 'FASTA')
    # Path where forward gmms are created with BuildSystem (based ont topology file)
    # Path to experimental gmms
    exp_gmm_dir = os.path.join(TOPDIR, 'modeling', 'Input_Information', 'ET_data', 'add_noise')

    # Setting model
    mdl = IMP.Model()

    # Read the topology file
    # t = IMP.pmi.topology.TopologyReader(topology_file, pdb_dir=pdb_dir, fasta_dir=fasta_dir, gmm_dir=forward_gmm_dir)
    t = IMP.pmi.topology.TopologyReader(topology_file, pdb_dir=pdb_dir, fasta_dir=fasta_dir)

    # Create a system from a topology file. Resolution is set on 1.
    bs = IMP.pmi.macros.BuildSystem(mdl, resolutions=1, name=f'Static_snapshots_{state}_{time}')
    bs.add_state(t)

    #  Macro execution: It gives hierarchy and degrees of freedom (dof).
    # In dof we define how much can each (super) rigid body translate and rotate between two adjacent Monte Carlo steps
    root_hier, dof = bs.execute_macro(max_rb_trans=1.0,
                                      max_rb_rot=0.5, max_bead_trans=2.0,
                                      max_srb_trans=1.0, max_srb_rot=0.5)

    # Adding Restraints
    # Empty list where the data from restraints should be collected
    output_objects = []

    # Two common restraints: ConnectivityRestraint and ExcludedVolumeSphere
    # ConnectivityRestraint is added for each "molecule" separately
    for m in root_hier.get_children()[0].get_children():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
        cr.add_to_model()
        output_objects.append(cr)

    # Add excluded volume
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[root_hier],
        resolution=1000)
    output_objects.append(evr)

    # Applying time-dependent EM restraint. Point to correct gmm / mrc file at each time point
    # Path to corresponding .gmm file (and .mrc file)
    em_map = exp_gmm_dir + f"/{time}_noisy.gmm"

    # Create artificial densities from hierarchy
    densities = IMP.atom.Selection(root_hier,
                                   representation_type=IMP.atom.DENSITIES).get_selected_particles()

    # Create EM restraint based on these densities
    emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,
        target_fn=em_map,  #
        slope=0.000001,
        scale_target_to_mass=True,
        weight=1000)
    output_objects.append(emr)

    # Generate random configuration
    if shuffle:
        IMP.pmi.tools.shuffle_configuration(root_hier,
                                        max_translation=50)

    # Add EM restraint and excluded volume restraint to the model
    evr.add_to_model()
    emr.add_to_model()

    # Perform replica exchange sampling
    rex = IMP.pmi.macros.ReplicaExchange(mdl,
                                         root_hier=root_hier,
                                         monte_carlo_sample_objects=dof.get_movers(),
                                         global_output_directory='output',
                                         output_objects=output_objects,
                                         monte_carlo_steps=200,
                                         number_of_best_scoring_models=0,
                                         number_of_frames=nframes)
    return rex

class static_snapshots(unittest.TestCase):
    def test_modeling_script(self):
        """Test the main modeling script runs with 3_1min"""
        with IMP.test.temporary_directory() as tmpdir:
            os.chdir(tmpdir)
            rep=run_sim(state='3', time='1min', nframes=1)
            rep.execute_macro()

    def test_score(self):
        """Test that independendent calls to the modeling script give the same score"""
        with IMP.test.temporary_directory() as tmpdir:
            # Set up modeling but don't run sampling
            os.chdir(tmpdir)
            rep=run_sim(state='3', time='1min', nframes=1, shuffle=False)

            # Run the original ReplicaExchange and get the final score
            IMP.random_number_generator.seed(99)
            rep.execute_macro()
            rs = IMP.pmi.tools.get_restraint_set(rep.model)
            original_score = rs.evaluate(False)
            del rep, rs

            # With the same random seed, we should get the exact same trajectory
            IMP.random_number_generator.seed(99)
            newrep = run_sim(state='3', time='1min', nframes=1, shuffle=False)
            newrep.execute_macro()
            rs = IMP.pmi.tools.get_restraint_set(newrep.model)
            new_score = rs.evaluate(False)
            self.assertAlmostEqual(original_score, new_score, delta=1e-4)


if __name__ == '__main__':
    unittest.main()