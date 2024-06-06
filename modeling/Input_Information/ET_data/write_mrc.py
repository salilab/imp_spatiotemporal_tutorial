"""
   Code to compute synthetic electron tomography (ET) data 
   for the spatiotemporal integrative modeling tutorial
   
   Takes 2 inputs: 1) filename for the input pdb file 2) filename of the output MRC file
"""
import sys
import os
import math
import numpy as np
import IMP
import IMP.core
import IMP.em

def write_mrc(sim_file,mrc_file,MRCresolution=10.0,voxel=5.0):
"""
	Function to write an MRC file from a PDB structure
	
	@param sim_file: PDB file name of a protein structure (string)
	@param mrc_file: MRC file to write (string)
	@param MRCresolution: half width the Gaussian (float)
	@param voxel: voxel size, in Angstroms (float)

	Writes an MRC file to mrc_file and returns the density map.
    """
    # Read in RMF file
    model = IMP.Model()
    heir=IMP.atom.read_pdb(sim_file,model)
    print('Number of leafs:')
    print(len(IMP.atom.get_leaves(heir)))
    ps=IMP.atom.get_leaves(heir)
    # calculate density map
    dmap = IMP.em.SampledDensityMap(ps, MRCresolution, voxel)
    dmap.calcRMS()
    dmap.set_was_used(True)
    IMP.em.write_map(dmap, mrc_file, IMP.em.MRCReaderWriter())
    return dmap
    
# Input filenames
pdb_file=sys.argv[1]
sim_mrc_file=sys.argv[2]

# call main function
model_density=write_mrc(pdb_file,sim_mrc_file)
print('Wrote density map!')
