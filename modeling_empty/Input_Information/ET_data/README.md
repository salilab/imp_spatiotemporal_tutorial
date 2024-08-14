### This directory code necessary for making synthetic electron tomography (ET) density maps from PDB structures.

#### Creation of ET data involved 3 steps:
1. write_mrc.py - script that converts a structure (PDB format) into a density map
2. make_gmm/make_gmm.sh - script that fits each density map to a gmm and writes a new density map. CC between the new and old map can be evaluated by calc_CC.py
3. add_noise/add_noise.py - script that takes a gmm, adds noise, and returns a new, noisy gmm. CC between the new and old map can be evaluated by calc_CC.py
