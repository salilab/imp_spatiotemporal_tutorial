### This directory code necessary for making synthetic fluorescence correlation spectroscopy (FCS) data.

#### FCS data is created in 3 steps, performed by gen_FCS.py:
1. Take 10 samples from a Gausian distribution with constant standard deviation and 
a mean given by a tanh function of based on time (For example: A.txt)
2. Compute the mean and standard deviation of these 10 samples (For example: A_v2.txt)
3. Prepare this data for IMP by extracting time points of interest 
(For example: exp_compA.txt)

This data can be plotted by plot_FCS.m.
