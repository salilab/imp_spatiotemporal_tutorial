"""
   Code to compute synthetic fluorescence correlation spectroscopy (FCS) data 
   for the spatiotemporal integrative modeling tutorial
"""

import random
import numpy as np
import sys
import os
import pandas as pd


def make_random_samples():
	"""
	Function to make random samples from a tanh distribution.
	The number of samples, starting time, ending time, step size, and 
	tanh parameters are set in the first few lines of the function.
	
	Returns samples of 3 different proteins, which are also saved to a file.
    """
    # starting and stopping time (in minutes)
    start=0
    end=2.6
    step=0.1
    time=np.arange(start,end,step)
    N=len(time)
    # Number of independent samples to generate
    samples=10
    # standard deviation of the experiment
    sigma=0.5
    # smoothness variable of tanh function
    eta=5
    # shift variable for tanh function
    shift=-0.2

    # initiate output files
    A=np.zeros((N,samples+1))
    B=np.zeros((N,samples+1))
    C=np.zeros((N,samples+1))
    A[:,0]=time
    B[:,0]=time
    C[:,0]=time

    # for each sample, at each time point
    for i in range(samples):
        for j in range(N):
            # generate a normal random variable with a mean determined by the tanh function and a standard deviation of sigma
            A[j,i+1]=random.gauss(0.5+0.5*np.tanh(eta*(time[j]-(2+shift))), sigma)
            B[j,i+1]=random.gauss(0.5+0.5*np.tanh(eta*(time[j]-shift)), sigma)
            C[j,i+1]=random.gauss(0.5+0.5*np.tanh(eta*(time[j]-(1+shift))), sigma)
    # save samples to text
    np.savetxt('A.txt',A)
    np.savetxt('B.txt',B)
    np.savetxt('C.txt',C)
    return A,B,C

def compute_sample_mean_and_std(A,B,C):
	"""
	Function to convert sampled distributions into mean / 
	standard deviation of the distribution. Assumes the first column is the time, 
	and subsequent columns are different copy number samples
	
	@param A: distribution for protein A (numpy matrix)
	@param B: distribution for protein B (numpy matrix)
	@param C: distribution for protein A (numpy matrix)
	
	Returns the mean and standard deviation of the distribution, 
	which are also saved to files.
    """
    N = len(A[:,0])
    # initiate A. Calculate mean and standard deviation
    A_v2=np.zeros((N,3))
    A_v2[:,0]=A[:,0]
    A_v2[:,1]=np.mean(A,axis=1)
    A_v2[:,2]=np.std(A,axis=1)
    # initiate B. Calculate mean and standard deviation
    B_v2=np.zeros((N,3))
    B_v2[:,0]=B[:,0]
    B_v2[:,1]=np.mean(B,axis=1)
    B_v2[:,2]=np.std(B,axis=1)
    # initiate C. Calculate mean and standard deviation
    C_v2=np.zeros((N,3))
    C_v2[:,0]=C[:,0]
    C_v2[:,1]=np.mean(C,axis=1)
    C_v2[:,2]=np.std(C,axis=1)
    # save A, B, and C to file
    np.savetxt('A_v2.txt',A_v2,header='time\t\t\t\tmean\t\t\t\tstd')
    np.savetxt('B_v2.txt',B_v2,header='time\t\t\t\tmean\t\t\t\tstd')
    np.savetxt('C_v2.txt',C_v2,header='time\t\t\t\tmean\t\t\t\tstd')

    return A_v2,B_v2,C_v2

def prepare_for_IMP(means_stds,fn):
"""
	Function to convert mean / standard deviation of the distribution into specific means 
	/ standard deviations at the time points that will be used for integrative modeling.
	
	
	@param means_stds: mean / standard deviation of the protein copy number 
	from experiment. Assumes the first column is the time, the second column is the mean, 
	and the third column is the standard deviation.
	@param fn: name of file to which the data will be saved in an 
	IMP-ready format (string)
	
	Returns the mean and standard deviation of the distribution, 
	which are also saved to files.
    """
    # save these times
    times = [0, 1, 2]
    # initate a data frame
    CN = pd.DataFrame(index=times,columns=['Time', 'mean', 'std'])
    # Loop over all time points
    for i in range(len(means_stds[:,0])):
        if means_stds[i,0] in times:
            CN['Time'][means_stds[i,0]] = str(int(means_stds[i,0]))+'min'
            CN['mean'][means_stds[i,0]] = means_stds[i, 1]
            CN['std'][means_stds[i,0]] = means_stds[i, 2]
    CN.to_csv(fn)


# run main
S1,S2,S3=make_random_samples()
S1_v2,S2_v2,S3_v2=compute_sample_mean_and_std(S1,S2,S3)
prepare_for_IMP(S1_v2,'exp_compA.txt')
prepare_for_IMP(S2_v2,'exp_compB.txt')
prepare_for_IMP(S3_v2,'exp_compC.txt')

