"""
Script to add noise to a gmm.
Takes 3 inputs: 1) old gmm file (string), 2) new gmm file (string), 
and 3) new mrc file (string)
"""
import sys
import os
import math
import numpy as np
import random
import IMP
import IMP.isd.gmm_tools as gmm_tools

# old file- original gmm
old_file=sys.argv[1]
# new file - new gmm
new_file=sys.argv[2]
# new file2 - new mrc
new_file2=sys.argv[3]

old=open(old_file,'r')
new=open(new_file,'w')

def change_mean(gmm_line,noise=0.1):
"""
Changes the mean values of the gmm.
@param gmm_line: a line of a gmm file (string)
@param noise: level of noise to add to the gmm. 
Each mean value is multiplied by 1+/-noise (float)
"""
    # split means column of the gmm
    means=gmm_line[3].split(' ')
    new_means=[]
    for mean in means:
        # add noise to means
        new_means.append(float(mean)*random.uniform(1-noise,1+noise))
    # replace portion of the gmm line
    gmm_line[3]=str(new_means[0])+' '+str(new_means[1])+' '+str(new_means[2])
    return gmm_line

line=old.readline()
# loop over the old file
while line:
    if len(line)>0:
        # make sure not commented
        if line[0]=='#':
            new_line=line
        # edit the lines that correspond to the gmm
        else:
            line_split=line.split('|')
            line_split=change_mean(line_split)
            new_line = '|'
            new_line=new_line.join(line_split)
    # write to new file
    new.write(new_line)
    line=old.readline()
old.close()
new.close()

# read in gmm and add noise
model=IMP.Model()
ps=[]
gmm_tools.decorate_gmm_from_text(new_file,ps,model)
gmm_data = [IMP.core.Gaussian(p) for p in ps]
# write noisy gmm to a new file
gmm_tools.write_gmm_to_map(gmm_data,new_file2,5.0)


