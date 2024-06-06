"""
Code to calculate CC coefficient between two sets of files. These are hard-set by end1 and end2.
"""
import sys
import os
import math
import numpy as np
import IMP
import IMP.em

times=['0min','1min','2min']

end1='_fitted.mrc'
end2='_gmm.mrc'

cc_list=[]
for time in times:
    exp_density1 = IMP.em.read_map(time+end1, IMP.em.MRCReaderWriter())
    exp_density2 = IMP.em.read_map(time+end2, IMP.em.MRCReaderWriter())
    cc1 = IMP.em.get_coarse_cc_coefficient(exp_density1, exp_density2, 0, True)
    cc_list.append(cc1)

new_file='cc.txt'
new=open(new_file,'w')
new.write('CC between '+end1+' and '+end2+'\n')
for i in range(len(cc_list)):
    new.write(times[i]+'\t'+str(cc_list[i])+'\n')
new.close()
