"""
This code plots the heterogeneity of the system.
"""

import numpy as np
import matplotlib.pyplot as plt

# colors
blue1=np.array([0,0,128])/255
grey=np.array([188,  188,  188])/255


if __name__ == "__main__":
    # Import experimental data
    expA=np.loadtxt("../../Input_Information/gen_FCS/A_v2.txt")
    expB=np.loadtxt("../../Input_Information/gen_FCS/B_v2.txt")
    expC=np.loadtxt("../../Input_Information/gen_FCS/C_v2.txt")

    # Import heterogeneity models
    model0=np.loadtxt("../Heterogeneity_Modeling/0min.txt")
    model1=np.loadtxt("../Heterogeneity_Modeling/1min.txt")
    model2=np.loadtxt("../Heterogeneity_Modeling/2min.txt")

    # Convert heterogeneity models from copy number as a function of time to copy number as a function of protein
    # There are 7 total models. First row is the copy number of that model. 2nd row is the time of that model
    # For protein A
    modelA=np.zeros((7,2))
    modelA[0:3,0] = model0[:,0]
    modelA[0:3,1] = 0
    modelA[3:6, 0] = model1[:, 0]
    modelA[3:6, 1] = 1
    modelA[6, 0] = model2[0]
    modelA[6, 1] = 2

    # For protein B
    modelB = np.zeros((7, 2))
    modelB[0:3, 0] = model0[:, 1]
    modelB[0:3, 1] = 0
    modelB[3:6, 0] = model1[:, 1]
    modelB[3:6, 1] = 1
    modelB[6, 0] = model2[1]
    modelB[6, 1] = 2

    # For protein C
    modelC = np.zeros((7, 2))
    modelC[0:3, 0] = model0[:, 2]
    modelC[0:3, 1] = 0
    modelC[3:6, 0] = model1[:, 2]
    modelC[3:6, 1] = 1
    modelC[6, 0] = model2[2]
    modelC[6, 1] = 2

    # Plot data
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(5,8))
    ax1.errorbar(expA[:,0], expA[:,1], yerr=expA[:,2], color=grey)
    ax1.plot(modelA[:,1], modelA[:,0], 'o', color=blue1)
    ax1.set_ylabel('Copy number')
    ax1.set_xlabel('Time [min]')
    ax1.set_xlim(-0.2,2.7)
    ax1.set_ylim(-1,2.5)
    ax1.legend(['model', 'experiment'], loc='upper left')
    ax2.errorbar(expB[:,0], expB[:,1], yerr=expB[:,2], color=grey)
    ax2.plot(modelB[:,1], modelB[:,0], 'o', color=blue1)
    ax2.set_ylabel('Copy number')
    ax2.set_xlabel('Time [min]')
    ax2.set_xlim(-0.2, 2.7)
    ax2.set_ylim(-0.5, 3)
    ax2.legend(['model', 'experiment'], loc='upper left')
    ax3.errorbar(expC[:,0], expC[:,1], yerr=expC[:,2], color=grey)
    ax3.plot(modelC[:,1], modelC[:,0], 'o', color=blue1)
    ax3.set_ylabel('Copy number')
    ax3.set_xlabel('Time [min]')
    ax3.set_xlim(-0.2, 2.7)
    ax2.set_ylim(-0.5, 3)
    ax3.legend(['model', 'experiment'], loc='upper left')
    plt.subplots_adjust(hspace=0.5)
    plt.show()

