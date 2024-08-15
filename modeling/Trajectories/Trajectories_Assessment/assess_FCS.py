"""
Script for the assessment of final trajectory model. This script assess model in four ways:
-calculation of temporal precision
-calculation of precision of the model
-comparison of the model to data used in modeling (EM)
-comparison of the model to data not used in modeling (SAXS, native pdb of final complex)
"""
"""import IMP
import IMP.atom
import RMF
import IMP.rmf
import IMP.algebra
import IMP.core
import IMP.pmi.macros
import shutil
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis"""
import os
import numpy as np
import matplotlib.pyplot as plt

def read_labeled_pdf(pdf_file):
    """
    Function to read in a labeled probability distribution file output by spatiotemporal.create_DAG.
    Used to determine protein copy numbers by forward_model_copy_number.
    :param pdf_file (str): sting for the path of the labeled probability distribution file output by
    spatiotemporal.create_DAG.
    :return prob_dict (dict): dictionary defining the spatiotemporal model. Each key is a state, and each value is the
    probability of that state.
    """
    # create blank dictonary to store the results
    prob_dict = {}
    # read in labeled pdf file
    old = open(pdf_file, 'r')
    line = old.readline()
    # store the path through various nodes, as well as the probability of that path
    while line:
        line_split = line.split()
        # assumes the first string is the trajectory string, the second string is the probability
        if len(line_split) > 1:
            # use # for comments
            if line_split[0]=='#':
                pass
            else:
                trj = line_split[0]
                prob = float(line_split[1])
                # store in dictionary
                prob_dict[trj] = prob
        line = old.readline()
    old.close()
    return prob_dict

def copy_number_from_state(prot_list,trj,custom_data_folder = None):
    """
    For a trajectory, returns an array of protein copy numbers as a function of time. Used by
    forward_model_copy_number().
    :param prot_list (list): list of proteins in the model. These proteins are searched for in each config file.
    :param trj (str): string defining a single trajectory.
    :param custom_data_folder (str, optional): path to custom data folder. Defaults to None, which points to '../data/'
    :return _prots (array): 2D array of protein copy numbers. The first index loops over the time,
    while the second index value loops over the protein (ordered as A, B, C).
    :return N (int): Number of time points in each trajectory.
    """
    # find folder with config files
    if custom_data_folder:
        data_folder = custom_data_folder
    else:
        data_folder = 'data/'

    # split the trajectory into a list of individual states
    state_list=trj.split('|')
    state_list=state_list[:-1]

    N = len(state_list)
    # Map from index to protein: 0 - A, 1- B, 2- C
    _prots = np.zeros((N, len(prot_list)))

    # Grab _prots from .config file
    for i in range(0, N):
        prot_file = data_folder + state_list[i] + '.config'
        to_read = open(prot_file, 'r')
        line = to_read.readline()
        while line:
            # for each line, check if the protein is in that line
            for prot_index in range(len(prot_list)):
                if prot_list[prot_index] in line:
                    _prots[i, prot_index] += 1
            line = to_read.readline()

    return _prots,N

def forward_model_copy_number(prot_list,custom_labeled_pdf=None):
    """
    Code to perform copy number analysis on each protein in the model. Writes output files where each row is ordered
    according to the time point in the model and the first column is the mean copy number, while the second column is
    the standard deviation in copy number.
    :param prot_list (list): list of proteins in the model. These proteins are searched for in each config file.
    :param custom_labeled_pdf (str, optional): path to custom labeled probability distribution file output by
    spatiotemporal.create_DAG.
    """
    # find folder with config files
    if custom_labeled_pdf:
        _labeled_pdf = custom_data_folder
    else:
        _labeled_pdf = '../Trajectories_Modeling/output/labeled_pdf.txt'

    # Read in labeled_pdf file into a dictionary. Each trajectory is listed as a dictionary,
    # with keys as the trajectory and the values as the probability of that trajectory
    prob_dict = read_labeled_pdf(_labeled_pdf)

    # Loop over the full dictionary. Create a list with 2 values:
    # 1) the probability of the state, 2) the protein copy number of that state.
    key_list = prob_dict.keys()
    prot_prob = []
    for key in key_list:
        CN,N_times = copy_number_from_state(prot_list,key)
        prot_prob.append([prob_dict[key], CN])

    # Construct the full path to the output directory
    dir_name = "forward_model_copy_number"
    full_path = os.path.join(main_dir, dir_name)
    os.makedirs(full_path, exist_ok=True)
    os.chdir(full_path)

    # Determine copy number from the prot_prob
    for index in range(len(prot_prob[0][1][0])):
        copy_number = np.zeros((N_times, 2))
        # calculate mean
        for state in prot_prob:
            for i in range(N_times):
                copy_number[i, 0] += state[0] * state[1][i][index]
        # calculate std deviation
        for state in prot_prob:
            for i in range(N_times):
                # Calculate variance
                copy_number[i, 1] += state[0] * ((state[1][i][index] - copy_number[i, 0]) ** 2)
        # Take square root to get the standard deviation
        copy_number[:, 1] = np.sqrt(copy_number[:, 1])
        # save to file
        np.savetxt('CN_prot_'+prot_list[index]+'.txt', copy_number, header='mean CN\tstd CN')


if __name__ == "__main__":
    # state_dict - universal parameter
    state_dict = {'0min': 3, '1min': 3, '2min': 1}
    # model
    #model = IMP.Model()

    # current directory
    main_dir = os.getcwd()


    # start calling codes
    ## 1 - calculation of temporal precision
    # copy all the relevant files
    #copy_files_for_data(state_dict)

    # create two independent DAGs
    expected_subcomplexes = ['A', 'B', 'C']
    exp_comp = {'A': 'exp_compA.csv', 'B': 'exp_compB.csv', 'C': 'exp_compC.csv'}
    input = "./data/"
    outputA = "../output_modelA/"
    outputB = "../output_modelB/"

    """# Output from sampling precision and model precision to be saved in united dir: analysis_output_precision
    analysis_output = "./analysis_output_precision/"
    os.makedirs(analysis_output, exist_ok=True)

    nodesA, graphA, graph_probA, graph_scoresA = spatiotemporal.create_DAG(state_dict, out_pdf=True, npaths=3,
                                                                           input_dir=input, scorestr='_scoresA.log',
                                                                           output_dir=outputA,
                                                                           spatio_temporal_rule=False,
                                                                           expected_subcomplexes=expected_subcomplexes,
                                                                           score_comp=True, exp_comp_map=exp_comp,
                                                                           draw_dag=False)

    os.chdir(main_dir)
    nodesB, graphB, graph_probB, graph_scoresB = spatiotemporal.create_DAG(state_dict, out_pdf=True, npaths=3,
                                                                           input_dir=input, scorestr='_scoresB.log',
                                                                           output_dir=outputB,
                                                                           spatio_temporal_rule=False,
                                                                           expected_subcomplexes=expected_subcomplexes,
                                                                           score_comp=True, exp_comp_map=exp_comp,
                                                                           draw_dag=False)

    ## 1 - analysis
    analysis.temporal_precision(outputA + 'labeled_pdf.txt', outputB + 'labeled_pdf.txt',
                                output_fn='.' + analysis_output + 'temporal_precision.txt')
    # here is some difficulty accessing this directory (additional dot for output_fn should be added as described above)
    os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
    print("Step 1: calculation of temporal precision IS COMPLETED")
    print("")
    print("")"""

    ## 2 - calculation of precision of the model

    # precision is calculated from .labeled_pdf.txt in Trajectories_Modeling dir
    trajectories_modeling_input_dir = "../Trajectories_Modeling/output/"

    #analysis.precision(trajectories_modeling_input_dir + 'labeled_pdf.txt', output_fn=analysis_output + 'precision.txt')

    #os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
    #print("Step 2: calculation of precision of the model IS COMPLETED")
    #print("")
    #print("")

    # 3 - comparison of the model to data used in modeling (EM)
    #exp_mrc_base_path = "../../data/ET_data/experimental"
    #ccEM(exp_mrc_base_path)

    forward_model_copy_number(expected_subcomplexes)


    ## 4 - comparison of the model to data used in modeling (SAXS, native pdb of final complex)
    # 4a - SAXS
    #convert_rmfs(state_dict, model)
    #copy_SAXS_dat_files()
    #process_foxs(state_dict)
    #print("Step 4a: SAXS validation IS COMPLETED")
    #print("")
    #print("")

    # 4b - RMSD
    #pdb_path = "../../snapshots/PDB/3rpg.pdb"
    #RMSD(pdb_path=pdb_path, custom_n_plot=20)
    #print("Step 4a: SAXS validation IS COMPLETED")
    #print("")
    #print("")




