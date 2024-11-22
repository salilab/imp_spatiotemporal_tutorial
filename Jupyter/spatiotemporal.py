#!/usr/bin/env python3


# General imports for the tutorial
import sys, os, glob, shutil
import IMP
import RMF
import IMP.rmf
from IMP.spatiotemporal import prepare_protein_library
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import numpy as np
import matplotlib.pyplot as plt

# parameters for prepare_protein_library:
times = ["0min", "1min", "2min"]
exp_comp = {'A': '../modeling/Input_Information/gen_FCS/exp_compA.csv',
            'B': '../modeling/Input_Information/gen_FCS/exp_compB.csv',
            'C': '../modeling/Input_Information/gen_FCS/exp_compC.csv'}
expected_subcomplexes = ['A', 'B', 'C']
template_topology = '../modeling/Heterogeneity/Heterogeneity_Modeling/spatiotemporal_topology.txt'
template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}
nmodels = 3

# calling prepare_protein_library
prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                template_topology=template_topology, template_dict=template_dict)

def merge_scores(fileA, fileB, outputFile):
    """
    For each function merges scoresA.txt and scoresB.txt into {state}_{time}_scores.log

    :param fileA: path to scoresA.txt
    :param fileB: path to scoresB.txt
    :param outputFile: path to output merged .log file named {state}_{time}_scores.log for each snapshot.
    This type of .log file is used in crete_DAG to generate trajectory model.
    """
    # open both files, so data can be extracted
    with open(fileA, 'r') as file_a:
        data_a = file_a.readlines()

    with open(fileB, 'r') as file_b:
        data_b = file_b.readlines()

    # Merge the content of both files
    merged_data = data_a + data_b

    # Write the merged content into the output file
    with open(outputFile, 'w') as output:
        output.writelines(merged_data)

def create_data_and_copy_files(state_dict, custom_source_dir1 = None, custom_source_dir2 = None, custom_source_dir3 = None):
    """
    Copies three types of files important to generate trajectory models:
    -.config files created with start_sim.py in Snapshot_Modeling (source_dir1)
    -time-dependent stoichiometry data for each timepoint. Data should be presented in .csv file. With this function all
    csv file in source_dir2 will be copied. These .csv files will be used in the exp_comp dictionary in create_DAG
    function
    -scoresA and scoresB for each snapshot created with imp sampcon exhaust
    (source_dir1 + snapshot + good_scoring_models) are merged into total score .txt using merge_scores helper function.
    All copied files are gathered in newly created './data/' directory, where everything is prepared for create_DAG
    function.


    :param state_dict (dict): state_dict: dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_source_dir1 (optional - str): Custom path to heterogeneity modeling dir (heterogeneity_modeling.py),
    to copy .config files
    :param custom_source_dir2 (optional - str): Custom path to stoichiometry data dir
    :param custom_source_dir3 (optional - str): Custom path to snapshot modeling dir (start_sim.py), to copy .config
    files and to access scoresA/scoresB (custom_source_dir3 + snapshot{state}_{time} + 'good_scoring_models')
    """

    # Create the destination directory if it does not exist (./data/). Here all the
    destination_dir = './data/'
    os.makedirs(destination_dir, exist_ok=True)

    # Path to heterogeneity modeling dir
    if custom_source_dir1:
        source_dir1 = custom_source_dir1
    else:
        source_dir1 = '../../Heterogeneity/Heterogeneity_Modeling/'

    # Path to stoichiometry data dir
    if custom_source_dir2:
        source_dir2 = custom_source_dir2
    else:
        source_dir2 = '../../Input_Information/gen_FCS/'

    # Path to snapshot modeling dir
    if custom_source_dir3:
        source_dir3 = custom_source_dir3
    else:
        source_dir3 = '../../Snapshots/Snapshots_Modeling/'

    # Copy all .config files from the first source directory to the destination directory
    try:
        for file_name in os.listdir(source_dir1):
            if file_name.endswith('.config'):
                full_file_name = os.path.join(source_dir1, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, destination_dir)
        print(".config files are copied")
    except Exception as e:
        print(f".config files cannot be copied. Try do do it manually. Reason for Error: {e}")

    # Copy all .csv stoichiometry files from the second source directory to the destination directory
    try:
        for file_name in os.listdir(source_dir2):
            if file_name.endswith('.csv'):
                full_file_name = os.path.join(source_dir2, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, destination_dir)
        print(".csv stoichiometry files are copied")
    except Exception as e:
        print(f".csv stoichiometry files cannot be copied. Try do do it manually. Reason for Error: {e}")

    # Copy scoresA and scoresB from the snapshot_{state}_{time} directories and first source directory path
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"snapshot{state}_{time}"
            good_scoring_path = "good_scoring_models"
            file_a = os.path.join(source_dir3, dir_name, good_scoring_path, "scoresA.txt")
            file_b = os.path.join(source_dir3, dir_name, good_scoring_path, "scoresB.txt")
            output_file = os.path.join(destination_dir, f"{state}_{time}_scores.log") # name of the output file

            try:
                # Ensure the directory exists before try to read/write files
                if os.path.exists(file_a) and os.path.exists(file_b):
                    merge_scores(file_a, file_b, output_file) # call helper function to merge files
                    print(f"Scores for snapshot{state}_{time} have been merged and saved")
                else:  # many things can go wrong here, so it is good to know where is the problem
                    print(f"Path doesn't exist: {source_dir3}")
                    print(f"Files not found in directory: {dir_name}")
                    print(f"Files not found in directory: {file_a}")
                    print(f"Files not found in directory: {file_b}")
                    print(f"Output directory doesn't exist: {destination_dir}")
            except Exception as e:
                print(f"total scores files cannot be copied of merged. Reason for Error: {e}")

# copy all the relevant files for create_DAG
# it is important that everything starts from main dir
main_dir = os.getcwd()
os.chdir(main_dir)
state_dict = {'0min': 3, '1min': 3, '2min': 1}
create_data_and_copy_files(state_dict, custom_source_dir1=main_dir, custom_source_dir2='../modeling/Input_Information/gen_FCS/', custom_source_dir3='../modeling/Snapshots/Snapshots_Modeling/')

# then trajectory model is created based on the all copied data
expected_subcomplexes = ['A', 'B', 'C']
exp_comp = {'A': 'exp_compA.csv', 'B': 'exp_compB.csv', 'C': 'exp_compC.csv'}
input = './data/'
output = "../output/"

nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(state_dict, out_pdf=True, npaths=3,
                                                                       input_dir=input, scorestr='_scores.log',
                                                                       output_dir=output, spatio_temporal_rule=True,
                                                                       expected_subcomplexes=expected_subcomplexes,
                                                                       score_comp=True, exp_comp_map=exp_comp,
                                                                       draw_dag=True)

## 1 - calculation of temporal precision

# 1 - copy_files_for_data (copy all relevant files into 'data' directory)
def copy_files_for_data(state_dict, custom_source_dir1 = None, custom_source_dir2 = None, custom_source_dir3 = None):
    """
    Copies three types of files important to generate trajectory models:
    -.config files created with start_sim.py in Snapshot_Modeling (source_dir1)
    -time-dependent stoichiometry data for each timepoint. Data should be presented in .csv file. With this function all
    csv file in source_dir2 will be copied. These .csv files will be used in the exp_comp dictionary in create_DAG
    function
    -scoresA and scoresB for each snapshot created with imp sampcon exhaust
    (source_dir1 + snapshot + good_scoring_models) are merged into total score .txt using merge_scores helper function.
    All copied files are gathered in newly created './data/' directory, where everything is prepared for create_DAG
    function.


    :param state_dict (dict): state_dict: dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_source_dir1 (optional - str): Custom path to heterogeneity modeling dir (heterogeneity_modeling.py),
    to copy .config files
    :param custom_source_dir2 (optional - str): Custom path to stoichiometry data dir
    :param custom_source_dir3 (optional - str): Custom path to snapshot modeling dir (start_sim.py), to copy .config
    files and to access scoresA/scoresB (custom_source_dir3 + snapshot{state}_{time} + 'good_scoring_models')
    """
    # Create the destination directory for all the data copied in this function
    destination_dir = './data/'
    os.makedirs(destination_dir, exist_ok=True)

    # path to snapshot modeling dir
    if custom_source_dir1:
        source_dir1 = custom_source_dir1
    else:
        source_dir1 = '../../Heterogeneity/Heterogeneity_Modeling/'

    # path to stoichiometry data dir
    if custom_source_dir2:
        source_dir2 = custom_source_dir1
    else:
        source_dir2 = '../../Input_Information/gen_FCS/'

    # path to snapshot modeling dir
    if custom_source_dir3:
        source_dir3 = custom_source_dir3
    else:
        source_dir3 = '../../Snapshots/Snapshots_Modeling/'

    # Copy all .config files from the first source directory to the destination directory
    try:
        for file_name in os.listdir(source_dir1):
            if file_name.endswith('.config'):
                full_file_name = os.path.join(source_dir1, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, destination_dir)
        print(".config files are copied")
    except Exception as e:
        print(f".config files cannot be copied. Try do do it manually. Reason for Error: {e}")

    # Copy all .csv stoichiometry files from the second source directory to the destination directory
    try:
        for file_name in os.listdir(source_dir2):
            if file_name.endswith('.csv'):
                full_file_name = os.path.join(source_dir2, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, destination_dir)
        print(".csv stoichiometry files are copied")
    except Exception as e:
        print(f".csv stoichiometry files cannot be copied. Try do do it manually. Reason for Error: {e}")

    # Copy scoresA and scoresB from the snapshot_{state}_{time} directories and first source directory path
    try:
        for time in state_dict.keys():
            for state in range(1, state_dict[time] + 1):
                snapshot_dir = os.path.join(source_dir3, f'snapshot{state}_{time}')
                good_scoring_models_dir = os.path.join(snapshot_dir, 'good_scoring_models')
                if os.path.isdir(good_scoring_models_dir):
                    for score_file in ['scoresA.txt', 'scoresB.txt']:
                        full_file_name = os.path.join(good_scoring_models_dir, score_file)
                        if os.path.isfile(full_file_name):
                            new_file_name = f'{state}_{time}_{os.path.splitext(score_file)[0]}.log'
                            shutil.copy(full_file_name, os.path.join(destination_dir, new_file_name))
                            print(f"Copied {full_file_name} to {os.path.join(destination_dir, new_file_name)}")
    except Exception as e:
        print(f"scoresA.txt and scoresB.txt cannot be copied. Try do do it manually. Reason for Error: {e}")

os.chdir(main_dir)
# copy all the relevant files
copy_files_for_data(state_dict, custom_source_dir1='../modeling/Heterogeneity/Heterogeneity_Modeling/',
                   custom_source_dir2='../modeling/Input_Information/gen_FCS/',
                   custom_source_dir3='../modeling/Snapshots/Snapshots_Modeling/')

# create two independent DAGs
expected_subcomplexes = ['A', 'B', 'C']
exp_comp = {'A': 'exp_compA.csv', 'B': 'exp_compB.csv', 'C': 'exp_compC.csv'}
input = "./data/"
outputA = "../output_modelA/"
outputB = "../output_modelB/"

# Output from sampling precision and model precision to be saved in united dir: analysis_output_precision
analysis_output = "./analysis_output_precision/"
os.makedirs(analysis_output, exist_ok=True)

nodesA, graphA, graph_probA, graph_scoresA = spatiotemporal.create_DAG(state_dict, out_pdf=True, npaths=3,
                                                                        input_dir=input, scorestr='_scoresA.log',
                                                                        output_dir=outputA,
                                                                        spatio_temporal_rule=True,
                                                                        expected_subcomplexes=expected_subcomplexes,
                                                                        score_comp=True, exp_comp_map=exp_comp,
                                                                        draw_dag=False)

os.chdir(main_dir)
nodesB, graphB, graph_probB, graph_scoresB = spatiotemporal.create_DAG(state_dict, out_pdf=True, npaths=3,
                                                                        input_dir=input, scorestr='_scoresB.log',
                                                                        output_dir=outputB,
                                                                        spatio_temporal_rule=True,
                                                                        expected_subcomplexes=expected_subcomplexes,
                                                                        score_comp=True, exp_comp_map=exp_comp,
                                                                        draw_dag=False)

## 1 - analysis
analysis.temporal_precision(outputA + 'labeled_pdf.txt', outputB + 'labeled_pdf.txt',
                            output_fn='.' + analysis_output + 'temporal_precision.txt')
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
print("Step 1: calculation of temporal precision IS COMPLETED")
print("")
print("")

## 2 - calculation of precision of the model

# precision is calculated from .labeled_pdf.txt in Trajectories_Modeling dir
trajectories_modeling_input_dir = "./output/"

analysis.precision(trajectories_modeling_input_dir + 'labeled_pdf.txt', output_fn=analysis_output + 'precision.txt')

os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
print("Step 2: calculation of precision of the model IS COMPLETED")
print("")
print("")

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
        _labeled_pdf = custom_labeled_pdf
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

# 3b - comparison of the model to data used in modeling (copy number)
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
forward_model_copy_number(expected_subcomplexes,custom_labeled_pdf='output/labeled_pdf.txt')
print("Step 3b: copy number validation IS COMPLETED")
print("")
print("")

# 4a - SAXS
"""
Comparing center models of the most dominant cluster for each snapshot (rmfs) to the SAXS data for each time point
 can be done in two steps:
-converting rmfs to pdb files
-comparing pdbs of each snapshot to experimental SAXS profile using FoXS
"""

def convert_rmfs(state_dict, model, custom_path=None):
    """
    The purpose of this function is to automate the conversion of RMF files into PDB files for all the states from
    state_dict. Created PDBs are further used in comparison of SAXS profiles using FoXS. Additionally, they can be
    used for comparison to native PDB if available.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param model (str): An IMP (Integrative Modeling Platform) model object.
    :param custom_path (optional - str): A custom path for the RMF file, allowing for flexibility in file location
    (should be compliant with stat_dict).
    """

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            if custom_path:
                sim_rmf = custom_path # option for custom path
            else:
                sim_rmf = f"../../modeling/Snapshots/Snapshots_Assessment/exhaust_{state}_{time}/cluster.0/cluster_center_model.rmf3"

            pdb_output = f"snapshot{state}_{time}.pdb" # define the output of converted .pdb file

            if os.path.exists(sim_rmf):
                try:
                    rmf_fh = RMF.open_rmf_file_read_only(sim_rmf) # open rmf file for reading
                    rmf_hierarchy = IMP.rmf.create_hierarchies(rmf_fh, model)[0] # extract 1st hierarchy
                    IMP.atom.write_pdb_of_c_alphas(rmf_hierarchy, pdb_output) # write coordinates of CA to .pdb
                    print(f"Finishing: snapshot{state}_{time}.pdb")
                except Exception as e:
                    print(f"{sim_rmf} is empty or there is another problem: {e}")


def copy_SAXS_dat_files(custom_src_dir = None):
    """
    Copies all files ending with .dat from the specified directory to the current directory.

    :param custom_src_dir (optional - str): Path to the source directory
    """
    if custom_src_dir:
        src_dir = custom_src_dir
    else:
        src_dir = '../../../Input_Information/gen_SAXS'
    try:
        files = os.listdir(src_dir) # Get the list of all files in the src_dir directory
        dat_files = [f for f in files if f.endswith('.dat')] # Filter out files that end with .dat

        # Copy each .dat file to the current directory, so FoXS can be used
        for file_name in dat_files:
            full_file_name = os.path.join(src_dir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
                # print(f"Copied: {full_file_name} to {main_dir}")

        print("All .dat files have been copied successfully...")

    except Exception as e:
        print(f"An error occurred: {e}")


def process_foxs(state_dict, custom_dat_file = None):
    """
    This function automates the FoXS analysis for all specified time points in a single execution. It processes PDB
    files generated by the convert_rmfs function and uses SAXS data copied with the copy_SAXS function. All of this
    data should be present in the current running directory.
    FoXS tutorial is available here: https://integrativemodeling.org/tutorials/foxs/foxs.html

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_dat_file (optional - str)): A custom name of SAXS files for each time point (should be
    compliant with stat_dict)
    """


    print("...lets proceed to FoXS")

    for time in state_dict.keys():
        try:
            if state_dict[time] > 1:
                # if there is more than one state in timepoint, FoXS creates fit.plt and it should be renamed
                if custom_dat_file:
                    dat_file = custom_dat_file
                else:
                    dat_file = f"{time}_exp.dat"

                pdb_files = " ".join([f"snapshot{state}_{time}.pdb" for state in range(1, state_dict[time] + 1)])

                command1 = f"foxs -r -g {pdb_files} {dat_file}"
                # example how FoXS command should look like: foxs -r -g snapshot1_0min.pdb snapshot2_0min.pdb snapshot3_0min.pdb 0min_exp.dat
                os.system(command1)
                print(f"FoXS for {time} is calculated and ready to create a plot. Nr of states is: {state_dict[time]}")

                command2 = f"gnuplot fit.plt" # create plot from .plt code
                os.system(command2)

                command3 = f"mv fit.plt {time}_FoXS.plt" # rename .plt to avoid to be overwritten
                os.system(command3)

                command4 = f"mv fit.png {time}_FoXS.png" # rename plot to avoid to be overwritten
                os.system(command4)

                print(f"Plot {time}_FoXS.png is created")

            elif state_dict[time] == 1:
                print(f"There is only one state in {time}")
                dat_file1 = f"{time}_exp.dat"
                pdb_file1 = f"snapshot1_{time}.pdb"

                command5 = f"foxs -r -g {pdb_file1} {dat_file1}"
                os.system(command5)
                print(f"FoXS for {time} is calculated and ready to create a plot. Nr of states is: {state_dict[time]}")

                command6 = f"gnuplot snapshot1_{time}_{time}_exp.plt"
                os.system(command6)

                command7 = f"mv snapshot1_{time}_{time}_exp.plt {time}_FoXS.plt"
                os.system(command7)

                command8 = f"mv snapshot1_{time}_{time}_exp.png {time}_FoXS.png"
                os.system(command8)
            else:
                print(f"There is no states in this timepoint. Check stat_dict.")

        except Exception as e:
            print(f"FoXS can not be executed properly due to following problem: {e}")


# 4a - SAXS
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
SAXS_output = "./SAXS_comparison/"
os.makedirs(SAXS_output, exist_ok=True)
os.chdir(SAXS_output)
model = IMP.Model()
convert_rmfs(state_dict, model)
copy_SAXS_dat_files(custom_src_dir='../../modeling/Input_Information/gen_SAXS')
process_foxs(state_dict)
print("Step 4a: SAXS validation IS COMPLETED")
print("")
print("")
