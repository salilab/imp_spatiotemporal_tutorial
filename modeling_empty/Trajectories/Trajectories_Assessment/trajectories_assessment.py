"""
Script for the assessment of final trajectory model. This script assess model in four ways:
-calculation of temporal precision
-calculation of precision of the model
-comparison of the model to data used in modeling (EM)
-comparison of the model to data not used in modeling (SAXS, native pdb of final complex)
"""
import IMP
import IMP.atom
import RMF
import IMP.rmf
import IMP.algebra
import IMP.core
import IMP.pmi.macros
import os
import shutil
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import numpy as np
import matplotlib.pyplot as plt


## 1 - calculation of temporal precision
"""
To calculate temporal precision, two independent trajectory models (A and B) should be generated from two separate 
samples (SampleA and SampleB) obtained from each snapshot. This is done most simply by referencing two 
independent sets of scores generated by sampcon_exhaust. The similarity of these independent sampled scores for 
each snapshot can be easily verified in the KS-test.txt file.

Temporal precision is calculated by using the sequence of three functions:
-copy_files_for_data (copy all relevant files into 'data' directory)
-create_DAG (from IMP.spatiotemporal)
-analysis (from IMP.spatiotemporal)
"""
# 1 - copy_files_for_data (copy all relevant files into 'data' directory)
def copy_files_for_data(state_dict, custom_source_dir1 = None, custom_source_dir2 = None):
    """
    Copies three types of files important to generate two independent trajectory models (A and B):
    -.config files created with start_sim.py in Snapshot_Modeling (source_dir1)
    -time-dependent stoichiometry data for each timepoint. Data should be presented in .csv file. With this function all
    csv file in source_dir2 will be copied. These .csv files will be used in the exp_comp dictionary in create_DAG
    function
    -scoresA and scoresB for each snapshot created with imp sampcon exhaust
    (source_dir1 + snapshot + good_scoring_models)

    :param state_dict (dict): state_dict: dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_source_dir1 (optional - str): Custom path to snapshot modeling dir (start_sim.py), to copy .config
    files and scoresA/scoresB (custom_source_dir1 + snapshot{state}_{time} + 'good_scoring_models')
    :param custom_source_dir2 (optional - str): Custom path to stoichiometry data dir
    """
    # Create the destination directory for all the data copied in this function
    destination_dir = './data/'
    os.makedirs(destination_dir, exist_ok=True)

    # path to snapshot modeling dir
    if custom_source_dir1:
        source_dir1 = custom_source_dir1
    else:
        source_dir1 = '../../Snapshots/Snapshots_Modeling/'

    # path to stoichiometry data dir
    if custom_source_dir2:
        source_dir2 = custom_source_dir1
    else:
        source_dir2 = '../../Input_Information/gen_FCS/'

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
                snapshot_dir = os.path.join(source_dir1, f'snapshot{state}_{time}')
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


##
def extract_coordinates(rmf_hierarchy,rmf_fh):
    """
    This function extracts coordinates from the rmf hierarchy.

    :param rmf_hierarchy (str): hierarchy of certain rmf file from where coordinates are extracted
    :param rmf_fh (str): rmf file opened in read only
    :return (list): list of extracted coordinates
    """
    coords = []
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(0)) # in each rmf there is only one frame
    for molecule in rmf_hierarchy.get_children()[0].get_children():
        for leaf in IMP.core.get_leaves(molecule):
            coords.append(leaf)

    return coords


# for each snapshot separately it creates mrc
def write_mrc(scoring_path, mrc_file, MRCresolution=10.0,voxel=5.0):
    """
    This function accesses rmfs from both samples (sample_A and sample_B) for each snapshot{state}_{time}, opens them in
    read only and extracts hierarchy. After that, helper function extract_coordinates is called to convert hierarchy to
    coordinates. Used by ccEM function.
    Based on coordinates, density maps are calculated (using IMP.em.SampledDensityMap) and added in .mrc file
    (IMP.em.write_map).

    :param scoring_path (str): path to rmfs
    :param mrc_file (str): path to created mrc file
    :param MRCresolution (float): set MRC resolution
    :param voxel (float): set voxel for .mrc
    :return dmap2: (IMP DensityMap class) Cumulative density map for certain snapshot{state}_{time}
    """
    count_frames = 0
    for sample in ['sample_A', 'sample_B']:
        sample_path = os.path.join(scoring_path, sample)
        if os.path.exists(sample_path):
            for file in os.listdir(sample_path):
                if file.endswith('.rmf3'):
                    sim_rmf = os.path.join(sample_path, file)

                    rmf_fh = RMF.open_rmf_file_read_only(sim_rmf)
                    model = IMP.Model()
                    rmf_hierarchy = IMP.rmf.create_hierarchies(rmf_fh, model)[0]

                    ps = extract_coordinates(rmf_hierarchy, rmf_fh)

                    # calculate density map
                    dmap = IMP.em.SampledDensityMap(ps, MRCresolution, voxel)
                    dmap.calcRMS()
                    dmap.set_was_used(True)
                    # dmap2 stores the overall density map for the snapshot{state}_{time}
                    # dmap is the density map for current rmf
                    # dmap3 is a temporary variable for combining the current dmap with dmap2

                    if count_frames == 0:
                        dmap2 = dmap

                    else:
                        bbox1 = IMP.em.get_bounding_box(dmap2)
                        bbox2 = IMP.em.get_bounding_box(dmap)
                        bbox1 += bbox2
                        dmap3 = IMP.em.create_density_map(bbox1, voxel)
                        dmap3.set_was_used(True)
                        dmap3.add(dmap)
                        dmap3.add(dmap2)
                        dmap2 = dmap3
                    count_frames = count_frames + 1

                if count_frames == 0:
                    print('Warning: RMF empty, no frames were read')
                else:
                    # density normalization and write overall density map
                    dmap2.multiply(1. / count_frames)
                    IMP.em.write_map(dmap2, mrc_file, IMP.em.MRCReaderWriter())
        else:
            print(f"{sample_path} is empty, invalid or path is incorrect")
    return dmap2

## 3 - comparison of the model to data used in modeling (EM)

def ccEM(exp_mrc_base_path, custom_output_directory = None, custom_base_path = None):
    """
    This function constructs paths to experimental .mrc and to 'good_scoring_models' directory for each
    snapshot{state}_{time}, where rmfs are saved in two independent samples (sample_A and sample_B). Density map (mrc)
    for certain snapshot{state}_{time} is calculated by calling write_mrc function. After that, cross-correlation
    between experimental density map and calculated density map of each snapshot is calculated. All calculations are
    gathered in combined .txt file, which is saved together with all the calculated .mrc files in output_directory.
    Calculated .mrc files and experimental .mrc files can be also visually compared in one "overlapping" ChimeraX
    session.

    :param exp_mrc_base_path (str): path to directory with time-dependent EM data
    :param custom_output_directory (optional - str): If desired, different name of output directory
    (where plots and .txt files are saved) can be set. Default name: ccEM_output
    :param custom_base_path (optional -str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are.
    """

    # create output directory
    if custom_output_directory:
        output_directory = custom_output_directory
    else:
        output_directory = "./ccEM_output"

    os.makedirs(output_directory, exist_ok=True)

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            # Construct directory path
            if custom_base_path:
                base_path = custom_base_path
            else:
                base_path = "../../Snapshots/Snapshots_Modeling"

            snapshot = f"snapshot{state}_{time}"  # directory snapshot{state}_{time} is created with start_sim.py
            print(f"Now we are extracting from {snapshot}")

            # Creating .mrc file for each snapshot, where average density map should be saved
            sim_mrc_file = f"MRC_{state}_{time}.mrc"
            sim_mrc_path = os.path.join(output_directory, sim_mrc_file)

            # construct path to the experimental .mrc files for each time point
            exp_mrc_file = f'{time}_noisy.mrc' # name can be changed, but it is important that time variable is included
            exp_mrc_path = os.path.join(exp_mrc_base_path, exp_mrc_file)

            # Continue path to rmf3 files
            directory_path = os.path.join(base_path, snapshot)
            # rmfs from sample_A and sample_B should be extracted
            good_scoring_path = os.path.join(directory_path, 'good_scoring_models')

            # preparing .mrc for CC calculation
            # calling a function to write a .mrc file for each snapshot
            model_density = write_mrc(good_scoring_path, sim_mrc_path)
            exp_density = IMP.em.read_map(exp_mrc_path, IMP.em.MRCReaderWriter())

            # calculation of cross-correlation between experimental density map and
            # calculated density map of each snapshot
            cc = IMP.em.get_coarse_cc_coefficient(exp_density, model_density, 0, True)

            # Create a .txt file to save all ccEM calculations together with mrc files
            txt_file = 'ccEM_calculations.txt'
            txt_file_path = os.path.join(output_directory, txt_file)
            header = 'snapshot, cc EM value\n'

            # Check if the file already exists
            try:
                with open(txt_file_path, 'x') as file:
                    file.write(header)
            except FileExistsError:
                pass  # File already exists, so we don't need to write the header again

            with open(txt_file_path, 'a') as file:
                file.write(f'{snapshot}, {cc}\n')
                print(f"Data for {snapshot} is saved: {txt_file_path}")


## 4 - comparison of the model to data used in modeling (SAXS, native pdb of final complex)
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
                sim_rmf = f"../snapshots/snapshot_analysis/exhaust_{state}_{time}/cluster.0/cluster_center_model.rmf3"

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
        src_dir = '../../Input_Information/gen_SAXS'
    try:
        files = os.listdir(src_dir) # Get the list of all files in the src_dir directory
        dat_files = [f for f in files if f.endswith('.dat')] # Filter out files that end with .dat

        # Copy each .dat file to the current directory, so FoXS can be used
        for file_name in dat_files:
            full_file_name = os.path.join(src_dir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
                print(f"Copied: {full_file_name} to {main_dir}")

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

                command7 = f"mv snapshot1_{time}_{time}_exp.plt {time}_FoXS.png"
                os.system(command7)
            else:
                print(f"There is no states in this timepoint. Check stat_dict.")

        except Exception as e:
            print(f"FoXS can not be executed properly due to following problem: {e}")


## 4b - RMSD calculation between filtered rmfs and native pdb of final complex

"""
This series of function calculates RMSD between (i) residues coordinates of filtered rmfs from both samples 
(sampleA and sampleB) created with imp sampcom exhaust and (ii) coordinates of CA from native pdb of final complex. 
NOTE: it is important for user to understand hierarchy of rmf, so code can be changed based on the specific example
"""

def RMSD_from_dict(dict1, dict2):
    """
    This function calculates rmsd between rmf coordinate dictionary and pdb coordinate dictionary. It is important that
    both dictionaries include same number of key:value pairs, otherwise dictionaries are not aligned.

    :param dict1 (dict): rmf coordinate dictionary created with append_rmf_coordinates
    :param dict2 (dict): pdb coordinate dictionary created with create_pdb_coordinates function
    :return (float): rmsd result
    """
    Rmf = list(dict1.keys())
    Pdb = list(dict2.keys())
    if len(Rmf) != len(Pdb):
        print("Cannot calculate RMSD properly as residues are not aligned. Check hierarchy in pdb, rmf...")
        return None
    coords1 = []
    coords2 = []
    for i in range(0, len(Rmf)):
        coords1.append(dict1[Rmf[i]])
        coords2.append(dict2[Pdb[i]])

    from IMP.algebra import get_rmsd
    rmsd = get_rmsd(coords1, coords2)
    return rmsd


def create_pdb_coordinates(coord_list, pdb_hierarchy, chains_to_include):
    """
    This function extract name:coordinate dictionary from native pdb of final complex.
    This dictionary is later used in the RMSD_from_dict function, which calculates RMSD.
    NOTE: This function is specifically designed to access the coordinates of CA for the different proteins named
    with "Chain X" in the pdb used in this tutorial. For different systems, alternative "hierarchical loops"  may be
    required to access the desired names or coordinates from the pdb hierarchy. Use muted (#) 'print statements'
    to explore and understand the hierarchy for your specific system.

    :param coord_list (list): empty list set in RMSD function where pdb dictionary should be saved
    :param pdb_hierarchy (str): native pdb hierarchy of final complex from where coordinates are extracted
    :param chains_to_include (list): list of names from the pdb to include in calculation (based on the .config file
    for each snapshot{state}_{time}). Default naming here is "Chain X".
    :return:
    """
    coords = {}
    for pdb_molecule in pdb_hierarchy.get_children():
        pdb_molecule_names = pdb_molecule.get_name()
        # use this to understand naming of proteins in pdb file
        # print(f"Names of molecules in pdb file are: {pdb_molecule_names}")
        # it should select only the chains from configuration file created
        if pdb_molecule_names in chains_to_include:
            for pdb_leaf in IMP.core.get_leaves(pdb_molecule):
                pdb_leaf_name = pdb_leaf.get_name()
                if pdb_leaf_name.startswith(f"Atom CA of residue"):
                    leaf_coords = IMP.core.XYZ(pdb_leaf).get_coordinates()
                    pdb_combined_name = f"{pdb_molecule_names} leaf: {pdb_leaf_name}"
                    # use this to understand naming and hierarchy of pdb
                    # print(f"Understanding hierarchy. Combined name is: {pdb_combined_name}")
                    coords[pdb_combined_name] = leaf_coords
                    # use this to understand how coordinates are extracted from pdb
                    # print(f"Coords pdb list: {coords}")
    coord_list.append(coords)
    return coord_list


def append_rmf_coordinates(coord_list, rmf_hierarchy, rmf_fh):
    """
    This function extract name:coordinate dictionary from rmf hierarchy for certain snapshot{state}_{time}.
    This dictionary is later used in the RMSD_from_dict function, which calculates RMSD.
    NOTE: This function is specifically designed to work with the rmf hierarchy for the provided tutorial example.
    For different systems, alternative "hierarchical loops"  may be required to access the desired coordinates from the
    hierarchy. Use muted (#) 'print statements' to explore and understand the hierarchy for your specific system.

    :param coord_list (list): list of list of dictionaries, where each item in the list is a dictionary for a
    structural state, the keys in that dictionary are the names of proteins and the values are the coordinates for the
    corresponding protein.
    :param rmf_hierarchy (str): hierarchy of certain rmf file from where coordinates are extracted
    :param rmf_fh (str): rmf file opened in read only
    :return:  coord_list (list): Updated coordinate list, with coordinates for the loaded state.
    """
    coords = {}
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(0))  # there is only one frame in each rmf3
    for molecule in rmf_hierarchy.get_children()[0].get_children():
        molecule_names = molecule.get_name()
        # use this to understand naming of proteins in rmf files
        # print(f"Names of molecules in pdb file are: {molecule_names}")
        # Naming should be same as it was set in topology file
        for leaf in IMP.core.get_leaves(molecule):
            leaf_name = leaf.get_name()
            leaf_coords = IMP.core.XYZ(leaf).get_coordinates()
            combined_name = f"{molecule_names} leaf:{leaf_name}"
            # use this to understand naming and hierarchy of rmf
            # print(f"Understanding hierarchy. Combined name is: {combined_name}")
            coords[combined_name] = leaf_coords
            # use this to understand how coordinates are extracted from rmf
            # print(f"Coords pdb list: {coords}")

    coord_list.append(coords)
    return coord_list


def RMSD(pdb_path, custom_n_plot=None, custom_output_directory=None, custom_base_path=None):
    """
    This function reads a .config file generated by the imp sampcon exhaust function
    for each state to determine which chains from the pdb  file should be included in the calculation.
    NOTE: Users must modify this section of the code if the chain names in the pdb file differ from the default "Chain X"
    for certain "sub part" of protein based on which RMSD is calculated

    Overall, this function accesses the pdb file and the filtered rmfs, then calls additional functions
    (create_pdb_coordinates and append_rmf_coordinates). This process extracts a dictionary for the pdb file and a list
    of dictionaries for all filtered rmfs. Using these dictionaries, the final RMSD is calculated by calling the RMSD_from_dict
    function and all the results for each snapshot are gathered in the rmsd_results list. Based on this list plot of RMSD values
    distribution for each snapshot is generated and also common .txt file with all the results.


    :param pdb_path (str): path to the pdb of final complex
    :param custom_n_plot (optional - int): Optional downsampling of points visualized in RMSD_results plot.
    For clear plot, it is suggested to have less than 1000 RMSD values plotted.
    :param custom_output_directory (optional - str): If desired, different name of output directory (where plots and .txt
    files are saved) can be set. Default name: "RMSD_calculation_output"
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with start_sim.py are.
    :return: ?? (maybe return should be here rmsd_results, but I am not sure)
    """
    # create output directory where all RMSD output is gathered
    if custom_output_directory:
        output_directory = custom_output_directory
    else:
        output_directory = "./RMSD_calculation_output"

    os.makedirs(output_directory, exist_ok=True)

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            # Construct directory path
            # Base path is path to the directory where all the snapshot{state}_{time} are saved and .config files (both created with start_sim.py)
            if custom_base_path:
                base_path = custom_base_path
            else:
                base_path = "../../Snapshots/Snapshots_Modeling"
            snapshot = f"snapshot{state}_{time}"  # directory snapshot{state}_{time} is created with start_sim.py
            print(f"Now we are extracting from {snapshot}")

            directory_path = os.path.join(base_path, snapshot)

            # Read the configuration file created with start_sim.py
            config_path = os.path.join(base_path,
                                       f"{state}_{time}.config")
            if not os.path.exists(config_path):
                print(f"Configuration file {config_path} not found. Check base path and output")
                continue
            with open(config_path, 'r') as file:
                chains_to_include = ["Chain " + line.strip() for line in file if
                                     line.strip()]  # "A" in .config file is equivalent to "Chain A" in pdb naming
                print(f"{snapshot} chains: {chains_to_include}")

            # for each snapshot separate RMSD should be calculated
            rmsd_results = []

            scoring_path = os.path.join(directory_path, 'good_scoring_models')
            for sample in ['sample_A',
                           'sample_B']:  # RMSD for both independent samples for each snapshot should be calculated
                sample_path = os.path.join(scoring_path, sample)
                if os.path.exists(sample_path):
                    for file in os.listdir(sample_path):
                        if file.endswith('.rmf3'):
                            sim_rmf = os.path.join(sample_path, file)  # this is path to rmf file
                            # print(sim_rmf) # use it in case you need to know if rmfs are accurately accessed

                            if os.path.exists(sim_rmf):  # continue only in the case that rmfs can be accessed
                                try:
                                    # Read in native rmf
                                    rmf_coord_list = []  # for each snapshot separately
                                    rmf_fh = RMF.open_rmf_file_read_only(sim_rmf)
                                    # Where do I need to set the model? Here or at the beginning?
                                    # model = IMP.Model()
                                    rmf_hierarchy = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
                                    rmf_coord_list = append_rmf_coordinates(rmf_coord_list, rmf_hierarchy, rmf_fh)
                                except Exception as e:
                                    print(f"{sim_rmf} is empty or invalid: {e}")

                                try:
                                    # Read pdb
                                    pdb_native_list = []
                                    # model = IMP.Model()
                                    pdb_hierarchy = IMP.atom.read_pdb(pdb_path, model)
                                    pdb_native_list = create_pdb_coordinates(pdb_native_list, pdb_hierarchy,
                                                                             chains_to_include)
                                    # print(pdb_native_list) # this works
                                except Exception as e:
                                    print(f"{pdb_path} is empty or invalid: {e}")

                                    # Call the function to calculate rmsd.
                                    # Using the loop to calculate RMSD for each rmf separately
                                for l in range(len(rmf_coord_list)):
                                    rmsd_value = RMSD_from_dict(rmf_coord_list[l], pdb_native_list[0])
                                    rmsd_results.append(rmsd_value)
                                    # use this if you want to check how the rmf coordinates are extracted
                                    # print(rmf_coord_list)

            print(f"Extracting rmsd for {snapshot} is finished ")

            # To handle exception for empty rmsd_results list - very unlikely
            try:
                number_of_frames = len(rmsd_results)
                print(f"number_of_frames {number_of_frames}")

                min_rmsd_value = min(rmsd_results)
                print(f"min_rmsd_value {min_rmsd_value}")

                average_rmsd_value = np.mean(rmsd_results)
                print(f"average_rmsd_value {average_rmsd_value}")

                # Create .txt file in the output directory

                txt_file_name = f"RMSD_analysis.txt"
                txt_file_path = os.path.join(output_directory, txt_file_name)

                # Define the header
                header = ['snapshot', 'min_value', 'number_of_frames', 'average_rmsd_value']

                try:
                    with open(txt_file_path, 'x', newline='') as file:
                        file.write('\t'.join(header) + '\n')
                except FileExistsError:
                    pass  # File already exists
                # Append the new row of data
                with open(txt_file_path, 'a', newline='') as file:
                    file.write(
                        f"{snapshot}\t{min_rmsd_value}\t{number_of_frames}\t{average_rmsd_value}\n")

                # Create plots:
                rmsd_frames = list(range(number_of_frames))

                # If there is too many points (rmfs) to make plot clear, downsampling with custom_n_plot param should be considered
                if custom_n_plot:
                    n_plot = custom_n_plot
                else:
                    n_plot = 1

                rmsd_frames_downsampled = rmsd_frames[::n_plot]
                rmsd_results_downsampled = rmsd_results[::n_plot]

                plt.figure(figsize=(10, 6))
                plt.plot(rmsd_frames_downsampled, rmsd_results_downsampled, marker='o', linestyle='-', color='b')
                plt.xlabel('Frame')
                plt.ylabel('RMSD')
                plt.title(f'RMSD for {state}_{time}')
                plt.grid(True)
                plot_filename = f"{output_directory}/rmsd_{state}_{time}.png"
                plt.savefig(plot_filename)
                plt.close()

            except ValueError as e:
                print(f"Cannot calculate statistics, rmsd_results is empty: {e}")




if __name__ == "__main__":
    # state_dict - universal parameter
    state_dict = {'0min': 3, '1min': 3, '2min': 1}
    # model
    model = IMP.Model()


    # start calling codes
    ## 1 - calculation of temporal precision
    # copy all the relevant files
    copy_files_for_data(state_dict)

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
    print("")

    ## 2 - calculation of precision of the model

    # precision is calculated from .labeled_pdf.txt in Trajectories_Modeling dir
    trajectories_modeling_input_dir = "../Trajectories_Modeling/output/"

    analysis.precision(trajectories_modeling_input_dir + 'labeled_pdf.txt', output_fn=analysis_output + 'precision.txt')

    os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
    print("Step 2: calculation of precision of the model IS COMPLETED")
    print("")
    print("")

    # 3 - comparison of the model to data used in modeling (EM)
    exp_mrc_base_path = "../../data/ET_data/experimental"
    ccEM(exp_mrc_base_path)


    ## 4 - comparison of the model to data used in modeling (SAXS, native pdb of final complex)
    # 4a - SAXS
    convert_rmfs(state_dict, model)
    copy_SAXS_dat_files()
    process_foxs(state_dict)
    print("Step 4a: SAXS validation IS COMPLETED")
    print("")
    print("")

    # 4b - RMSD
    pdb_path = "../../snapshots/PDB/3rpg.pdb"
    RMSD(pdb_path=pdb_path, custom_n_plot=20)
    print("Step 4a: SAXS validation IS COMPLETED")
    print("")
    print("")



