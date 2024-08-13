"""
This function has two main parts:
-first all the required data is copied to created ./data/ repository
-then by using create_DAG trajectory model is generated.
Follow this for more examples how to use create_DAG code:
https://github.com/salilab/imp/tree/develop/modules/spatiotemporal/examples/toy
"""
import IMP.spatiotemporal as spatiotemporal
import os
import shutil

state_dict = {'0min': 3, '1min': 3, '2min': 1}

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

def create_data_and_copy_files(state_dict, custom_source_dir1 = None, custom_source_dir2 = None):
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
    :param custom_source_dir1 (optional - str): Custom path to snapshot modeling dir (start_sim.py), to copy .config
    files and to access scoresA/scoresB (custom_source_dir1 + snapshot{state}_{time} + 'good_scoring_models')
    :param custom_source_dir2 (optional - str): Custom path to stoichiometry data dir
    """

    # Create the destination directory if it does not exist (./data/). Here all the
    destination_dir = './data/'
    os.makedirs(destination_dir, exist_ok=True)

    # Path to snapshot modeling dir
    if custom_source_dir1:
        source_dir1 = custom_source_dir1
    else:
        source_dir1 = '../../Snapshots/Snapshots_Modeling/'

    # Path to stoichiometry data dir
    if custom_source_dir2:
        source_dir2 = custom_source_dir2
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
                print(f"File name to copy: {full_file_name}")
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, destination_dir)
            else:
                print(".csv file cannot be accessed")
        print(".csv stoichiometry files are copied")
    except Exception as e:
        print(f".csv stoichiometry files cannot be copied. Try do do it manually. Reason for Error: {e}")

    # Copy scoresA and scoresB from the snapshot_{state}_{time} directories and first source directory path
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"snapshot{state}_{time}"
            good_scoring_path = "good_scoring_models"
            file_a = os.path.join(source_dir1, dir_name, good_scoring_path, "scoresA.txt")
            file_b = os.path.join(source_dir1, dir_name, good_scoring_path, "scoresB.txt")
            output_file = os.path.join(destination_dir, f"{state}_{time}_scores.log") # name of the output file

            try:
                # Ensure the directory exists before try to read/write files
                if os.path.exists(file_a) and os.path.exists(file_b):
                    merge_scores(file_a, file_b, output_file) # call helper function to merge files
                    print(f"Scores for snapshot{state}_{time} have been merged and saved")
                else:  # many things can go wrong here, so it is good to know where is the problem
                    print(f"Path doesn't exist: {source_dir1}")
                    print(f"Files not found in directory: {dir_name}")
                    print(f"Files not found in directory: {file_a}")
                    print(f"Files not found in directory: {file_b}")
                    print(f"Output directory doesn't exist: {destination_dir}")
            except Exception as e:
                print(f"total scores files cannot be copied of merged. Reason for Error: {e}")



if __name__ == "__main__":
    # copy all the relevant files for create_DAG
    # it is important that everything starts from main dir
    main_dir = os.getcwd()
    os.chdir(main_dir)
    state_dict = {'0min': 3, '1min': 3, '2min': 1}
    custom_source_dir1 = '../../snapshots/snapshots'
    custom_source_dir2 = '../../snapshots/snapshots/stoichiometry_data/'

    create_data_and_copy_files(state_dict, custom_source_dir1 = custom_source_dir1, custom_source_dir2= custom_source_dir2)

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








