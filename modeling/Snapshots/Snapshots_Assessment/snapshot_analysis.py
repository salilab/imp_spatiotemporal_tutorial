"""
Analyzing snapshots consists of series of different functions:
- extracting_stat_files: extracts and collects all the relevant scores from all the runs
- general_rule_calculation: median value calculation of EM cross correlation for each snapshot
(general rule for filtering good scoring snapshots)
- general_rule_filter_independent_samples: filtering frames based on general rule (median value of EM cross correlation)
 and creating two independent samples for each snapshot from filtered frames
- create_histogram: create different "score histograms" to visualize score distributions after filtering
- count_rows_and_generate_report: counts number of rmfs in each of two independent samples to check if distribution is
similar
- !! create_density_dictionary (still needs to be written): based on .config file, it creates density dictionaries used
in exhaust function
- exhaust (most important part of this code): based on all filtered frames for each snapshot it checks if sampling of
snapshots is exhaustive (output plots) and generate centroid models of the most dominant cluster for each snapshot
- extract_exhaust_data: from exhaust output for each snapshot KS values and sampling precision (two most important
parameters) are gathered in one .txt file
- save_exhaust_data_as_png: convert output .txt file generated with extract_exhaust_data function into .png table for
better visualization
"""

import pandas as pd
import os
import ast
import re
import matplotlib.pyplot as plt


def extract_values_from_file(file_path, keys):
    """
    This function is a helper function which extract from stat file (stat.X.out) all the required values associated
    with keys provided in the 'keys_to_extract' parameter.

    :param file_path (str): path to each stat.X.out from which data should be extracted
    :param keys (list of int): 'keys_to_extract' parameter from the extracting_stat_files function.
    Look at the header of stat.X.out file to find which dictionary key corresponds to which variable.
    :return (dict): for each stat.X.out dictionary is created (keys and corresponding list of values). For our example
    there is only one key that we are interested in (key 3 - 'GaussianEMRestraint_None_CCC')
    """

    # Open the file and reads all the frames into "content"
    with open(file_path, 'r') as file:
        content = file.readlines()

    # Initializes an empty list 'data' to store dictionaries obtained from the file.
    data = []
    # To skip first dictionary
    skip_first_dict = True
    dict_lines = []

    # The purpose of this part of the code is to skip the first dictionary in the file
    for line in content:
        if skip_first_dict:
            dict_lines.append(line.strip())
            # Check if we've reached the end of the first dictionary
            if dict_lines.count('{') == dict_lines.count('}'):
                # When we reach the end of first dictionary, skipping is turned off
                skip_first_dict = False
            continue

        # Tries to parse the stripped line into a dictionary using ast.literal_eval
        # If successful, appends the dictionary to data. If it fails, it continues to the next line.
        try:
            dictionary = ast.literal_eval(line.strip())
            data.append(dictionary)
        except:
            continue

    # New dictionary where the extracted data from each dictionary in the 'data' list should be stored
    extracted_data = {key: [] for key in keys}

    # Extracting values from the 'data' list and combining them with corresponding keys
    for dictionary in data:
        for key in keys:
            extracted_data[key].append(dictionary.get(key, ''))

    return extracted_data


def extracting_stat_files(state_dict, runs_nr, replica_nr, replica_output_name, keys_to_extract, decimals_nr, custom_base_path = None):
    """
    This function extracts and collects all the relevant scores from all the runs. To understand this function, it is
    important to understand stat files (stat.X.out) created with pmi snapshot.py and IMP.pmi.macros.ReplicaExchange
    command in snapshot.py.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param runs_nr (int): number of runs generated with start_sim.py for each snapshot{state}_{time}
    :param replica_nr (int): number of processor required for ReplicaExchange set in Job.sh
    :param replica_output_name (str): name of output file set in IMP.pmi.macros.ReplicaExchange command
    in pmi snapshot.py
    :param keys_to_extract (list): List of strings which represents keys in stat file dictionary. For example: number 3
    corresponds to 'GaussianEMRestraint_None_CCC' in first dictionary of our model.
    Therefore, in each further dictionary in stat.X.out values related to key '3' corresponds to
    GaussianEMRestraint_None_CCC for certain frame.
    :param decimals_nr (int): Number of decimals that should be extracted from dictionaries of stat.X.out.
    For example: GaussianEMRestraint_None_CCC values have 16 decimals
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    :return: {state}_{time}_stat.txt is created in each snapshot{state}_{time} directory
    """
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"snapshot{state}_{time}"
            dir_path = os.path.join(base_path, dir_name)

            if not os.path.exists(dir_path):
                print(f"Directory {dir_path} does not exist. Skipping.")
                continue

            # Create the output file once for each directory
            output_file_path = os.path.join(dir_path, f"{state}_{time}_stat.txt")
            with open(output_file_path, 'w') as out_file:
                out_file.write("\t".join(map(str, keys_to_extract)) + "\n")

                for run in range(1, runs_nr + 1):
                    run_dir = os.path.join(dir_path, f"run{run}")

                    if not os.path.exists(run_dir):
                        print(f"Run directory {run_dir} does not exist. Skipping.")
                        continue

                    combined_data = {key: [] for key in keys_to_extract}

                    for stat_file_num in range(replica_nr):
                        stat_file_path = os.path.join(run_dir, replica_output_name, f"stat.{stat_file_num}.out")

                        if not os.path.exists(stat_file_path):
                            print(f"Stat file {stat_file_path} does not exist. Skipping.")
                            continue

                        extracted_data = extract_values_from_file(stat_file_path, keys_to_extract)

                        for key in keys_to_extract:
                            combined_data[key].extend(extracted_data[key])

                    # Transpose and write the data
                    rows = zip(*[combined_data[key] for key in keys_to_extract])
                    for row in rows:
                        out_file.write(
                            "\t".join(map(lambda x: f"{x:.{decimals_nr}f}" if isinstance(x, float) else str(x),
                                          row)) + "\n")

            print(f"Stat file {state}_{time}_stat.txt created.")


def general_rule_calculation(state_dict, general_rule_column, custom_general_rule_file = None, custom_base_path = None):
    """
    From each {state}_{time}_stat.txt created with extracting_stat_files function median value of desired column is
    calculated and results are gathered in general_rule_file.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param general_rule_column (int): Column from {state}_{time}_stat.txt on which general rule should be applied.
    For example: column '3' corresponds to 'GaussianEMRestraint_None_CCC' values in our example.
    :param custom_general_rule_file (optional - str): Custom name of output file created with this function
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    :return: general_rule_file .txt saved in the directory where this code is

    NOTE: if desired, name 'Median_Value_ccEM' can be manually changed in the 'columns' list
    """

    # optional parameters
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"

    if custom_general_rule_file:
        general_rule_file = custom_general_rule_file
    else:
        general_rule_file = 'general_rule.txt'

    output_data = [] # Here all the data is temporarily saved
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"snapshot{state}_{time}" # dir name is the first column in the general_rule_file .txt
            file_path = os.path.join(base_path, dir_name, f'{state}_{time}_stat.txt')

            # From {state}_{time}_stat.txt different following things are extracted
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, delim_whitespace=True)

                num_rows = df.shape[0]  # Number of rows, meaning number of generated rmfs for certain snapshot
                median_value = df[f'{general_rule_column}'].median() # calculation of general rule
                output_data.append([dir_name, num_rows, median_value])
            print(f'Extracting data for {dir_name}')

    output_df = pd.DataFrame(output_data,
                             columns=['Directory', 'Number_of_Generated_Frames', 'Median_Value_ccEM']) # header
    output_df.to_csv(general_rule_file, index=False, sep='\t') # create an output file: general_rule_file



def general_rule_filter_independent_samples(state_dict, main_dir, custom_general_rule_file = None, custom_base_path = None):
    """
    This function has two main roles: (i) filtering frames based on general rule and (ii) creating two independent
    samples for each snapshot from filtered frames. Both steps can be achieved by submitting the imp_sampcon select_good
    command in each snapshot{state}_{time} directory. More about imp_sampcon select_good can be found here:
    https://integrativemodeling.org/tutorials/actin/analysis.html
    (i) median value calculated in general_rule_file .txt is set as a lower filtering threshold. Filtered frames (rmfs)
    are then saved in  each snapshot{state}_{time} directory, where new directory 'good_scoring_models' is created.
    (ii) The flag -e in imp_sampcon select_good command results in creating 'good_scoring_models', a folder with the
    models that passed through our filtering criteria. Data from this directory is later used for exahust function as
    well as in the 'Trajectories_Assessment' step.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param main_dir (str): directory where this code is
    :param custom_general_rule_file (optional - str): Custom name of output file created with this function
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    """
    # optional parameters
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"

    if custom_general_rule_file:
        general_rule_file = custom_general_rule_file
    else:
        general_rule_file = 'general_rule.txt'

    # Read the text file and create temporary dictionary directory (keys) vs median values (values)
    median_values = {}
    with open(general_rule_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file: # in each line we are only interested in the 1st and 3rd element
            parts = line.split()  # Split the line into a list of strings
            directory = parts[0]  # First part is the directory name
            general_rule = parts[2]  # Third part is the median value
            median_values[directory] = general_rule  # Store in dictionary

    # Execute command in each directory with the corresponding median value
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"snapshot{state}_{time}"
            if dir_name in median_values:
                median_value = median_values[dir_name]
                command = f"imp_sampcon select_good -rd . -rp run \
       -sl 'GaussianEMRestraint_None_CCC'\
       -pl ConnectivityRestraint_Score \ExcludedVolumeSphere_Score \GaussianEMRestraint_None \
       Total_Score -alt {median_value} -aut 1.0 -e"

                print(f"Executing in {dir_name}: {command}")
                os.chdir(os.path.join(base_path, dir_name))  # Navigate to each snapshot directory where the command (-rd) is executed
                os.system(command)  # Execute the command
                os.chdir(main_dir)  # Change back to the universal directory to restart the loop


def create_histograms(state_dict, main_dir, score_list, custom_base_path = None):
    """
    This function creates different "score histograms" to visualize score distributions after filtering based on
    general rule. Histograms are created in each newly created histograms{state}_{time} directory by submitting
    imp_sampcon plot_score command there. More about imp_sampcon plot_score can be found here:
    https://integrativemodeling.org/tutorials/actin/analysis.html

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param main_dir (str): directory where this code is
    :param score_list (list of str): columns in 'model_ids_scores.txt' for which histograms should be created.
    For example: for this tutorial we created histograms for all scores defined in the imp_sampcon select_good
    (inside the general_rule_filter_independent_samples function)
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    """

    # optional parameter
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            # For each snapshot separate directory should be created to avoid overlapping of the plots, also this makes eveything more clear
            dir_name = f"histograms{state}_{time}"
            # Construct the full path to plots{state}_{time} directory where command is executed
            full_path = os.path.join(main_dir, dir_name)
            os.makedirs(full_path, exist_ok=True)
            os.chdir(full_path)
            for score in score_list:
                # command is executed in newly created histogram{state}_{time} directory, therefore additional '../' should be added
                command = f"imp_sampcon plot_score  ../{base_path}/snapshot{state}_{time}/good_scoring_models/model_ids_scores.txt \
      {score}"
                os.system(command) # Execute the command
                print(f"Executing command: {command}")
            os.chdir(main_dir) # Change back to the universal directory to restart the loop after all the histograms are created

def count_rows_and_generate_report(state_dict, custom_output_file_path = None, custom_base_path = None):
    """
    This function counts number of rmfs in each of two independent samples to check if distribution is similar. Such
    brief check can predict if we can expect some problems (especially with K-S test) in achieving sampling
    exhaustiveness in the next step (exhaust function). Writes a file to output_file_path.txt with the number of
    samples that passed filtering in each independent set of sampling runs.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_output_file_path: Custom name of output_file_path .txt
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    """

    # optional parameters
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"

    if custom_output_file_path:
        output_file_path = custom_output_file_path
    else:
        output_file_path = 'independent_samples_stat.txt'

    output_data = []
    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            snapshot_dir = f'snapshot{state}_{time}'
            good_scoring_models_dir = os.path.join(base_path, snapshot_dir, 'good_scoring_models')

            # Initialize row counts
            rows_A = 0
            rows_B = 0

            # Count rows in scoresA.txt
            scoresA_path = os.path.join(good_scoring_models_dir, 'scoresA.txt')
            if os.path.exists(scoresA_path):
                with open(scoresA_path, 'r') as file:
                    rows_A = sum(1 for line in file)

            # Count rows in scoresB.txt
            scoresB_path = os.path.join(good_scoring_models_dir, 'scoresB.txt')
            if os.path.exists(scoresB_path):
                with open(scoresB_path, 'r') as file:
                    rows_B = sum(1 for line in file)

            # Calculate total rows
            total_rows = rows_A + rows_B

            # Append data to the output list
            output_data.append([snapshot_dir, rows_A, rows_B, total_rows])

    # Write the output data to a text file
    with open(output_file_path, 'w') as output_file:
        # Write header
        output_file.write('Directory\tRows in scoresA.txt\tRows in scoresB.txt\tTotal Rows\n')

        # Write data rows
        for data in output_data:
            output_file.write('\t'.join(map(str, data)) + '\n')

    print(f"Output written to {output_file_path}") # print to check that function is completed

# !! create_density_dictionary function (still needs to be written)

def exhaust(state_dict, main_dir, custom_base_path = None):
    """
    Based on all filtered frames for each snapshot this function runs 'imp_sampcon exhaust' command in each newly
    created exhaust_{state}_{time} directory. This command checks if independently sampled runs come from the same
    distribution, uses clustering to determine the precision at which sampling is exhaustive,
    generates centroid models of the most dominant cluster for each snapshot, and writes these outputs out to files.
    More about imp_sampcon exhaust can be found here:
    https://integrativemodeling.org/tutorials/actin/analysis.html

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param main_dir (str): directory where this code is
    :param custom_base_path (optional - str): Custom path to the directory where snapshot{state}_{time} created with
    start_sim.py are
    """

    # optional parameter
    if custom_base_path:
        base_path = custom_base_path
    else:
        base_path = "../Snapshots_Modeling"

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            # Create directory
            dir_name = f"exhaust_{state}_{time}"
            full_path = os.path.join(main_dir, dir_name)
            os.makedirs(full_path, exist_ok=True)

            # command is executed in newly created exhaust_{state}_{time} directory,
            # therefore additional '../' should be added
            command = f"imp_sampcon exhaust \
       -n snapshot_{state}_{time} -p ../{base_path}/snapshot{state}_{time}/good_scoring_models/ -d ../{state}_{time}_density_ranges.txt \
       -m cpu_omp -c 8 -a -g 1.0 -gp"

            print(f"Executing in {dir_name}: {command}")
            os.chdir(full_path)
            os.system(command)
            os.chdir(main_dir)



def extract_exhaust_data(state_dict, custom_KS_sampling_precision_output = None):
    """
    This function gathers from exhaust output for each snapshot KS values and sampling precision into one .txt file.
    This .txt file (together with output plots) can serve as an overview for determining
    sampling exhaustiveness.

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param custom_KS_sampling_precision_output (optional -str): Custom name of output KS_sampling_precision_output .txt
    """

    if custom_KS_sampling_precision_output:
        KS_sampling_precision_output = custom_KS_sampling_precision_output
    else:
        KS_sampling_precision_output = 'KS_sampling_precision_output.txt'

    output_lines = ["snapshot KS_D_value KS_p_value sampling_precision"]

    for time in state_dict.keys():
        for state in range(1, state_dict[time] + 1):
            dir_name = f"./exhaust_{state}_{time}"
            snapshot = f'{state}_{time}'
            # from each exhaust_{state}_{time} created with exhaust function KS data and sampling precision is extracted
            KS_test_file_path = os.path.join(dir_name, f"snapshot_{state}_{time}.KS_Test.txt")
            sampling_precision_file_path = os.path.join(dir_name,
                                                        f"snapshot_{state}_{time}.Sampling_Precision_Stats.txt")

            if os.path.exists(KS_test_file_path) and os.path.exists(sampling_precision_file_path):
                with open(KS_test_file_path, 'r') as KS_file:
                    line = KS_file.readline().strip()
                    if line:
                        D_value, p_value = line.split()

                # best way to extract sampling precision values is directly from the text
                with open(sampling_precision_file_path, 'r') as sp_file:
                    content = sp_file.read()
                    # extracting precision values is directly from the text (.* is different name for each snapshot)
                    match = re.search(r"The sampling precision for our .* modeling is ([\d\.]+)", content)
                    if match:
                        sampling_precision = match.group(1)
                        print(f"Sampling precision for snapshot{snapshot} is: {sampling_precision}")

                output_lines.append(f"{snapshot} {D_value} {p_value} {sampling_precision}")

    with open(KS_sampling_precision_output, 'w') as output_file:
        output_file.write("\n".join(output_lines))


def save_exhaust_data_as_png(custom_input_file = None, custom_output_png = None):
    """
    This function converts output .txt file generated with extract_exhaust_data function into .png table for
    better visualization. Created .png is saved in the running directory next to KS_sampling_precision_output .txt.
    :param custom_input_file:
    :param custom_output_png:
    """
    # Custom parameters
    if custom_input_file:
        input_file = custom_input_file
    else:
        input_file = "KS_sampling_precision_output.txt"

    if custom_output_png:
        output_png = custom_output_png
    else:
        output_png = "KS_sampling_precision_output.png"

    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, delim_whitespace=True)

    # Format the D_value column to have 5 decimal places
    df['KS_D_value'] = df['KS_D_value'].apply(lambda x: f"{float(x):.5f}")

    # Format the sampling_precision column to have 3 decimal places
    df['sampling_precision'] = df['sampling_precision'].apply(lambda x: f"{float(x):.3f}")

    # Create the table using matplotlib and saving it as a PNG image
    fig, ax = plt.subplots(figsize=(8, len(df) * 0.4))  # Adjust the size based on the number of rows
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Save the table as a PNG image
    plt.savefig(output_png, bbox_inches='tight', dpi=300)
    print(f".png is created")



# Example usage
if __name__ == "__main__":
    # state_dict - universal parameter
    state_dict = {'0min': 3, '1min': 3, '2min': 1}
    # current directory
    main_dir = os.getcwd()

    # 1 calling extracting_stat_files function and related parameters
    keys_to_extract = [3]
    runs_nr = 50
    replica_nr = 16
    replica_output_name = 'output'
    decimals_nr = 16

    extracting_stat_files(state_dict, runs_nr, replica_nr, replica_output_name, keys_to_extract, decimals_nr)
    print("extracting_stat_files is DONE")
    print("")
    print("")

    # 2 calling general_rule_calculation and related parameters
    general_rule_column = '3'

    general_rule_calculation(state_dict, general_rule_column)

    print("general_rule_calculation is DONE")
    print("")
    print("")


    # 3 calling general_rule_filter_independent_samples
    general_rule_filter_independent_samples(state_dict, main_dir)
    print("general_rule_filter_independent_samples is DONE")
    print("")
    print("")

    # 4 calling create_histograms and related parameters
    score_list = [
        'Total_Score',
        'ConnectivityRestraint_Score',
        'ExcludedVolumeSphere_Score',
        'GaussianEMRestraint_None',
        'GaussianEMRestraint_None_CCC'
    ] # list of histograms we want to create in each histograms{state}_{time} directory

    create_histograms(state_dict, main_dir, score_list)
    print("create_histograms is DONE")
    print("")
    print("")

    # 5 calling count_rows_and_generate_report
    count_rows_and_generate_report(state_dict)
    print("calling count_rows_and_generate_report is DONE")
    print("")
    print("")

    # 6 !! create_density_dictionary (still needs to be written)

    # 7 calling exhaust
    exhaust(state_dict, main_dir)
    print("exhaust is DONE")
    print("")
    print("")

    # 8 calling extract_exhaust_data
    extract_exhaust_data(state_dict)
    print("extract_exhaust_data is DONE")
    print("")
    print("")

    # 9 calling save_exhaust_data_as_png
    save_exhaust_data_as_png()
    print("save_exhaust_data_as_png is DONE")
    print("")
    print("")









