'''
This code has two main functions:
-prepare_protein_library (IMP.spatiotemporal.prepare_protein_library.prepare_protein_library): based on
time - dependent stoichiometry data this function calculates configurations of possible snapshot models for
each time together with corresponding pmi topology file.
For example: for our system for first two time points we calculated three most possible configurations
(check state_dict - important dictionary that describe our snapshot model).
For more information regarding prepare_protein_library check: (add link!)

-generate_all_snapshots: (after configuration and topology files are created) this function automatically creates pmi
 model representation script for each of the snapshots by incorporating corresponding EM data and topology file
'''
import IMP
import numpy as np
import itertools
import sys
import os
import math
import pandas as pd
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import prepare_protein_library
import os
import shutil


def generate_all_snapshots(state_dict, main_dir, items_to_copy, job_template, number_of_runs):
    '''
    NOTE: Before running generate_all_snapshots function, main_dir should contain:
    -ET_data: This directory should be manually created based on desired time-dependent EM restraint for pmi representation.
    For example: For our system, from ../../Input_Information/ET_data only 'noisy' experimental data (.gmm and .rmc)
    was copied into ./experimental directory. This ./experimental subdirectory is created as .gmm densities generated
    can be stored here for each snapshot{state}_{time}.
    -{state}_{time}_topol.txt topology files created with IMP.spatiotemporal.prepare_protein_library.prepare_protein_library
    (for this function ./gen_FCS directory should be copied from Input_Information directory and 'spatiotemporal_topology.txt'
    should be written)
    -static_snapshot.py
    For more information regarding static_snapshot.py and spatiotemporal_topology.txt check this tutorial:
    https://integrativemodeling.org/tutorials/actin/pmidesign.html

    The purpose of this function is that it can automatically create pmi and run model representation script for each
    of the snapshots by incorporating corresponding EM data and topology file at once. This can be achieved with following
    steps:
    a) Copy all three directories/files to each snapshot{state}_{time} directory and Job.sh is generated
    b) Create run directories and copy files in each "child" run directory inside "parent" snapshot{state}_{time} directory
    c) The last step is to submit Job.sh using qsub inside each run directory <state

    :param state_dict (dict): dictionary that defines the spatiotemporal model.
           The keys are strings that correspond to each time point in the
           stepwise temporal process. Keys should be ordered according to the
           steps in the spatiotemporal process. The values are integers that
           correspond to the number of possible states at that timepoint.
    :param main_dir (str): directory where this function is
    :param items_to_copy (list): List to copy in each snapshot{state}_{time} directory (except {state}_{time}_topol.txt)
    :param job_template (str - text): Job.sh that allows to run desired number of runs for all the snapshots across
    all timepoints simultaneously. The purpose of Job.sh is to submit static_snapshot.py for each snapshot{state}_{time}
    with running parameters <state> <time>. In this way snapshots (static_snapshot.py) are generated with corresponding
    EM restraint (ET_data) and topology file ({state}_{time}_topol.txt).
    :param number_of_runs (int): Desired number of runs for each snapshot{state}_{time}.
    :return: qsub Job.sh - to run desired number of runs for all the snapshots across all timepoints (snapshot{state}_{time})
    '''

    for time in state_dict.keys():
        for state in range(1,state_dict[time]+1):
            new_dir = f"snapshot{state}_{time}" # Name of each new directory
            # IMPORTANT: name snapshot{state}_{time} for each snapshot directory is used in all further functions
            os.makedirs(new_dir, exist_ok=True)

            # a) Copy directories and files to each snapshot{state}_{time} directory
            for item in items_to_copy:
                item_path = os.path.join(main_dir, item)
                if os.path.isdir(item_path):
                    shutil.copytree(item_path, os.path.join(new_dir, item))
                else:
                    shutil.copy2(item_path, new_dir)

            topol_file = f"{state}_{time}_topol.txt"
            shutil.copy2(os.path.join(main_dir, topol_file), new_dir)

            # Create Job.sh
            job_script_path = os.path.join(new_dir, "Job.sh")
            with open(job_script_path, "w") as job_script:
                job_script.write(job_template.format(state=state, time=time))
            # After this, in each snapshot{state}_{time} directory we have:
            # -ET_data directory
            # -newly created Job.sh
            # -static_snapshot.py code
            # -specific topology file

            # Create run directories
            # b) Copy files in each "child" run directory inside "parent" snapshot{state}_{time} directory
            for run in range(1, number_of_runs+1):
                run_dir = os.path.join(new_dir, f"run{run}")
                os.makedirs(run_dir, exist_ok=True)
                shutil.copy2(job_script_path, run_dir)
                shutil.copy2(os.path.join(new_dir, 'static_snapshot.py'), run_dir)
                # After this, in each run directory we have:
                # -newly created Job.sh
                # -static_snapshot.py code

                # c) The last step is to submit Job.sh using qsub inside each run directory
                os.system(f"cd {run_dir} && qsub Job.sh")
                print(f"All Job.sh for snapshot{state}_{time} have been successfully submitted")
                # Check then qstat what is running


if __name__ == "__main__":
    # state_dict - universal parameter
    state_dict = {'0min': 3, '1min': 3, '2min': 1}

    # 1a - parameters for prepare_protein_library:
    times = ["0min", "1min", "2min"]
    exp_comp = {'A': './gen_FCS/exp_compA.csv', 'B': './gen_FCS/exp_compB.csv',
                'C': './gen_FCS/exp_compC.csv'}
    expected_subcomplexes = ['A', 'B', 'C']
    template_topology = 'spatiotemporal_topology.txt'
    template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}
    nmodels = 3

    # 1b - calling prepare_protein_library
    prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                    template_topology=template_topology, template_dict=template_dict)


    # 2a - parameters for generate_all_snapshots

    main_dir = os.getcwd()
    items_to_copy = ['ET_data', 'static_snapshot.py']  # additionally we need to copy only specific topology file
    job_template = """#!/bin/bash
    #$ -S /bin/bash
    #$ -cwd
    #$ -r n
    #$ -j y
    #$ -N test1
    #$ -pe smp 16
    #$ -l h_rt=48:00:00

    module load Sali
    module load mpi/openmpi-x86_64
    module load imp
    module load python3/scikit/0.21.3

    mpirun -np $NSLOTS python3 static_snapshot.py {state} {time}
    """
    number_of_runs = 50

    # 2b - calling generate_all_snapshots
    generate_all_snapshots(state_dict, main_dir, items_to_copy, job_template, number_of_runs)



