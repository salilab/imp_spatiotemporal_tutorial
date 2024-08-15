"""
This code has two main functions:
-prepare_protein_library (IMP.spatiotemporal.prepare_protein_library.prepare_protein_library): based on
time - dependent stoichiometry data this function calculates configurations of possible snapshot models for
each time together with corresponding pmi topology file.
For example: for our system for first two time points we calculated three most possible configurations
(check state_dict - important dictionary that describe our snapshot model).
For more information regarding prepare_protein_library check:
https://integrativemodeling.org/nightly/doc/ref/namespaceIMP_1_1spatiotemporal_1_1prepare__protein__library.html

-generate_all_snapshots: (after configuration and topology files are created) this function automatically creates pmi
 model representation script for each of the snapshots by incorporating corresponding EM data and topology file.
 The current job script is currently configured for our cluster and will need to be edited for other systems.
"""
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
    """
    NOTE: Before running generate_all_snapshots function, main_dir should contain:
    -{state}_{time}_topol.txt topology files created with
    IMP.spatiotemporal.prepare_protein_library.prepare_protein_library
    -static_snapshot.py
    For more information regarding static_snapshot.py and spatiotemporal_topology.txt check this tutorial:
    https://integrativemodeling.org/tutorials/actin/pmidesign.html

    The purpose of this function is that it can automatically create pmi and run model representation script for each
    of the snapshots by incorporating corresponding EM data and topology file at once.
    This can be achieved with following steps:
    a) Copy all necessary files to each snapshot{state}_{time} directory and Job.sh is generated
    b) Create run directories and copy files in each "child" run directory inside "parent"
    snapshot{state}_{time} directory
    c) The last step is to submit Job.sh using qsub inside each run directory <state>

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
    EM restraint and topology file ({state}_{time}_topol.txt).
    :param number_of_runs (int): Desired number of runs for each snapshot{state}_{time}.
    """

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
            # Make directory for forward densities
            os.mkdir(os.path.join(new_dir, "forward_densities"))

            # find topology file
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
    # 1a - parameters for prepare_protein_library:
    times = ["0min", "1min", "2min"]
    exp_comp = {'A': '../../Input_Information/gen_FCS/exp_compA.csv',
                'B': '../../Input_Information/gen_FCS/exp_compB.csv',
                'C': '../../Input_Information/gen_FCS/exp_compC.csv'}
    expected_subcomplexes = ['A', 'B', 'C']
    template_topology = 'spatiotemporal_topology.txt'
    template_dict = {'A': ['Ubi-E2-D3'], 'B': ['BMI-1'], 'C': ['E3-ubi-RING2']}
    nmodels = 3

    # 1b - calling prepare_protein_library
    prepare_protein_library.prepare_protein_library(times, exp_comp, expected_subcomplexes, nmodels,
                                                    template_topology=template_topology, template_dict=template_dict)


    # 2a - parameters for generate_all_snapshots
    # state_dict - universal parameter
    state_dict = {'0min': 3, '1min': 3, '2min': 1}

    main_dir = os.getcwd()
    items_to_copy = ['static_snapshot.py']  # additionally we need to copy only specific topology file
    # jobs script will likely depend on the user's cluster / configuration
    job_template = ("#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -N Tutorial\n#$ -pe smp 16\n"
                    "#$ -l h_rt=48:00:00\n\nmodule load Sali\nmodule load imp\nmodule load mpi/openmpi-x86_64\n\n"
                    "mpirun -np $NSLOTS python3 static_snapshot.py {state} {time}")
    number_of_runs = 50

    # 2b - calling generate_all_snapshots
    generate_all_snapshots(state_dict, main_dir, items_to_copy, job_template, number_of_runs)



