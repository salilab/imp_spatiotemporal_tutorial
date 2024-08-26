Snapshot modeling step 5: assessment {#snapshot_assess}
====================================

Now, we have a variety of alternative snapshot models. In general, we would like to assess these models in at least 4 ways: estimate the sampling precision, compare the model to data used to construct it, validate the model against data not used to construct it, and quantify the precision of the model. Here, we will focus specifically on estimating the sampling precision of the model, while quantitative comparisons between the model and experimental data will be reserved for the final step, when we assess [trajectories](https://integrativemodeling.org/tutorials/spatiotemporal/trajectory_assess.html).

# Filtering good scoring models

Initially, we want to the various alternative models to select those that meet certain parameter thresholds. In this case, we filter the structural models in each snapshot by the median cross correlation with EM data. This involves three steps. In the first step, we look through the `stat.*.out` files to select the cross correlation with EM data, which, in this case, is labeled column 3, `GaussianEMRestraint_None_CCC`. In other applications, the column that corresponds to each type of experimental data may change, depending on the scoring terms for each model. For each snapshot, a new file is written with this data (`{state}_{time}_stat.txt`).

\code{.py}
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
\endcode

In the second step, we want to determine the median value of EM cross correlation for each snapshot. The `general_rule_calculation` looks through the `general_rule_column` for each `{state}_{time}_stat.txt` file and determines both the median value and the number of structures generated to determine that median value.

\code{.py}
# 2 calling general_rule_calculation and related parameters
general_rule_column = '3'

general_rule_calculation(state_dict, general_rule_column)

print("general_rule_calculation is DONE")
print("")
print("")
\endcode

In the third step, we use the `imp_sampcon select_good` tool to filter each snapshot, according to the median value determined previously. More information on `imp_sampcon` is available in the [actin tutorial](https://integrativemodeling.org/tutorials/actin/).

\code{.py}
# 3 calling general_rule_filter_independent_samples
general_rule_filter_independent_samples(state_dict, main_dir)
print("general_rule_filter_independent_samples is DONE")
print("")
print("")
\endcode

# Plotting data, clustering models, and determining sampling precision

# Visualizing models