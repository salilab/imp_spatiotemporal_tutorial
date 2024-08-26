Snapshot modeling step 5: assessment {#snapshot_assess}
====================================

Now, we have a variety of alternative snapshot models. In general, we would like to assess these models in at least 4 ways: estimate the sampling precision, compare the model to data used to construct it, validate the model against data not used to construct it, and quantify the precision of the model. Here, we will focus specifically on estimating the sampling precision of the model, while quantitative comparisons between the model and experimental data will be reserved for the final step, when we assess [trajectories](https://integrativemodeling.org/tutorials/spatiotemporal/trajectory_assess.html). To assess these snapshot models, we navigate to the `Snapshots/Snapshots_Assessment` folder and run `snapshot_assessment.py`. This script performs the following analysis.

# Filtering good scoring models

Initially, we want to filter the various alternative models to select those that meet certain parameter thresholds. In this case, we filter the structural models in each snapshot by the median cross correlation with EM data. We note that this filtering criteria is subjective, and developing a Bayesian method to objectively way different restraints for filtering remains an interesting future development in integrative model.

The current filtering procedure involves three steps. In the first step, we look through the `stat.*.out` files to write out the cross correlation with EM data for each model, which, in this case, is labeled column `3`, `GaussianEMRestraint_None_CCC`. In other applications, the column that corresponds to each type of experimental data may change, depending on the scoring terms for each model. For each snapshot, a new file is written with this data (`{state}_{time}_stat.txt`).

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

In the third step, we use the `imp_sampcon select_good` tool to filter each snapshot, according to the median value determined in the previous step. For each snapshot, this function produces a file, `good_scoring_models/model_ids_scores.txt`, which contains the run, replicaID, scores, and sampleID for each model that passes filtering. It also saves RMF files with each model from two independent groups of sampling runs from each snapshot to `good_scoring_models/sample_A` and `good_scoring_models/sample_B`, writes the scores for the two independent groups of sampling runs to `good_scoring_models/scoresA.txt` and `good_scoring_models/scoresB.txt`, and writes `good_scoring_models/model_sample_ids.txt` to connect each model to its division of sampling run. More information on `imp_sampcon` is available in the analysis portion of the [actin tutorial](https://integrativemodeling.org/tutorials/actin/analysis.html).

\code{.py}
# 3 calling general_rule_filter_independent_samples
general_rule_filter_independent_samples(state_dict, main_dir)
print("general_rule_filter_independent_samples is DONE")
print("")
print("")
\endcode

# Plotting data, clustering models, and determining sampling precision

Next, scores can be plotted for analysis. Here, the `create_histograms` function runs `imp_sampcon plot_score` to plot distributions for various scores of interest. Each of these plots are saved to `histograms{state}_{time}/{score}.png`, where score is an object listed in the `score_list`. These plots are useful for debugging the modeling protocol, and should appear roughly Gaussian.

\code{.py}
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
\endcode

We then check the number of models in each sampling run by running `count_rows_and_generate_report`, which writes the `independent_samples_stat.txt` file. Empirically, we have found checking the overall number of models in each independent sample that pass filtering serves a good first check on sampling convergence.

\code{.py}
# 5 calling count_rows_and_generate_report
count_rows_and_generate_report(state_dict)
print("count_rows_and_generate_report is DONE")
print("")
print("")
\endcode

Next, we calculate the density range dictionaries, which are output as `{state}_{time}_density_ranges.txt`. These dictionaries label each protein in each snapshot model, which will be passed into `imp_sampcon` to calculate the localization density of each protein.

\code{.py}
# 6 calling create_density_dictionary:
create_density_dictionary_files(state_dict, main_dir)
print("create_density_dictionary is DONE")
print("")
print("")
\endcode

Finally, we run `imp_sampcon exhaust` on each snapshot. This code performs checks on the exhaustiveness of the sampling. Specifically it analyzes the convergence of the model score, whether the two model sets were drawn from the same distribution, and whether each structural cluster includes models from each sample proportionally to its size. The output for each snapshot is written out to the `exhaust_{state}_{time}` folder.

\code{.py}
# 7 calling exhaust
exhaust(state_dict, main_dir)
print("exhaust is DONE")
print("")
print("")
\endcode

Plots for determining the sampling precision are shown below for a single snapshot, 1_2min. (a) Tests the convergence convergence of the lowest scoring model (`snapshot_{state}_{time}.Top_Score_Conv.pdf`). Error bars represent standard deviations of the best scores, estimated by selecting different subsets of models 10 times. The light-blue line indicates a lower bound reference on the total score. (b) Tests that the scores of two independently sampled models come from the same distribution (`snapshot_{state}_{time}.Score_Dist.pdf`). The difference between the two distributions, as measured by the KS test statistic (D) and KS test p-value (p) indicates that the difference is both statistically insignificant (p>0.05) and small in magnitude (D<0.3). (c) Determines the structural precision of a snapshot model (`snapshot_{state}_{time}.ChiSquare.pdf`). RMSD clustering is performed at 1 Å intervals until the clustered population (% clustered) is greater than 80%, and either the χ<sup>2</sup> p-value is greater than 0.05 or Cramer’s V is less than 0.1. The sampling precision is indicated by the dashed black line. (d) Populations from sample 1 and sample 2 are shown for each cluster (`snapshot_{state}_{time}.Cluster_Population.pdf`).

\image html Snapshot_Exhaust.png width=1200px

Further structural analysis can be calculated by using the `cluster.*` files. The `cluster.*.{sample}.txt` files contain the model number for the models in that cluster, where `{sample}` indicates which round of sampling the models came from. The `cluster.*` folder contains an RMF for centroid model of that cluster, along with the localization densities for each protein. The localization densities of each protein from each independent sampling can be compared to ensure independent samplings produce the same results.

Ideally, each of these plots should be checked for each snapshot. As a way to summarize the output of these checks, we can gather the results of the KS test and the sampling precision test for all snapshots. This is done by running `extract_exhaust_data` and `save_exhaust_data_as_png`, which write `KS_sampling_precision_output.txt` and `KS_sampling_precision_output.png`, respectively.

\code{.py}
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
\endcode

These codes write a table that include the KS two sample test statistic (D), the KS test p-value, and the sampling precision for each snapshot model, which is replotted below.

\image html Snapshot_sampling.png width=600px

# Visualizing models

The resulting RMF files and localization densities from this analysis can be viewed in [UCSF Chimera](https://www.rbvi.ucsf.edu/chimera/) (version>=1.13) or [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/).

Here, we plotted each centroid model (A - blue, B - orange, and C - purple) from clustering and compared that model to the experimental EM profile (gray).

\image html static_snapshots_noCC.png width=600px

Finally, now that snapshot models have been assessed, we can perform [Trajectory modeling steps 2-4: representation, scoring, and search process.] (@ref trajectory1)

