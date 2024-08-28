Trajectory modeling step 5: assessment {#trajectory_assess}
====================================

Now that the set of spatiotemporal models has been constructed, we must evaluate these models. We can evaluate these models in at least 4 ways: estimate the sampling precision, compare the model to data used to construct it, validate the model against data not used to construct it, and quantify the precision of the model.

Navigate to `Trajectories/Trajectories_Assessment` and run `trajectories_assessment.py`. This code will perform the following steps to assess the model.

# Sampling precision

To begin, we calculate the sampling precision of the models. The sampling precision is calculated by using `spatiotemporal.create_DAG` to reconstruct the set of trajectory models using 2 independent sets of samplings for snapshot models. Then, the overlap between these snapshot models is evaluated using `analysis.temporal_precision`, which takes in two `labeled_pdf` files.

\code{.py}
# state_dict - universal parameter
state_dict = {'0min': 3, '1min': 3, '2min': 1}
# current directory
main_dir = os.getcwd()
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
\endcode

The output of `analysis.temporal_precision` is written in `analysis_output_precision/temporal_precision.txt`, shown below. The temporal precision can take values between 1.0 and 0.0, and indicates the overlap between the two models in trajectory space. Hence, values close to 1.0 indicate a high sampling precision, while values close to 0.0 indicate a low sampling precision. Here, the value close to 1.0 indicates that sampling does not effect the weights of the trajectory models.

\code{.txt}
Temporal precision between ../outputA/labeled_pdf.txt and ../outputB/labeled_pdf.txt:
1.0
\endcode

# Model precision

Next, we calculate the precision of the model, using `analysis.precision`. Here, the model precision calculates the number of trajectories with high weights. The precision ranges from 1.0 to 1/d, where d is the number of trajectories. Values approaching 1.0 indicate the model set can be described by a single trajectory, while values close to 1/d indicate that all trajectories have similar weights.

\code{.py}
## 2 - calculation of precision of the model

# precision is calculated from .labeled_pdf.txt in Trajectories_Modeling dir
trajectories_modeling_input_dir = "../Trajectories_Modeling/output/"

analysis.precision(trajectories_modeling_input_dir + 'labeled_pdf.txt', output_fn=analysis_output + 'precision.txt')

os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
print("Step 2: calculation of precision of the model IS COMPLETED")
print("")
print("")
\endcode

The `analysis.precision` function reads in the `labeled_pdf` of the complete model, and writes the output file to `analysis_output_precision/precision.txt`, shown below. The value close to 1.0 indicates that the set of models can be sufficiently represented by a single trajectory.

\code{.txt}
Precision of ../Trajectories_Modeling/output/labeled_pdf.txt:
1.0
\endcode

# Comparison against data used in model construction

We then evaluate the model against data used in model construction. First, we calculate the cross-correlation between the original EM map and the forward density projected of each snapshot model. We wrote the `ccEM` function to perform this comparison for all snapshots.

\code{.py}
# 3a - comparison of the model to data used in modeling (EM)
exp_mrc_base_path = "../../Input_Information/ET_data/add_noise"
ccEM(exp_mrc_base_path)
print("Step 3a: ET validation IS COMPLETED")
print("")
print("")
\endcode

The output of `ccEM` is written in `ccEM_output/`. It contains forward densities for each snapshot model (`MRC_{state}_{time}.mrc`) and `ccEM_calculations.txt`, which contains the cross-correlation to the experimental EM profile for each snapshot.

After comparing the model to EM data, we aimed to compare the model to copy number data, and wrote the `forward_model_copy_number` function to evaluate the copy numbers from our set of trajectory models.

\code{.py}
# 3b - comparison of the model to data used in modeling (copy number)
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
forward_model_copy_number(expected_subcomplexes)
print("Step 3b: copy number validation IS COMPLETED")
print("")
print("")
\endcode

The output of `forward_model_copy_number` is written in `forward_model_copy_number/`. The folder contains `CN_prot_{prot}.txt` files for each protein, which have the mean and standard deviation of protein copy number at each time point.

Here, we plot the comparison between the experimental data used in model construction and the set of trajectory models. This analysis includes comparisons between experimental and modeled protein copy numbers for A (a), B (b), and C (c), as well as the cross-correlation coefficient between the experimental EM density and the forward density of the set of sufficiently good scoring modeled structures in the highest weighted trajectory (d). Here, we see the model is in good agreement with the data used to construct it.

\image html Spatiotemporal_Assessment_Used.png width=1200px

# Validation against data not used in model construction

\code{.py}
## 4 - comparison of the model to data used in modeling (SAXS, native pdb of final complex)
# 4a - SAXS
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
SAXS_output = "./SAXS_comparison/"
os.makedirs(SAXS_output, exist_ok=True)
os.chdir(SAXS_output)
convert_rmfs(state_dict, model)
copy_SAXS_dat_files()
process_foxs(state_dict)
print("Step 4a: SAXS validation IS COMPLETED")
print("")
print("")
\endcode

\code{.py}
# 4b - RMSD
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
pdb_path = "../../Input_Information/PDB/3rpg.pdb"
RMSD(pdb_path=pdb_path, custom_n_plot=20)
print("Step 4b: RMSD validation IS COMPLETED")
print("")
print("")
\endcode

\image html Spatiotemporal_Assessment_Unused.png width=1200px

\image html Chi2_Table.png width=600px
