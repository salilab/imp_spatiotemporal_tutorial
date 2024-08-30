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

The output of `analysis.temporal_precision` is written in `analysis_output_precision/temporal_precision.txt`, shown below. The temporal precision can take values between 1.0 and 0.0, and indicates the overlap between the two models in trajectory space. Hence, values close to 1.0 indicate a high sampling precision, while values close to 0.0 indicate a low sampling precision. Here, the value close to 1.0 indicates that sampling does not affect the weights of the trajectory models.

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

We then evaluate the model against data used in model construction. First, we calculate the cross-correlation between the original EM map and the forward density projected from each snapshot model. We wrote the `ccEM` function to perform this comparison for all snapshots.

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

Here, we plot the comparison between the experimental data used in model construction and the set of trajectory models. This analysis includes the cross-correlation coefficient between the experimental EM density and the forward density of the set of sufficiently good scoring modeled structures in the highest weighted trajectory (a), as well as comparisons between experimental and modeled protein copy numbers for proteins A (b), B (c), and C (d). Here, we see the model is in good agreement with the data used to construct it.

\image html Spatiotemporal_Assessment_Used.png width=1200px

# Validation against data not used in model construction

Finally, we aim to compare the model to data not used in model construction. Specifically, we reserved SAXS data for model validation. We aimed to compare the forward scattering profile from the centroid model of each snapshot to the experimental profile. To make this comparison, we wrote functions that converted each centroid RMF to a PDB (`convert_rmfs`), copied the experimental SAXS profiles to the appropriate folder (`copy_SAXS_dat_files`), and ran [FoXS](https://integrativemodeling.org/tutorials/foxs/foxs.html) on each PDB to evaluate its agreement to the experimental profile (`process_foxs`).

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

The output of this analysis is written to `SAXS_comparison`. Standard FoXS outputs are available for each snapshot (`snapshot{state}_{time}.*`). In particular, the `.fit` files include the forward and experimental profiles side by side, with the \f$\chi^2\f$ for the fit. Further, the `{time}_FoXS.*` files include the information for all snapshots at that time point, including plots of each profile in comparison to the experimental profile (`{time}_FoXS.png`).

As our model was generated from synthetic data, the ground truth structure is known at each time point. In addition to validating the model by assessing its comparison to SAXS data, we could approximate the model accuracy by comparing the snapshot model to the PDB structure, although this comparison is not perfect as the PDB structure was used to inform the structure of *rigid bodies* in the snapshot model. To do so, we wrote a function (`RMSD`) that calculates the RMSD between each structural model and the orignal PDB.

\code{.py}
# 4b - RMSD
os.chdir(main_dir)  # it is crucial that after each step, directory is changed back to main
pdb_path = "../../Input_Information/PDB/3rpg.pdb"
RMSD(pdb_path=pdb_path, custom_n_plot=20)
print("Step 4b: RMSD validation IS COMPLETED")
print("")
print("")
\endcode

The output of this function is written in `RMSD_calculation_output`. The function outputs `rmsd_{state}_{time}.png` files, which plots the RMSD for each model within each snapshot. This data is then summarized in `RMSD_analysis.txt`, which includes the minimum RMSD, average RMSD, and number of snapshots in each snapshot model.

Here, we plot the results for assessing the spatiotemporal model with data not used to construct it. Comparisons are made between the centroid structure of the most populated cluster in each snapshot at each time point and the experimental SAXS profile for 0 (a), 1 (b), and 2 (c) minutes. Further, we plot both the sampling precision (dark red) and the RMSD to the PDB structure (light red) for each snapshot in the highest trajectory model (d).

\image html Spatiotemporal_Assessment_Unused.png width=1200px

To quantitatively compare the model to SAXS data, we used the \f$\chi^2\f$ to compare each snapshot model to the experimental profile. We note that the \f$\chi^2\f$ are substantially lower for the models along the highest trajectory (1_0min, 1_1min, and 1_2min) than for other models, indicating that the highest weighted trajectory is in better agreement with the experimental SAXS data than other possible trajectories.

\image html Chi2_Table.png width=600px

Next, we can evaluate the accuracy of the model by comparing the RMSD to the PDB to the sampling precision of each snapshot model. We note that this is generally not possible, because in most biological applications the ground truth is not known. In this case, if the average RMSD to the PDB structure is smaller than the sampling precision, the PDB structure lies within the precision of the model. We find that the RMSD is within 1.5 Å of the sampling precision at all time points, indicating that the model lies within 1.5 Å of the ground truth.

# Next steps

After assessing our model, we can must decide if the model is sufficient to answer biological questions of interest. If the model does not have sufficient precision for the desired application, assessment of the current model can be used to inform which new experiments may help improve the next iteration of the model. The [integrative spatiotemporal modeling procedure](@ref steps) can then be repeated iteratively, analogous to [integrative modeling of static structures](@ref procedure).

If the model is sufficient to provide insight into the biological process of interest, the user may decide that it is ready for publication. In this case, the user should create an [mmCIF file](https://mmcif.wwpdb.org/) to deposit the model in the [PDB-dev database](https://pdb-dev.wwpdb.org/). This procedure is explained in the [deposition tutorial](https://integrativemodeling.org/tutorials/deposition/develop/).
