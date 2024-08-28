Trajectory modeling steps 2-4: representation, scoring, and search process {#trajectory1}
====================================

Now, we have snapshot models for various intermediate states along our process of interest. In this step, we will combine these snapshot models to produce trajectories, as described in full [here](https://www.biorxiv.org/content/10.1101/2024.08.06.606842v1.abstract).

# Background behind integrative spatiotemporal modeling

## Representing the model

We choose to represent dynamic processes as a trajectory of snapshot models, with one snapshot at each time point. In this case, we modeled snapshots at 3 time points (0, 1, and 2 minutes), so a single trajectory model will consist of 3 snapshots, one at each 0, 1, and 2 minutes. The modeling procedure described here will produce a set of scored trajectory models, which can be displayed as a directed acyclic graph, where nodes in the graph represent the snapshot model and edges represent connections between snapshots at neighboring time points.

## Scoring the model

To score trajectories, we incorporate both the scores of individual snapshots, as well as the scores of transitions between them. Under the assumption that the process is Markovian (*i.e.* memoryless), the weight of a trajectory takes the form:

\f[
W(\Chi) \propto   \displaystyle\prod^{T}_{t=0} P( X_{N,t}, N_{t} | D_{t}) \cdot\\ \displaystyle\prod^{T-1}_{t=0} W(X_{N,t+1},N_{t+1} | X_{N,t},N_{t}, D_{t,t+1}),
\f]

where \f$t\f$ indexes times from 0 until the final modeled snapshot (\f$T\f$); \f$P(X_{N,t}, N_{t} | D_{t})\f$ is the snapshot model score; and \f$W(X_{N,t+1},N_{t+1} | X_{N,t},N_{t}, D_{t,t+1})\f$ is the transition score. Transition scores are currently based on a simple metric that either allows or disallows a transition. Transitions are only allowed if all proteins in the first snapshot model are included in the second snapshot model. In the future, we hope to include more detailed transition scoring terms, which may take into account experimental information or physical models of macromolecular dynamics.

## Searching for good scoring models

Trajectories are constructed by enumerating all connections between adjacent snapshots and scoring these trajectories according to the equation above. This procedure results in a set of weighted trajectories.

# Code for integrative spatiotemporal modeling

Navigate to `Trajectories/Trajectories_Modeling`. The `create_trajectories.py` script runs the above steps to create a spatiotemporal integrative model from the previously sampled snapshots.

First, we copy all relevant files to a new directory, `data`. These files are (i) `{state}_{time}.config` files, which include the subcomplexes that are in each state, (ii) `{state}_{time}_scores.log`, which is a list of all scores of all models in that snapshot, and (iii) `exp_comp{prot}.csv`, which is the copy number for each protein (`{prot}`) as a function of time.

\code{.py}
# copy all the relevant files for create_DAG
# it is important that everything starts from main dir
main_dir = os.getcwd()
os.chdir(main_dir)
state_dict = {'0min': 3, '1min': 3, '2min': 1}
create_data_and_copy_files(state_dict)
\endcode

We then build the spatiotemporal graph by running `spatiotemporal.create_DAG`, [documented here](https://integrativemodeling.org/nightly/doc/ref/namespaceIMP_1_1spatiotemporal_1_1create__DAG.html). This function represents, scores, and searches for trajectory models.

\code{.py}
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
\endcode

The inputs we included are:
- state_dict (dict): a dictionary that defines the spatiotemporal model. Keys are strings for each time point in the spatiotemporal process and values are integers corresponding to the number of snapshots modeled at that time point
- out_pdf (bool): whether to write the probability distribution function (pdf).
- npaths (int): Number of states two write to a file (path*.txt).
- input_dir (str): directory with the input information.
- scorestr (str): final characters at the end of the score files.
- output_dir (str): directory to which model will be written. Will be created if it does not exist.
- spatio_temporal_rule (bool): whether to include our transition scoring term, which enforces that all proteins in the first snapshot model are included in the second snapshot model.
- expected_subcomplexes (list): list of string objects, which is the subcomplexes to look when enforcing the spatiotemporal rule. Strings should be substrings of those in `{state}_{time}.config` files.
- score_comp (bool): whether to score the composition of each snapshot.
- exp_comp_map (dictionary): key is a string with the name of each protein that will undergo composition scoring, value is the `.csv` file with the copy number data for that protein.
- draw_dag (bool): whether to write out an image of the directed acyclic graph.

After running `spatiotemporal.create_DAG`, a variety of outputs are written:
- `cdf.txt`: the cumulative distribution function for the set of trajectories.
- `pdf.txt`: the probability distribution function for the set of trajectories.
- `labeled_pdf.txt`: Each row has 2 columns and represents a different trajectory model. The first column labels a single trajectory as a series of snapshots, where each snapshot is written as `{state}_{time}|` in sequential order. The second column is the probability distribution function corresponding to that trajectory.
- `dag_heatmap.eps` and `dag_heatmap`: image of the directed acyclic graph from the set of models.
- `path*.txt`: files where each row includes a `{state}_{time}` string, so that rows correspond to the states visited over that trajectory model. Files are numbered from most likely to least likely.

Now that we have a trajectory model, we can plot the directed acyclic graph (left) and the series of centroid models from each snapshot along the most likely trajectory (right). Each row corresponds to a different time point in the assembly process (0 min, 1 min, and 2 min). Each node is shaded according to its weight in the final model (\f$W(X_{N,t}N_{t})\f$). Proteins are colored as A - blue, B - orange, and C - purple.

\image html Spatiotemporal_Model.png width=600px

Finally, now that the set of spatiotemporal models has been constructed, we can [assess the spatiotemporal integrative models.] (@ref trajectory_assess)

