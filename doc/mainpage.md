Integrative spatiotemporal modeling in IMP {#mainpage}
====================================

[TOC]

# Introduction {#introduction}

Biomolecules are constantly in motion; therefore, a complete depiction of their function must include their dynamics instead of just static structures.
We have developed an integrative spatiotemporal approach to model dynamic systems.
Our approach begins by computing a series of static "snapshot models" at discrete time points, which is a set of structural models with the same composition, with each subcomplex assigned to a specific location in the fully assembled complex.
Then, we connect snapshots at neighboring time points to produce trajectories. These trajectories can be scored based on both the scores of static structures and the transitions between them, allowing for the creation of trajectories that are in agreement with the input information by construction.

If you use this tutorial or its accompanying method, please site the corresponding publications:

- Latham, A.P.; Tempkin, J.O.B.; Otsuka, S.; Zhang, W.; Ellenberg, J.; Sali, A. bioRxiv, 2024, https://doi.org/10.1101/2024.08.06.606842.
- Rožič et al. (tutorial)

# Integrative spatiotemporal modeling workflow {#steps}

In general, integrative modeling proceeds through five steps (i. gathering information, ii. choosing the model representation, iii. scoring alternative models, iv. searching for good scoring models, and v. assessing the models). In integrative spatiotemporal modeling, these five steps are repeated for each the modeling of static snapshots and the modeling of trajectories.

\image html Overview.png width=600px

This tutorial will walk you through the links below, which contain a breakdown of each of these modeling steps.

- \subpage gatherinfo

- \subpage snapshot1

- \subpage snapshot_assess

- \subpage trajectory1

- \subpage trajectory_assess

To work through this example, a variety of python packages will be necessary in addition to [IMP](https://integrativemodeling.org/). These packages are [numpy](https://numpy.org/), [os](https://docs.python.org/3/library/os.html), [warnings](https://docs.python.org/3/library/warnings.html), [sys](https://docs.python.org/3/library/sys.html), [itertools](https://docs.python.org/3/library/itertools.html), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [pyRMSD](https://pypi.org/project/pyRMSD/), and [graphviz](https://graphviz.org/). [Gnuplot](http://www.gnuplot.info/) is also used for plotting SAXS profiles.

The tutorial is available through its [GitHub link](https://github.com/salilab/imp_spatiotemporal_tutorial). A complete, worked example of the tutorial is available in the `modeling` folder; meanwhile, the necessary code and input information without the resulting model are available in the `modeling_empty` folder, from which the user can run the code for themselves.

Our tutorial begins by [gathering information for both snapshot and trajectory modeling.] (@ref gatherinfo)