Integrative spatiotemporal modeling in IMP {#mainpage}
====================================

[TOC]

# Introduction {#introduction}

Biomolecules are constantly in motion; therefore, a complete depiction of their function must include their dynamics instead of just static structures.
We have developed an integrative spatiotemporal approach to model dynamic systems.
Our approach begins by computing a series of static "snapshot models" at discrete time points, which is a set of structural models with the same composition, with each subcomplex assigned to a specific location in the fully assembled complex.
Then, we connect snapshots at neighboring time points to produce trajectories. These trajectories can be scored based on both the scores of static structures and the transitions between them, allowing for the creation of trajectories that are in agreement with the input information by construction.

If you use this tutorial or its accompanying method, please site the corresponding publications:

- Latham et al. (methods)
- Rožič et al. (tutorial)

# Integrative spatiotemporal modeling workflow {#steps}

In general, integrative modeling proceeds through five steps. In integrative spatiotemporal modeling, these five steps are repeated for each the modeling of static snapshots and the modeling of trajectories.

\image html Overview.png width=600px

This tutorial will walk you through the links below, which contain a breakdown of each of these modeling steps.

- \subpage gatherinfo

- \subpage snapshot1

- \subpage snapshot_assess

- \subpage trajectory1

- \subpage trajectory_assess

First, we begin by [gathering information for both snapshot and trajectory modeling.] (@ref gatherinfo)