Step 1: gather information for both snapshot and trajectory modeling {#gatherinfo}
====================================

# Gathering information

In the first step, information about the process of interest is gathered. This information can theoretically be applied to either modeling static snapshots or modeling trajectories, and can be used for representing the model, scoring the model, sampling alternative models, filtering sampled models, or validating the models.

For this tutorial, we used the X-ray crystal structure of the complete Bmi1/Ring1b-UbcH5c complex (a), synthetically generated electron tomography (ET) density maps during the assembly process (b), synthetically generated small-angle X-ray scattering (SAXS) profiles during the assembly process (c), and synthetically generated protein copy numbers during the assembly process (d). The crystal structure of the complex informs the final state of our model as well as the structure of the individual proteins. The time-dependent ET and SAXS data give two inputs to inform the size and shape of the assembling complex. The protein copy number data informs the stoichiometry of the complex during assembly.

\image html Input.png width=600px

In addition to the four types of data used here, a variety of data could be useful for the spatiotemporal modeling of protein complexes. IMP currently features  restraints for a variety of experimental data or prior models, including chemical cross-links, Förster resonance energy transfer, comparative structural models, and deep-learning structural models, all of which could inform spatiotemporal modeling through a procedure similar to the one presented here.

Next, we will demonstrate how to perform \subpage snapshot1