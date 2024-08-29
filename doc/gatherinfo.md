Step 1: gather information for both snapshot and trajectory modeling {#gatherinfo}
====================================

# Gathering information

In the first step, information about the process of interest is gathered. This information can theoretically be applied to either modeling static snapshots or modeling trajectories, and can be used for representing the model, scoring the model, sampling alternative models, filtering sampled models, or validating the models.

For this tutorial, we used the X-ray crystal structure of the complete Bmi1/Ring1b-UbcH5c complex (a), synthetically generated electron tomography (ET) density maps during the assembly process (b), synthetically generated protein copy numbers during the assembly process, which can be calculated from experiments such as fluorescence correlation spectroscopy (FCS) (c), and synthetically generated small-angle X-ray scattering (SAXS) profiles during the assembly process (d). The crystal structure of the complex informs the final state of our model as well as the structure of the individual proteins. The time-dependent ET and SAXS data give two inputs to inform the size and shape of the assembling complex. The protein copy number data informs the stoichiometry of the complex during assembly.

\image html Input.png width=600px

These pieces of information are stored in the `Input_Information` folder. In addition to containing the raw data used for the tutorial, this folder contains the code necessary to generate the synthetic data. This code is described in `README` files in each directory, but is not the focus of our tutorial.

The `FASTA` folder contains `3rpg.fasta.txt`, which provides the sequence information for each protein in the Bmi1/Ring1b-UbcH5c complex. The `PDB` folder contains the PDB structure for the fully assembled Bmi1/Ring1b-UbcH5c complex, [3RPG](https://www.rcsb.org/structure/3rpg).

The `gen_FCS` folder contains protein copy number data for each protein as a function of time. Our code will use the `exp_comp{prot}.csv` files, where {prot} is the protein corresponding to that copy number data. Each csv file has 3 rows, which correspond to the time at which the data was taken ("Time"), the mean protein copy number at that time ("mean"), and the standard deviation in protein copy number at that time ("std").

The `ET_data` folder contains the time-dependent ET data. Briefly, at each time point, a subset of Bmi1/Ring1b-UbcH5c complex proteins were used to compute a density map at each time point, and then random noise was added to this true density profile. The results of this computation are stored as `add_noise/{time}_noisy.mrc` and `add_noise/{time}_noisy.gmm`, where {time} is the time point in which the time dependent ET data was calculated.

The `gen_SAXS` folder contains the time-dependent SAXS data. Experimental SAXS profiles are forward profiles calculated from the true structure by [FoXS](https://modbase.compbio.ucsf.edu/foxs/), and are stored as `{time}_exp.dat`, where {time} is the time point in which the time dependent ET data was calculated.

In addition to the four types of data used here, a variety of data could be useful for the spatiotemporal modeling of protein complexes. IMP currently features  restraints for a variety of experimental data or prior models, including chemical cross-links, Förster resonance energy transfer, comparative structural models, and deep-learning structural models, all of which could inform spatiotemporal modeling through a procedure similar to the one presented here.

Next, we will demonstrate how to perform [Snapshot modeling steps 2-4: representation, scoring, and search process](@ref snapshot1).

