This directory contains an example to compute the principal components of a holliday junction as shown in figure 10 of the paper. Click [here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2FPCA%2Fholliday.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2FPCA%2Fholliday.top) to preview the structure in oxView

1. Run an oxDNA simulation using the provided input file.

   `/path/to/oxDNA input_dna`

   **A couple of notes on this simulation:** This simulation is small and can be run on a laptop or non-high end computer.
     The input file here will run for 1e8 steps.  For production runs we usually recomend 1e9 steps, however 1e8 will be fine for an example. This simulation uses the sequence-dependent model of DNA.  By default, oxDNA uses an average sequence model, but in cases where the different strengths of interactions between the four nucleobases may be important, it is more appropriate to use sequence dependence (see Snodin et. al. (2015)).

2. Compute the mean structure via SVD

   `oat mean -o mean.dat trajectory_trap.dat

   PCA is defined by deviations from some reference.  In this case, we will use the mean structure as a reference.

3. Compute the principal components

   `oat pca input_dna trajectory_trap.dat mean.dat pca.json`

   **A few notes on the PCA script:**
   If you are running on a computer with multiple CPUs, you can additionally specify the `-p <number>` option to compute principal components in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus
   The `-c` option will run the clustering algorithm on reconstructions of each configuration in principal component space.  By default this uses all components, but if you want to use only the top few components, uncomment the line that truncates the linear terms in the last block of code in the script.

This will produce an oxView json file that plots arrows on the structure which corresond to the weighted sum of the first n components.  n can be set by modifying the SUM variable in the pca.py script.  To view in the viewer, drag and drop the mean structure, topology and PCA json files onto the viewer window. 

This script also creates a scree plot showing the weights of components.
