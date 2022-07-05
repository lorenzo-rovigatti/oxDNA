This directory contains an example to analyze distances between the mobile and immoble units of the tethered multiflourophore (tmf) structure from Schickinger, Zacharias, and Dietz (2018). 

This example will use oxDNA to run simulations of the [open](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fdistances%2Ftmf_open.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fdistances%2Ftmf.top) and [closed](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fdistances%2Ftmf_closed.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fdistances%2Ftmf.top) states and then generate graphics showing the difference in spacing between the two states.

1. Run oxDNA simulations of the open and closed states using the provided input files
    ```
    /path/to/oxDNA input_closed
    ```

    ```
    /path/to/oxDNA input_open
    ```

   **A couple of notes on this simulation:**
     This is a large structure, and as such use the CUDA implementation of oxDNA.  In order to run this simulation you must be on a computer that has a GPU and the proper NVIDIA drivers.
     The input file here will run for 1e8 steps.  For production runs we usually recommend 1e9 steps, however 1e8 will be fine for an example.

2. Compute the distance between the nucleotides at the end of the tether.

   ```
   oat distance -o distance.png -f both -n closed open -i trajectory_closed.dat 6453 9630 -i trajectory_open 6453 9630
   ```

   A couple notes on the distance script:
     * `-o` names the output files.  If the format flag is set to "both" it will append the type of graph to the end of the name.
     * `-f` specifies the format of the graph.  The options are "histogram", "trajectory" and "both"
     * `-i` preceeds an input set.  This must be an oxDNA input file, a trajectory and the pairs of particle ids to calculat the distance between.
     * `-d` will dump the data used to construct the graphs as a plaintext file with the specified name
     * `-n` will name the data series on the plots, otherwise it defaults to particle ids.
     * If running on only a single trajectory, the `-c` option will run the clustering algorithm on the distance data and produce separate clusters with similar distances between the provided particles.
