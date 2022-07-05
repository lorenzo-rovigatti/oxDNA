This folder contains example files to perform angle comparison between two DNA wireframes from Han et. al. (2017). 

The structures can be previewed in oxView:
 * [Design 20](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fangles%2F20.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fangles%2F20.top)
 * [Design 23](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fangles%2F23.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fangles%2F23.top).

This example will use oxDNA to run a simulation and then use duplex_angle_finder.py to extract duplex information from the trajectory and duplex_angle_plotter.py to make a plot like Figure 7 in the paper.

1. Run an oxDNA simulation using the provided input file
   
   For design 20 (hexagonal lattice):
   `/path/to/oxDNA input_20`

   For design 23 (square lattice):
   `/path/to/oxDNA input_23`

   **A couple of notes on these simulations:** Both of these are large structures, and as such use the CUDA implementation of oxDNA.  In order to run these simulations you must be on a computer that has a GPU and the proper NVIDIA drivers. The input files here will run for 1e8 steps. For production runs we usually recommend 1e9 steps, however 1e8 will be fine for an example.

2. Fit vectors to the duplexes in the simulation

   `oat duplex_finder -o angles_20.txt input_20 trajectory_20.dat`
   
   `oat duplex_finder -o angles_23.txt input_23 trajectory_23.dat`

   **A couple notes on the duplex_angle_finder script:** This script uses the hydrogen bonding potential to determine if two nucleotides are bonded, which is unaffected by fix diffusion. You do not need to set the reference particle for this script. If you are running on a computer with multiple CPUs, you can additionally specify the `-p <number>` option to fit duplexes in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus.

4. Plot the angle between two adjacent duplexes using duplex_angle_plotter.py

   `oat duplex_angle_plotter -o angle.png -f both -n design20 design23 -i angles_20.txt 2045 1950 -i angles_23.txt 1225 1299`

   **A couple notes on the duplex_angle_plotter script:**
     The console output from this script not only gives the median, mean and standard deviation of the angle distribution, but also a representation score.  This tells you how frequently the duplex starts on the specified nucleotide.  Check a few entries in the angle file to make sure that the ID you're specifying is the duplex you want as junctions are known to shift.
     Use the selection feature of oxView to determine the nucleotide IDs to use with this script.
     * `-o` names the output files.  If the format flag is set to "both" it will append the type of graph to the end of the name.
     * `-f` specifies the format of the graph.  The options are "histogram", "trajectory" and "both".
     * `-n` names the data series in the plots. If not set it will default to particle IDs.
     * `-i` preceeds an input triplet.  This must be an angle file and two nucleotide IDs that mark the start of the duplex.  These must be the lowest number ID in the duplex.
