This directory contains and example to calculate the hydrogen bond occupancies compared with a reference design.  This example will used oxDNA to run a simulation and then use bond_analysis.py to compute the bond occupancy fraction.  The example structure used here is a yet-unpublished RNA tile design that the authors are working with, you can [click here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fh_bonds%2Frna_tile.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fh_bonds%2Frna_tile.top) to preview it in oxView.

1. Run an oxDNA simulation using the provided input file
   ```
   /path/to/oxDNA input_rna
   ```
   **A couple of notes on this simulation:**
   This simulation is small and can be run on a laptop or non-high end computer.
   The input file here will run for 1e8 steps.  For production runs we usually recomend 1e9 steps, however 1e8 will be fine for an example.

2. Generate the pairs list from the initial structure
   ```
   oat generate_force -o forces.txt -f pairs.txt input_rna rna_tile.dat
   ```
   **A couple of notes on the generator script:**
   This script handles both pairfile and forcefile generation.  When relaxing structures or joining structures together, it is frequently beneficial to hold structures together using mutual traps between the designed pairs that you then release for the production run.
   If you already have a force file containing all pairs (for example, from the Tiamat converter on TacoxDNA), you can convert that to the pairsfile used here via `oat forces2pairs`

3. Compute the hydrogen bond occupancy
   ```
   oat bond_analysis input_rna trajectory_trap.dat pairs.txt h_bonds.json
   ```
   **A couple of notes on the bond analysis script:**
   Nucleotides that have no pair defined by the design will have an occupancy of 0.
   If you are running on a computer with multiple CPUs, you can additionally specify the `-p <number>` option to calculate bonding in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus.

3. Visualize the occupancy by dragging and dropping `last_conf_trap.dat`, `rna_tile.top`, and `h_bonds.json` onto an open oxView window.

   **Notes on visualization:**
   The default colormap is a red -> blue map called `'cooltowarm'`.  If you want to switch the colormap, open the developer console (`Ctrl-shift-J` on Chrome) and type `api.changeColormap('colormapName')`.  The colormap used in the figure is `'viridis'`.
   Pressing the `P` key will take a screencap of the current scene.
   The colorbar is dynamic, if you drag and drop multiple mean+topology+occupancy file triplets, the upper and lower bounds will automatically rescale to accomodate the new data.
