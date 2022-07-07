This folder contains example files to perform the svd-based mean structures shown in figure 4 of the paper.  This example will use oxDNA to run a simulation and then use the `compute_mean.py` and `compute_deviations.py` scripts to compute the average structure and per-nucleotide deviation.  There are two examples here:

 * The first is design 19, a DNA wireframe from  Jun et. al. (2019). Click [here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fsvd_mean%2F19.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fsvd_mean%2F19.top) to preview it in oxView
 * The second is a single-stranded RNA origami from Han et. al. (2017). Click [here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fsvd_mean%2Frna_rectangle.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fsvd_mean%2Frna_rectangle.top) to preview it in oxView.
 
---

1. Run an oxDNA simulation using the provided input file

   For the DNA wireframe:
   ```
   /path/to/oxDNA input_19
   ```
   For the RNA origami:
   ```
   /path/to/oxDNA input_rna
   ```
   **A couple of notes on these simulations:**
     Both of these are large structures, and as such use the CUDA implementation of oxDNA. In order to run these simulations you must be on a computer that has a GPU and the proper NVIDIA drivers. The input files here will run for 1e8 steps.  For production runs we usually recommend 1e9 steps, however 1e8 will be fine for an example. Both input files will write to the same output file names.  If you want to run both simulations, we recommend either copying them to their own directories or changing the ouput file names in the input file.

2. Compute the mean structure and the per-nucleotide deviations

   For the DNA wireframe:
   ```
   oat mean -o mean.dat -d devs.json trajectory_trap.dat
   ```
   For the RNA origami:
   ```
   oat mean -o mean.dat -d devs.json trajectory_trap.dat
   ```
   **A couple of notes on the compute_mean script:**
     
   * `-o` sets the ouput file name  
   * `-d` tells the script to automatically run compute_deviations.py using the provided filename as the output file
   * `-i` Specifies an index file, which is a space-separated list of particle IDs to perform alignment on a subset of particles
   * `-a` Specifies the configuration in the trajectory to use as a reference.  If not set a random configuraiton is used.
   * If you are running on a computer with multiple CPUs, you can additionally specify the `-p <number>` option to compute the mean structure in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus.

3. Visualize the mean structure and deviations by dragging and dropping the `mean.dat`, the topology for the structure you're studying, and `devs.json` onto an open oxView window.

   **Notes on visualization:**
   The default colormap is a red -> blue map called cooltowarm.  If you want to switch the colormap, open the developer console (`Ctrl-shift-J` on Chrome) and type `api.changeColormap('colormapName')`.  The two colormaps used in the figure are `'viridis'` and `'rainbow'`
   Pressing the `P` key will take a screencap of the current scene.
   The colorbar is dynamic, if you drag and drop multiple mean+topology+deviations file triplets, the upper and lower bounds will automatically rescale to accomodate the new data.
