This directory contains an example to produce a mean structure based on local particle distances reconstructed with multidimensional scaling as shown in figure 5 of the paper.  This example will use oxDNA to run a simulation and then use multidimensional_scaling_mean.py to compute the average structure and the per-nucleotide local distance deviation.  The example here is the same RNA structure from Han et. al. (2017) used in the SVD mean example, click [here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fmds_mean%2Frna_rectangle.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fmds_mean%2Frna_rectangle.top) to preview it in oxView

1. Run an oxDNA simulation using the provided input file
   ```
   /path/to/oxDNA input_rna
   ```
   **A couple of notes on this simulation:**
     This is a large structure, and as such use the CUDA implementation of oxDNA.  In order to run this simulation you must be on a computer that has a GPU and the proper NVIDIA drivers. The input file here will run for 1e8 steps.  For production runs we usually recommend 1e9 steps, however 1e8 will be fine for an example.

2. Compute the mean structure and the per-nucleotide deviations
   ```
   oat multidimensional_scaling_mean -o mean.dat -d devs.json trajectory_trap.dat
   ```
   **A couple of notes on the multidimensional_scaling_mean script:**
     Multidimensional scaling fails to compute solutions for structures larger than a few thousand nucleotides.  This structure, at 2,018 nucleotides, is the largest that has been successful.  It does not work above 10,000.  We are interested in feedback if users find this script works on structures larger than this one.
     If you are running on a computer with multiple CPUs, you can additionally specify the `-p <number>` option to compute the mean structure in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus.

3. Visualize the mean structure and deviations by dragging and dropping the `mean.dat`, `rna_rectangle.top`, and `devs.json` onto an open oxView window.

   **Notes on visualization:**
   The default colormap is a red -> blue map called `'cooltowarm'`.  If you want to switch the colormap, open the developer console (`Ctrl-shift-J` on Chrome) and type `api.changeColormap('colormapName')`.  The colormap used in the figure is 'viridis'.
   Pressing the "P" key will take a screencap of the current scene.  
   The colorbar is dynamic, if you drag and drop multiple mean+topology+deviations file triplets, the upper and lower bounds will automatically rescale to accomodate the new data.

4. In figure 6b, the svd and mds means are shown overlaid.  To reproduce this figure, compute the mean structure using both scripts and then use `oat superimpose` to generate the alignment:
   ```
   oat superimpose mean.dat ../svd_mean/mean.dat
   ```
This will generate a file called `aligned0.dat` which is the svd mean structure aligned to the mean structure generated via mds.
