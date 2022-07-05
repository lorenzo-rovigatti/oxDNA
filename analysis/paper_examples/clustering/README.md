This folder contains example files to perform unsupervised clustering using DBSCAN on different structural states visited during a simulation of a single-stranded RNA tile (ongoing work not yet published) as seen in figure 11.  As the simulation is very long, difficult to execute, and still under analysis, we provide a trajectory file instead of having the user run the simulation as in the other examples. Click [here](https://sulcgroup.github.io/oxdna-viewer/?configuration=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fclustering%2Ftrajectory.dat&topology=https%3A%2F%2Fraw.githubusercontent.com%2Fsulcgroup%2Foxdna_analysis_tools%2Fmaster%2Fpaper_examples%2Fclustering%2Frna_tile.top) to view the trajectory in oxView.

Clustering can be done on any data set that can position each configuration to a point in n-dimemsional space.  In this package, it is currently implemented for distances between particles and for linear combinations of principal components.  In this example, we will use the principal components as this size of system is amenible to calculation of components and it is the most information-rich example.

1. Compute the mean structure of the trajectory to use as a reference for pca.

   `oat mean -f oxDNA -o mean_all.dat trajectory.dat`

If you have multiple CPUs available, this can be run with the `-p <n_cpus>` option to speed up computation.  We do not need to calculate deviations for this example.

2. Compute principal components and pass the output to the clustering algorithm

   `oat pca -c input_rna trajectory.dat mean_all.dat pca.json`

The `-c` option will tell the script to consider every configuration as a linear combination of the principal components and then calculate clusters based on the weights of each component.

This will produce a number of files.  There are three produced by the pca script:
 * scree.png - the standard scree plot showing the relative weights of each component
 * coordinates.png - A plot showing reconstructions in 3D space without clustering.  Can tell you if a particular simulation is amenible to clustering.
 * pca.json - An oxView overlay showing the first principal component as arrows.

The next two are for information about clustering:
 * cluster_data.json is the output from pca.py serialized to make it faster to adjust clustering parameters.  If clustering fails for a design, adjust the eps and min_samples parameters in clustering.py and re-run it with:

   `oat clustering cluster_data.json`

 * animated.mp4 - a video of coordinates.png spinning around, now with colors corresponding to the clusters.  If you would rather view the plot interactivley, uncomment the code block relating to plotting 3D plots in clustering.py

Finally, the clusterer will split the trajectory file into the respective clusters.  In this example the first three clusters correspond to the ones in figure 11.  The remaining two clusters are entirely linear structures that we do not believe are physically relevant.  In the figure, each cluster was then run through multidimensional_scaling_mean.py, bond_analysis.py, and duplex_angle_finder.py to obtain the data shown.
