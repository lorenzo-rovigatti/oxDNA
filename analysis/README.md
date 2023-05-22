# oxDNA Analysis Tools v 2.0

A suite of command line Python tools for performing generic structural analyses of oxDNA simulations.
Our goal in developing these tools is to provide a foundation of common analyses applicable to most simulations and to provide code examples for researchers looking to implement tools to meet their own research needs.
Running any script without arguments will print a brief description and list of required files.

An overarching description can be found in this paper: https://academic.oup.com/nar/article/48/12/e72/5843822.

Since the original publication all scripts have been completley rewritten to be 10-100x faster, however in the process most of the names have changed and the coding style is very different, please make sure you update your CLI autocompletes and take a look at `mean.py` for an example of how to use the new framework. 

## Dependencies and installation

OxDNA analysis tools will be automatically installed when you install oxDNA with `oxpy` Python bindings enabled.  For complete information, please see the [installation](https://lorenzo-rovigatti.github.io/oxDNA/install.html) section of the documentation.  For the impatient, these are the commands to run from the `oxDNA` root directory:

```bash
mkdir build && cd build
cmake -DPython=1 ..      #Also include -DCUDA=1 if compiling for GPU and `-DOxpySystemInstall=0` if on a personal machine or if using Python managed by Conda.
make -j4                 #-j specifies how many threads to use while compiling
make install             #installs both oxpy and oxDNA_analysis_tools
oat config               #verifies installation
```

### Numpy header error
If when get the following error when running `OAT` scripts it means that your version of Numpy differes from the one used to build the Cython portions of the code
```
numpy.ndarray size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject
```

To solve this problem, install `OAT` without pip's environment isolation by running the following in the root `oxDNA` directory:
```
python -m pip install ./analysis --no-build-isolation
```

### Setting up Bash autocompletes
The invocation `oat` is calling a Python script which then handles calling the other available scripts.  If you would like autocompletes for the specific script names (and are using a Unix command line), these are provided by `oat-completion.sh` which can also be found in the repository.  To add autocompletes to your system, either append it to your local `.bash_completion` file with:  

`cat oat-completion.sh >> ~/.bash_completion`

Or add it to your global completions with:  

`sudo cp oat-completion.sh /etc/bash_completion.d/`

Running `oat` with no arguments will list all available scripts.

### Pip installation
oxDNA Analysis Tools is also available as a standalone package through PyPi and can be installed through pip with:

`pip install oxDNA-analysis-tools`

This will also install all dependencies.  Bash autocompletes will not be set up, see below for setting up autocompletion.

**Installation from PyPi will only work on Linux and OSX.  If you want to use OAT on Windows, please install from source (or install on WSL)**

Please note that `oxpy` is a dependency for any script that needs to interface with the oxDNA energy function and therefore installing `OAT` in this way is not recommended.

### Build and installation from source
`OAT` can also be installed from the [GitHub repository](https://github.com/lorenzo-rovigatti/oxDNA) or the zip file of the source code available on PyPi via the following method:  

1. Clone the repository or download and inflate the zip file.  
2. Navigate to the analysis directory:
   `cd analysis`
2. Build the Cython code for your system:  
   `python -m build`
3. Run one of the following commands (pip to automatically install dependencies or setup.py if you would like to manage them yourself):  
   `python -m pip install .`  
   `python setup.py install`  

If you are not installing via pip, the following dependencies are required and can all be obtained from either pip or conda:  
[Python](https://www.python.org/): >=3.9,<br/>
[oxDNA compiled with oxpy](https://github.com/lorenzo-rovigatti/oxDNA): (minimum version >=3.2.2)<br/>
[NumPy](https://numpy.org/): >=1.16,<br/>
[MatPlotLib](https://matplotlib.org/index.html): >=3.0.3 (minimum version 3.0),<br/>
[BioPython](https://biopython.org/): >=1.73,<br/>
[Scikit-Learn](https://scikit-learn.org/stable/): >=0.21.2,<br/>

### Test your installation
To check your installation run:  
   `oat config`

-------------------------------------------------------------------

## Using oxDNA analysis tools
Once installed, all standalone scripts can be called from the command line via the following invocation:  
`oat <script name> <script arguments>`  

For example, to compute the mean structure and deviations of a file called `trajectory.dat` using 4 CPUs and outputting to files called `mean.dat` and `devs.json`, you would run:  
`oat mean -p 4 -o mean.dat -d devs.json trajectory.dat`

To see a detailed description of the script command line arguments, run the script with the `-h` flag.

For a full list of scripts and arguments, please see the [documentation](https://lorenzo-rovigatti.github.io/oxDNA/oat/cli.html).

These scripts are intended to be extensible and re-used for custom analysis by users.  The functions in this library can be imported into your Python scripts via:  
`from oxDNA_analysis_tools.<script name> import <object name>`

So for example, if you would like to use the file reader, you would include the following in your imports:  
`from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs`

For a complete list of importable functions, please see the [documentation](https://lorenzo-rovigatti.github.io/oxDNA/oat/index.html#scripting-interface)

-------------------------------------------------------------------

## Brief script descriptions

Running instructions can be obtained for all scripts by running them with no arguments or the -h flag.

 * `align (-p <n_cpus> -i <index_file> -r<reference_conf_file>) <trajectory_file> <output_file>` Aligns all configurations in the trajectory to the first configuration.  Produces an oxDNA trajectory that is a copy of the provided trajectory with all translations and rotations removed.<br/>
 * `anm_parameterize <index_file> <mean_file> <trajectory_file> <output_file>` Computes the fluctuations for 'superparticles' made up of multiple nucleotides defined by lists of particles in the index file.<br/>
 * `backbone_flexibility (-p <n_cpus> -o <output_file> -d <data_file>) <trajectory> <topology>` Produces a Ramachandran plot of backbone torsions vs dihedrals for the given trajectory.<br/>
 * `bond_analysis (-p <n_cpus>) <input> <trajectory> <designed_pairs_file> <output_file>`  Calculates the hydrogen bond occupancy compared with the intended design.  Produces an oxView json file to `<output_file>` that contains a color overlay corresponding to the occupancy of each pair observed in the simulation.<br/>
 * `centroid (-p <n_cpus> -o <centroid file> -i <index file>) <reference_structure> <trajectory>` Takes a reference structure (usually a mean structure) and a trajectory and returns the structure in the trajectory with the lowest RMSD to the reference as an oxDNA configuration file. This is not a true centroid which has the lowest RMSD to **all** structures, but it's a good approximation.<br/>
 * `clustering (-e <eps> -m <min_points>) <serialized data input>` Takes a set of configuration coordinates from other scripts and performs a DBSCAN clustering.  Produces trajectory files for each cluster (note that these trajectories do not necessarily contain contiguous timesteps) and a visual representation of the clusters in either a 2D or 3D plot (2D if only 2 dimensions of input data are given, 3D if 3 or more are given with the first three displayed).  The -c option on pca.py and distance.py will call this script. Clustering.py serializes its own data to a file called `cluster_data.json` so you can re-launch the script to modify clustering parameters without re-running the analysis.<br/>
 * `config (-n <chunk_size>)` Performs dependency checks and sets system parameters. The chunk size determines how many configurations are loaded into memory at a time (ncpus * chunk_size configurations will be loaded).  For small systems this number can be increased to improve performance.<br/>
 * `contact_map (-p <n_cpus> -g <graph_file> -d <data_file>) <trajectory>` produces a contact map of internucleotide distances.<br/>
 * `db_to_force (-o <output_file> -s <stiffness>) <db_file>` Convert a dot-bracket file to an oxDNA mutual trap file.  Useful for assembling small, unknotted structures from NUPACK predictions.<br/>
 * `deviations (-i<index_file> -p <n_cpus> -o <deviations_file> -r <rmds_plot> -d <order_parameter_file>) <mean_structure> <trajectory file>` Computes the per-nucleotide RMSF and per-configuration RMSD from the mean structure.  Can be called automatically by `mean` with the -d option. Produces an oxView json flie that colors each particle based on its RMSF and an oxView order parameter file which overlays the RMSD on the trajectory.<br/>
 * `distance (-o <graph_file> -f <histogram/trajectory/both> -d <data_output_file> -n <data series names> -p <n_cpus> -c) -i <<trajectory> <particleID 1> <particleID 2> (<particleID 1> <particleID 2> ...)>` Computes the distance between provided particle pairs. The `-i` option can be called multiple times to overlay data from multiple trajectories.  Additional calls will be overlaid on the same graph. Produces the user's choice of histograms, timeseries or text outputs (set by `-f` and `-d` options).  The `-c` option will run the output of the distances through the `clustering` script<br/>
* `duplex_angle_plotter (-o <output> -f <histogram/trajectory/both> -d <data_output_file> -n <data series names>) -i <<angle file> <particleID 1> <particleID 2> (<particleID 1> <particleID 2> ...)>` Reads the angle file produced by `duplex_finder` and produces either histograms or timeseries of the angle between specified duplexes.  The `-i` option can be called multiple times to overlay data from multiple trajectories.  Additional calls will be overlaid on the same graph.<br/>
* `duplex_angle_finder (-p <n_cpus> -o <output>) <input> <trajectory>` Produces an angle file containing identification information for all duplexes at each configuration in the trajectory.  This file is visualized by `duplex_angle_plotter`<br/>
* `forces2pairs -o <output_file> <force file>` Takes an oxDNA external forces file (which can be generated for all pairs in a structure using `generate_force` or with oxView) and produces a pairs file.  This output is used as an input for `bond_analysis`. <br/>
* `generate_force (-o <output> -f <pairs file> -s <stiffness>) <input> <configuration>` Produces an external force file enforcing the current base pair arrangement. The `-f` option will automatically call `forces2pairs`.<br/>
 * `mean (-p <n_cpus> -o <mean_file> -d <deviations_file> -i <index_file> -a <align_conf_id>) <trajectory file>` Computes the mean structure of a trajectory using single value decomposition.  If the `-i` flag is added with an index file containing a list of particle IDs, the mean structure will be calculated based only on the subset of particles included in the list.  These lists can be created from selected bases in oxView using the "Download selected base list" button.  By default, this script aligns to a random configuration in the trajectory.  However, if you would like to align to a specific configuration, you can specify its position in the trajectory with the `-a` flag.  The `-d` flag will automatically run `deviations` from the mean structure.<br/>
 * `minify (-p <n_cpus> -d <precision> -a) <trajectory> <output_file>` Reduces the file size of a trajectory by dropping the velocities.  If called with `-a`, the a1 and a3 vectors will also be dropped.  If called with `-d`, the position, a1 and a3 vectors will be truncated to the specified precision.<br/>
 * `multidimensional_scaling_mean (-p <n_cpus> -o <mean_structure_file> -d <deviations_file>) <trajectory>` Computes the mean structure based on local pairwise distances between particles.  An alternative to `mean` that works better for highly flexible structures.  Produces an oxDNA configuration and an oxView json file showing the per-particle deviation in distance to neighboring particles.<br/>
 * `output_bonds (-v <output> -u <oxDNA/pNnm> -p <n_cpus>) <input> <trajectory>` Lists all the interactions between nucleotides.  The output is the same as the `pair_energy` observable from oxDNA.  The `-v` option will instead create an oxView color overlay with the average energy per nucleotide for each potential.<br/>
 * `oxDNA_PDB topology configuration_file direction  <list of protein pdb files in system in order of occurence> (-o <output_file> -r <rmsf_file> -H -u -1)`  Converts either oxDNA DNA files to pdb, or oxDNA DNA/Protein hybrids to pdb format. Note that the PDB files must contain only the proteins, even if the structure was generated from a file which also contained DNA. The `-H` flag will add hydrogens to the PDB file, the `-u` flag will use uniform residue names in the PDB file, the `-1` flag will generate a separate PDB file for each strand in the structure.  Does not work with RNA structures.<br/>
 * `pca (-p <n_cpus> -c) <trajectory> <mean file> <output>` Computes the principal components of deviations away from the mean.  The principal components are written as an oxView json overlay file that will show the direction of the top mode.  More components can be summed and added to the overlay by modifying the `N` variable in the script.  If the `-c` flag is used it will also run the clustering script on each configuration's position in principal component space. <br/>
 * `persistence_length (-p <n_cpus> -d <data_file> -n <plot_name>) <trajectory> <input> <nucid1> <nucid2>` Calculate the persistence length of a continuously paired segment of DNA. <br/>
 * `plot_energy (-o <output_file> -f <histogram/trajectory/both>) <energy_file>` Plot the energy file from an oxDNA simulation <br/>
 * `subset_trajectory (-p <n_cpus>) -i <index_file output_file> <trajectory> <topology>` Split a the trajectory of a configuration into trajectories containint only a subset of the particles given in the index file.  `-i` can be called multiple times to produce multiple subsets from the same trajectory.
 * `superimpose (-i <index_file>) <configuration> <configuration> (<configuration> <configuration> ...)` Superimposes all further configurations onto the first configuration file provided.  It is expected that all referenced nucleotides are the same, so the configurations either must share the same topology, or you must only align to the shared particles using an index file. An index file can be downloaded from oxView using the "Download selected base list" button.
 
### UTILS
The UTILS directory contains utility modules used by other scripts in this package.

* `chunksize.py` Sets the size of chunks read in by the `oat_multiprocesser`.  The chunk size determines how many configurations are loaded into memory at a time (ncpus * chunk_size configurations will be loaded).  For small systems this number can be increased to improve performance. <br/>
* `data_structures.py` Contains definitions for common data structures used in the scripts.  Includes definitions such as `Trajinfo`, `TopInfo`, and `System`
* `dd12_na.pdb` Used during pdb conversion script
* `geom.py` A set of algorithms to find various geometric parameters of DNA/RNA helices.<br/>
* `oat_multiprocessor.py` Parallelization method which uses partial function composition to distribute functions and configuration blocks to processors and accumulate the results.
* `pdb.py` Helper Functions/Classes for pdb conversion <br/>
* `protein_to_pdb` Contains protein specific functions for protein to pdb conversion<br/>
* `RyeReader.py` File handling functions. `describe` is used to extract metadata from oxDNA files while `get_confs` is used to read trajectories.
* `utils.py` Contains utility functions for pdb conversion<br/>

### Cython utils
The file reader is written in Cython for maximum speed.  Did you know the built-in Python file readers are **really** slow??
* `copy_build.sh` Local build pipeline used for development.
* `get_confs.pyx` Optimized file reader for oxDNA trajectories
* `get_confs.c` get_confs transpiled to C.  Distributed for compatibility purposes.
* `setup.py` Cython build instructions.

### External Force Utils (WIP)
The external_force_utils directory contains function definitions for working with external forces.

* `exclude_force <force_file> <index_file>` A script which removes forces on the particles specified in the index file.
* `force_reader` Reader functions for reading oxDNA external force files.
* `forces` Definitions for all force types used by oxDNA.

## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

By default scripts in this package that produce graphs save them as .png files.  All graphing is done using the Matplotlib interface and users are encouraged to make modifications to the graph styles to fit their unique needs.
   
## File formats

This package mostly uses the oxDNA files as described in [the oxDNA documentation](https://dna.physics.ox.ac.uk/index.php/Documentation).  A brief description of each file is provided here for easy reference:  
* **trajectory** - A file containing a sequence of oxDNA configurations.  Each configuration starts with a three line header containing the timestep, box size and energy information.  There is then one line per particle with 15 values corresponding to the position, orientation and velocity of each particle. 
* **topology** - A file containing sequence and connectivity information of the simulated structure.  The first line defines the number of particles and the number of strands.  There is then one line per particle with 4 values corresponding to the strand ID, the base type, the 3' neighbor and the 5' neighbor of each particle.  Note that oxDNA files are written 3'-5' rather than the traditional 5'-3'.  
* **input** - The input file used to run oxDNA.  This contains simulation information such as number of steps, simulation method and temperature as well as I/O information.  Example files can be found in the "example_input_files" and "paper_examples" directories.  
* **force file**: An oxDNA mutual trap file that defines an external force between two particles in a simulation.  This is also defined in [the oxDNA documentation](https://dna.physics.ox.ac.uk/index.php/Documentation).  

The Following files are unique to this package:  
* **oxView json file**: This file contains overlay information that can be read by [oxView](https://github.com/sulcgroup/oxdna-viewer).  There are two different formats that are produced by these scripts.  The first is a one-value-per-particle file that creates a colormap overlay with extreme colors corresponding to the minimum and maximum values in the file.  The second is a three-values-per-particle file that oxView uses to draw arrows on the scene.  OxView automatically differentiates files based on the number of values corresponding to each particle.  
* **designed pairs file**: This file contains a list of particle pairs in the intended design.  Each line corresponds to a pair and each pair is a space-separated pair of particle IDs.  Designed pairs files can be generated by `forces2pairs` and `generate_force`.  
* **angle file**: The output file generated by `duplex_finder`.  Details on the format can be found in a comment in the `duplex_angle_plotter` script, but briefly each line contains starting and ending nucleotides and orientation data for a duplex in the structure.  Like trajectories, this contains information for every configuration in a trajectory.  
* **index file**: A space-seperated list of particle IDs used for subset alignment.  It can be generated by the "Download Selected Base List" button in oxView.  
* **serialized data input**: To make it easy to adjust clustering parameters, the clustering script serializes its input in json format so the script can be re-launched quickly with this file as the only input. 

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
Erik Poppleton, Joakim Bohlin, Michael Matthies, Shuchi Sharma, Fei Zhang, Petr Sulc: Design, optimization, and analysis of large DNA and RNA nanostructures through interactive visualization, editing, and molecular simulation. (2020) Nucleic Acids Research e72. https://doi.org/10.1093/nar/gkaa417
