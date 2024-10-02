# oxDNA Analysis Tools v 2.1

OxDNA analysis tools (`oat`) is both a command line program as well as an importable Python library for structural analysis of oxDNA simulation trajectories.

The documentation for using `oat` as either a command-line-interface or as an imported Python library can be found [here](https://lorenzo-rovigatti.github.io/oxDNA/oat/index.html)

An overarching description can be found in this paper: https://academic.oup.com/nar/article/48/12/e72/5843822.

Note, however, that since the original publication all scripts have been completley rewritten to be 10-100x faster, however in the process many of the names have changed and the coding style is very different, please make sure you update your CLI autocompletes and take a look at `mean.py` for an example of how to use the new framework. 

## Dependencies and installation

OxDNA analysis tools will be automatically installed when you install oxDNA with `oxpy` Python bindings enabled.  For complete information, please see the [installation](https://lorenzo-rovigatti.github.io/oxDNA/install.html) section of the documentation.  For the impatient, these are the commands to run from the `oxDNA` root directory:

```bash
mkdir build && cd build
cmake -DPython=1 ..      #Also include -DCUDA=1 if compiling for GPU and `-DOxpySystemInstall=1` if on a personal machine or if using Python managed by Conda.
make -j4                 #-j specifies how many threads to use while compiling
make install             #installs both oxpy and oxDNA_analysis_tools
oat config               #verifies installation
```

`oat` can also be installed (or updated) separate from oxDNA via `pip` by running the following in the `oxDNA` directory:
```
python -m pip install ./analysis
```

The following dependencies are required by `oat`:<br/>
[Python](https://www.python.org/): >=3.9,<br/>
[oxDNA compiled with oxpy](https://github.com/lorenzo-rovigatti/oxDNA): (minimum version >=3.2.2)<br/>
[NumPy](https://numpy.org/): >=1.16,<br/>
[MatPlotLib](https://matplotlib.org/index.html): >=3.0,<br/>
[Scikit-Learn](https://scikit-learn.org/stable/): >=0.21,<br/>

### Numpy header error
If when get the following error when running `oat` scripts it means that your version of Numpy differs from the one used to build the Cython portions of the code
```
numpy.ndarray size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject
```

To solve this problem, install `oat` without pip's environment isolation by running the following in the root `oxDNA` directory:
```
python -m pip install ./analysis --no-build-isolation
```

### Setting up Bash autocompletes
The invocation `oat` is calling a Python script which then handles calling the other available scripts.  If you would like autocompletes for the specific script names (and are using a Unix command line), these are provided by `oat-completion.sh` which can also be found in the repository.  To add autocompletes to your system, either append it to your local `.bash_completion` file with:  

`cat oat-completion.sh >> ~/.bash_completion`

Or add it to your global completions with:  

`sudo cp oat-completion.sh /etc/bash_completion.d/`

Running `oat` with no arguments will list all available scripts.

### Build and installation from source 
If you make changes to the Cython files in `src/oxDNA_analysis_tools/cython_utils`, you must rebuild the library to have them available.  This requires `Cython`, `setuptools`, `numpy` and `python-build`.

1. Navigate to the analysis directory:
   `cd analysis`
2. Build the Cython code for your system:  
   `python -m build`
3. Install `oat`:  
   `python -m pip install .`  

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

## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

By default scripts in this package that produce graphs save them as .png files.  All graphing is done using the Matplotlib interface and users are encouraged to make modifications to the graph styles to fit their unique needs.
   
## File formats

This package mostly uses the oxDNA files as described in [the oxDNA documentation](https://lorenzo-rovigatti.github.io/oxDNA/configurations.html).  A brief description of each file is provided here for easy reference:  
* **trajectory** - A file containing a sequence of oxDNA configurations.  Each configuration starts with a three line header containing the timestep, box size and energy information.  There is then one line per particle with 15 values corresponding to the position, orientation and velocity of each particle. 
* **topology** - A file containing sequence and connectivity information of the simulated structure.  The first line defines the number of particles and the number of strands.  There is then one line per particle with 4 values corresponding to the strand ID, the base type, the 3' neighbor and the 5' neighbor of each particle.  Note that oxDNA files are written 3'-5' rather than the traditional 5'-3'.  
* **input** - The input file used to run oxDNA.  This contains simulation information such as number of steps, simulation method and temperature as well as I/O information.  Example files can be found in the "example_input_files" and "paper_examples" directories.  
* **force file**: An oxDNA mutual trap file that defines an external force between two particles in a simulation.  This is also defined in [the oxDNA documentation](https://lorenzo-rovigatti.github.io/oxDNA/forces.html).  

The Following files are unique to this package:  
* **oxView json file**: This file contains overlay information that can be read by [oxView](https://github.com/sulcgroup/oxdna-viewer).  There are two different formats that are produced by these scripts.  The first is a one-value-per-particle file that creates a colormap overlay with extreme colors corresponding to the minimum and maximum values in the file.  The second is a three-values-per-particle file that oxView uses to draw arrows on the scene.  OxView automatically differentiates files based on the number of values corresponding to each particle.  
* **designed pairs file**: This file contains a list of particle pairs in the intended design.  Each line corresponds to a pair and each pair is a space-separated pair of particle IDs.  Designed pairs files can be generated by `forces2pairs` and `generate_force`.  
* **angle file**: The output file generated by `duplex_finder`.  Details on the format can be found in a comment in the `duplex_angle_plotter` script, but briefly each line contains starting and ending nucleotides and orientation data for a duplex in the structure.  Like trajectories, this contains information for every configuration in a trajectory.  
* **index file**: A space-seperated list of particle IDs used for subset alignment.  It can be generated by the "Download Selected Base List" button in oxView.  
* **serialized data input**: To make it easy to adjust clustering parameters, the clustering script serializes its input in json format so the script can be re-launched quickly with this file as the only input. 

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
Erik Poppleton, Joakim Bohlin, Michael Matthies, Shuchi Sharma, Fei Zhang, Petr Sulc: Design, optimization, and analysis of large DNA and RNA nanostructures through interactive visualization, editing, and molecular simulation. (2020) Nucleic Acids Research e72. https://doi.org/10.1093/nar/gkaa417  
Erik Poppleton, Michael Matthies, Debesh Mandal, Flavio Romano, Petr Sulc, Lorenzo Rovigatti: oxDNA: coarse-grained simulations of nucleic acids made simple. (2023) Journal of Open Source Software 8(81), 4693. https://doi.org/10.21105/joss.04693
