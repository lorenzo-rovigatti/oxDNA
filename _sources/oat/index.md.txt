# oxDNA Analysis Tools

oxDNA Analysis Tools (*oat* in short) is a suite of command line Python tools for performing generic structural analyses of oxDNA simulations. The package functions both as a series of standalone command-line scripts and as importable modules for building custom analysis tools.

## Command Line Interface

The command line version of the scripts can be run via:
```
oat <script name> <script arguments>
```

There are bash autocompletes avilable for the script names, which can be activated by copying the file `/oxDNA/analysis/oat-completion.sh` to your autocompletes file (generally `~/.bash_completion`) and restarting your command line session.

Documentation for individual scripts:

```{eval-rst}
.. toctree::
   :maxdepth: 2

   cli.md
```

## Scripting interface

The scripting API for `oat` can be imported into Python scripts through the `oxDNA_analysis_tools` package.  For example, to use the optimized oxDNA trajectory reader you would include the following lines in your script:

```python
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs

# get the top and traj names from the command line or hardcode them

top_info, traj_info = describe(top, traj)
confs = get_confs(top_info, traj_info, start_conf, n_confs)
```

Full API documentation:

### Top-level API

```{eval-rst}
.. toctree::
   :maxdepth: 2
   
   api.md
```

### Utility API
```{eval-rst}
.. toctree::
   :maxdepth: 2

   utils.md
```

### Forces API
```{eval-rst}
.. toctree::
   :maxdepth: 2

   forces.md
```

## Analysis notebooks
As simulations become more prevalent and analysis pipelines more complicated, it can be beneficial for researchers to define their entire simulation/analysis pipeline in Python notebooks, much as is often done in the machine learning community.  The modular nature of `oxpy` and `oxDNA_analysis_tools` lends itself to composing analyses together to ask specific scientific questions.

This section represents a minimal example of how to run, analyze and visualize a simulation from within a [Jupyter Notebook](https://jupyter.org/).  For a more comprehensive example which includes live plotting of values while the simulation is running, see the [OXPY_JUPYTER](https://github.com/lorenzo-rovigatti/oxDNA/blob/master/examples/OXPY_Jupyter/literate_sim.ipynb) example in the `oxDNA/examples` directory.

From a directory containing an oxDNA input file, and initial configuration/topology files:

**First, we will set up some functions for running simulations with oxpy**
Instead of modifying an existing input file, you can also simply provide the entire input file as a `dict` to kwargs and remove the `init_from_filename` call.
```python
import oxpy
import multiprocessing

# Running oxpy multiple times from the same jupyter kernel crashes the kernel for some reason
# So we will be starting simulations from a separate thread.
def spawn(f, kwargs = {}):
    p = multiprocessing.Process(target=f, kwargs = kwargs)
    p.start()
    return p

# Sets up and runs the simulation
def run(**kwargs):
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename("input")

        # Modify inputfile at runtime
        # all input parameters provided to oxDNA need to be strings 
        for k,v in kwargs.items(): 
            inp[k] = v
        
        manager = oxpy.OxpyManager(inp)
        manager.run_complete()

p = None
```

**Next, we will initiate a simulation.**  
For this example, we set the number of steps and print intervals to be smaller such that this simulation finishes in a reasonable amount of time.

```python
# Kill any currently running simulations
try:
    p.terminate()
except:
    pass

# You can pass modifications to a default input file as dict of key : value pairs
# Both key and values must be strings
input_mods = {
    "steps" : "1e6",
    "print_conf_interval" : "1e4",
    "print_energy_every" : "1e4"
}

# Run a simulation and obtain a reference to the background process
p = spawn(run, input_mods)

# The input file used to run the simulation gets lost in the subprocess, 
# so make it again so we have access to it for analysis
with oxpy.Context():
    inp = oxpy.InputFile()
    inp.init_from_filename("input")
    for k,v in input_mods.items(): 
        inp[k] = v
```

**You can stop a running simulation by terminating its process:**
```python
# You can stop the simulation at any time by running this cell.
p.terminate()
```

**Use `oat` to compute the mean structure and RMSF**
```python
from oxDNA_analysis_tools.mean import mean
from oxDNA_analysis_tools.deviations import deviations
from oxDNA_analysis_tools.UTILS.RyeReader import describe

# Get the trajectory and topology file names from the oxpy.InputFile object
top = inp["topology"]
traj = inp["trajectory_file"]
top_info, traj_info = describe(top, traj)

# Compute the mean structure and RMSFs
mean_conf = mean(traj_info, top_info, ncpus=4)
RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=4)

#They come out as numpy arrays, need to be a dict with a list for visualization
RMSFs = {"RMSF": RMSFs.tolist()}
```

**oxView iframes can be embedded in Jupyter Notebooks with `oat`'s oxView lib**
```python
from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf

# Display python objects in an oxview iframe
oxdna_conf(top_info, mean_conf, RMSFs)
```

Which will produce an oxView iframe:
<script>
    function handle(){
        inbox_settings = ["Monomer", "Origin"]
        frame_id = "1"
        let t_files = ["../_static/ico.top", "../_static/mean.dat", "../_static/devs.json"];
        let t_blobs = []
        for (let i = 0; i < t_files.length; i++){
            let f = new XMLHttpRequest();
            f.open("GET", t_files[i], false);
            f.onreadystatechange = function () {
                t_blobs.push(new Blob([f.responseText], {type : 'text/plain'}));
            }
            f.send(null)
        }
        let t_ext = ["top", "dat", "json"];
        const frame = document.getElementById('oxview-frame-1');
        frame.contentWindow.postMessage({message : 'iframe_drop',files: t_blobs, ext: t_ext, inbox_settings : inbox_settings}, "https://sulcgroup.github.io/oxdna-viewer/");
    }
</script>

<iframe width="99%" height="500"  src="https://sulcgroup.github.io/oxdna-viewer/" id="oxview-frame-1" onload="handle()"></iframe>


## Citation
If you find oxDNA Analysis Tools and oxView useful, please cite our paper:

Erik Poppleton, Joakim Bohlin, Michael Matthies, Shuchi Sharma, Fei Zhang, Petr Šulc: Design, optimization and analysis of large DNA and RNA nanostructures through interactive visualization, editing and molecular simulation, *Nucleic Acids Research*, Volume 48, Issue 12, Page e72 (2020). DOI: [10.1093/nar/gkaa417](https://doi.org/10.1093/nar/gkaa417)

And additionally for oxView, please also cite:

Joakim Bohlin, Michael Matthies, Jonah Procyk, Erik Poppleton, Aatmik Mallya, Hao Yan, Petr Šulc: Design and simulation of DNA, RNA and hybrid protein–nucleic acid nanostructures with oxView. *Nature Protocols*, (2022). DOI: [10.1038/s41596-022-00688-5](https://doi.org/10.1038/s41596-022-00688-5) 