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