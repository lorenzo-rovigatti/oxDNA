# Oxpy

Oxpy is a Python3 library that makes it possible to use oxDNA from Python.

## A simple example

The following snippet imports the `oxpy` module, initialises the simulation machinery, runs a short simulation and computes the average position of the final configuration:

	import numpy as np
	import oxpy
	
	with oxpy.Context():
	    # init the manager with the given input file
	    manager = oxpy.OxpyManager(["input"])
	    manager.load_options()
	    manager.init()
	
	    # run 1k steps
	    manager.run(1000)
	
	    # run 10k steps more
	    manager.run(10000)
	
	    # do some computation with the current configuration
	    particles = manager.config_info().particles()
	    avg_pos = np.average(list(map(lambda p: p.pos, particles)), axis=0)
	    print("Average final position:", avg_pos)
	
## Library API

```eval_rst
.. toctree::
   :maxdepth: 2
   
   core/core.md
```

## Extending Oxpy

