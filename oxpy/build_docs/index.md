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
	    
	    # compute the average position of the particles' backbones
	    avg_pos = np.average(list(map(lambda p: p.backbone_site(), particles)), axis=0)
	    print("Average final position:", avg_pos)
	    
	    # and the interaction energy between the first two particles
	    print("Interaction energy between particle 0 and particle 1:", manager.config_info().interaction.pair_interaction(particles[0], particles[1]))
	
## Library API

```eval_rst
.. toctree::
   :maxdepth: 2
   
   core/core.md
```

## Exceptions

The oxDNA code raises `oxDNAException`s when the simulation cannot be correctly initialised or when it incurs in an unrecoverable error. These exceptions are automatically translated into Python exceptions of type `oxpy.core.OxDNAError`, which can then be handled in a regular [`try except` block](https://docs.python.org/3/tutorial/errors.html).  

## Extending Oxpy

To be written.
