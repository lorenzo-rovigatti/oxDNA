# Oxpy

Oxpy is a Python3 library that makes it possible to use oxDNA from Python.

## An example of a simple simulation

The following snippet imports the `oxpy` module, initialises the simulation machinery, runs a short simulation using the input file `input`, changes the temperature, runs more simulations steps and computes the average position of the final configuration:

```python
import numpy as np
import oxpy

with oxpy.Context():
	# init the manager with the given input file
	manager = oxpy.OxpyManager("input")
	
	# run 1k steps
	manager.run(1000)
	
	# change the temperature
	manager.update_temperature(0.11)
	
	# run 1k steps more
	manager.run(1000)
	
	# do some computation with the current configuration
	particles = manager.config_info().particles()
	
	# compute the average position of the particles' backbones
	avg_pos = np.average(list(map(lambda p: p.backbone_site(), particles)), axis=0)
	print("Average final position:", avg_pos)
	
	# and the interaction energy between the first two particles
	print("Interaction energy between particle 0 and particle 1:", manager.config_info().interaction.pair_interaction(particles[0], particles[1]))
```
	    
If you want, you can initialise the input file yourself and change some of the options before initialising the manager:

```python
my_input = oxpy.InputFile()
my_input.init_from_filename("input")
my_input["backend"] = "CUDA"
my_input["steps"] = "1e9"
manager = oxpy.OxpyManager(my_input)
```
	
You can also use {func}`oxpy.utils.generate_default_input()` to generate the following basic input file:

```python
backend = CPU
sim_type = MD

verlet_skin = 0.2
dt = 0.001

T = 0.1

steps = 10000
print_energy_every = 1000
print_conf_interval = 100000
restart_step_counter = yes
refresh_vel = true
time_scale = linear

topology = topology.top
conf_file = init_conf.dat
trajectory_file = trajectory.dat
energy_file = energy.dat
```
	
## An example of a simple analysis

Here we loop over all the configurations stored in an oxDNA trajectory file, printing the position of the first particle. 

```python
import numpy as np
import oxpy

with oxpy.Context():
    inp = oxpy.InputFile()
    inp.init_from_filename("input")
    # this object will make it possible to access the trajectory data
    backend = oxpy.analysis.AnalysisBackend(inp)

    # loop over all the configurations stored in the trajectory file
    while backend.read_next_configuration():
        # you can access the particles' details from BaseParticle instances
        print(backend.particles[0].pos)
        # or from the flattened_conf object, which exposes the simulation data as vectors that can be converted to numpy arrays
        numpy_positions = np.array(backend.flattened_conf.positions, copy=False)
        print(numpy_positions[0])
```
	
## Library API

```{eval-rst}
.. toctree::
   :maxdepth: 2
   
   modules/core/core.md
   modules/analysis.md
   modules/utils.md
```

## Exceptions

The oxDNA code raises `oxDNAException`s when the simulation cannot be correctly initialised or when it incurs in an unrecoverable error. These exceptions are automatically translated into Python exceptions of type {class}`oxpy.core.OxDNAError`, which can then be handled in a regular [`try except` block](https://docs.python.org/3/tutorial/errors.html).  

## Extending Oxpy

```{eval-rst}
.. toctree::
   :maxdepth: 2
   
   extending/observables.md
   extending/forces.md
```
