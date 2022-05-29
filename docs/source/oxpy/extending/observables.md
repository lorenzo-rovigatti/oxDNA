# Writing observables in Python

You can write custom observables that can analyse configurations on the fly (that is, while the simulation is running) by subclassing {class}`~oxpy.core.BaseObservable` and overloading its {meth}`~oxpy.core.BaseObservable.get_output_string` method, which takes a single parameter (the current simulation step) and returns the string to be output. Use the {attr}`~oxpy.core.BaseObservable.config_info` attribute to access the simulation data required for the analysis (particle positions, velocities, *etc.*).

Use {class}`~oxpy.core.OxpyManager`'s {meth}`~oxpy.core.OxpyManager.add_output` method to add a new output file and associate observables to it. Use its `print_every` parameter to set the output frequency (in number of simulation time steps).

Here is an example of a custom observable that prints the position of the last particle to `my_obs_output_file.dat` every 100 time steps: 

```python
import oxpy

class MyObs(oxpy.observables.BaseObservable):
    def get_output_string(self, curr_step):
        # take the position of the last particle
        pos = self.config_info.particles()[-1].pos
        # use it to build the output string
        return "%lf %lf %lf" % (pos[0], pos[1], pos[2])

with oxpy.Context():
    manager = oxpy.OxpyManager("input")

    my_obs = MyObs()
    manager.add_output("my_obs_output.dat", print_every=100, observables=[my_obs, ])

    # run 1k steps
    manager.run(1000)
```
