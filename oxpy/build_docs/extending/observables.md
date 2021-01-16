# Writing observables in Python

You can write custom observables that can analyse configurations on the fly (that is, while the simulation is running) by subclassing `oxpy.BaseObservable` and overloading its [`get_output_string`](../modules/core.html#oxpy.core.BaseObservable.get_output_string) method, which takes a single parameter (the current simulation step) and returns the string to be output. Use [`oxpy.BaseObservable`s `config_info`](../modules/core.html#oxpy.core.BaseObservable.config_info) to access the simulation data required for the analysis (particle positions, velocities, *etc.*).

Use [`OxpyManager`'s `add_output`](../modules/core.html#oxpy.core.OxpyManager.add_output) method to add a new output file and associate observables to it. Use its `print_every` parameter to set the output frequency (in number of simulation time steps).

Here is an example of a custom observable that prints the position of the last particle to the `my_obs_output_file` every 100 time steps: 

	import oxpy
	
	class MyObs(oxpy.BaseObservable):
	    def get_output_string(self, curr_step):
	        # take the position of the last particle
	        pos = self.config_info.particles()[-1].pos
	        # use it to build the output string
	        return "%lf %lf %lf" % (pos[0], pos[1], pos[2])
	
	with oxpy.Context():
	    manager = oxpy.OxpyManager("input")
	    manager.load_options()
	    manager.init()
	
	    my_obs = MyObs()
	    manager.add_output("my_obs_output.dat", print_every=100, observables=[my_obs, ])
	
	    # run 1k steps
	    manager.run(1000)