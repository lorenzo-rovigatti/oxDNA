import numpy as np
import oxpy

class MyObs(oxpy.observables.BaseObservable):
    def get_output_string(self, curr_step):
        # take the position of the last particle
        pos = self.config_info.particles()[-1].pos
        # use it to build the output string
        return "%lf %lf %lf" % (pos[0], pos[1], pos[2])

with oxpy.Context():
    my_input = oxpy.InputFile()
    my_input.init_from_filename("input")
    
    # we overwrite some options: note that the original input's conf_file does 
    # not exist and the simulation wouldn't start if this was not done
    my_input["conf_file"] = "init.dat"
    my_input["seed"] = "123456"
    
    manager = oxpy.OxpyManager(my_input)
    
    # we add a custom observable to the simulation
    my_obs = MyObs()
    manager.add_output("last_nucleotide.dat", print_every=100, observables=[my_obs, ])

    # run 1k steps
    manager.run(1000)
    
with oxpy.Context():
    inp = oxpy.InputFile()
    inp.init_from_filename("input_analysis")
    # this object will make it possible to access the trajectory data
    backend = oxpy.analysis.AnalysisBackend(inp)

    with open("first_nucleotide.dat", "w") as f:
        # loop over all the configurations stored in the trajectory file
        while backend.read_next_configuration():
            # print the position of the first nucleotide
            print(backend.config_info().particles()[0].pos, file=f)
