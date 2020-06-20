from .core import InputFile

def generate_default_input(options=[]):
    default_input = InputFile()
    
    default_input["backend"] = "CPU"
    default_input["sim_type"] = "MD"
    
    default_input["verlet_skin"] = "0.2"
    default_input["dt"] = "0.001"
    
    default_input["T"] = "0.1"

    default_input["steps"] = "10000"
    default_input["print_energy_every"] = "1000"
    default_input["print_conf_interval"] = "100000"
    default_input["restart_step_counter"] = "yes"
    default_input["refresh_vel"] = "true"
    default_input["time_scale"] = "linear"
    
    default_input["topology"] = "topology.top"
    default_input["conf_file"] = "init_conf.dat"
    default_input["trajectory_file"] = "trajectory.dat"
    default_input["energy_file"] = "energy.dat"
    
    return default_input
