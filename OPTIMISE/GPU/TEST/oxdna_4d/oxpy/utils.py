from .core import InputFile


def generate_default_input(options=[]):
    '''Build and return a minimal input file
    
    The returned input file contains the minimum amount of options required to start a simulation. 
    You can list them by printing the input file after generation::
    
        my_input = oxpy.utils.generate_default_input()
        print(my_input)
        
    Returns
    -------
    :class:`~oxpy.core.InputFile`
        An input file object that can be used to set up an oxDNA simulation. 
    
    '''
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


def Kelvin_to_oxDNA(T):
    '''Convert the temperature given in Kelvin to oxDNA simulation units
    
    Parameters
    ----------
        T : float
            The temperature (in degrees Kelvin) to be converted
            
    Returns
    -------
    float
        The temperature in oxDNA units
    '''
    return float(T) / 3000.0


def Celsius_to_oxDNA(T):
    '''Convert the temperature given in Celsius to oxDNA simulation units
    
    Parameters
    ----------
        T : float
            The temperature (in degrees Celsius) to be converted
            
    Returns
    -------
    float
        The temperature in oxDNA units
    '''
    return Kelvin_to_oxDNA(float(T) + 273.15)

