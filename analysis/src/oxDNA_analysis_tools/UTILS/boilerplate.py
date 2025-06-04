import oxpy
from os.path import join, abspath, dirname, exists, basename
from os import mkdir,chdir, getcwd
from multiprocessing import Process
from shutil import copyfile, rmtree
from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
from json import dumps
import matplotlib.pyplot as plt
from copy import deepcopy
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs,write_conf
import ipywidgets as widgets
from IPython.display import display, IFrame
import numpy as np
from functools import wraps
from typing import List, Union, Dict
import contextlib 
from collections.abc import MutableMapping
#oat analysis
from oxDNA_analysis_tools.mean import mean
from oxDNA_analysis_tools.deviations import deviations
from oxDNA_analysis_tools.align import align


default_input_file = {
    "T" :"20C",
    "steps" :"1e9",
    "salt_concentration" :"1",
    "backend" :"CUDA",
    "interaction_type" :"DNA2",
    "print_conf_interval" :"1e5",
    "print_energy_every"  :"1e4",
    "dt" :"0.005",
    "sim_type" :"MD",
    "max_density_multiplier" :"10",
    "verlet_skin" :"0.5",
    "time_scale" :"linear",
    "ensemble" :"NVT",
    "thermostat" :"john",
    "diff_coeff" :"2.5",
    "backend_precision" :"mixed",
    "refresh_vel" :"1",
    "restart_step_counter" :"1",
    "newtonian_steps" :"103",
    "CUDA_list" :"verlet",
    "CUDA_sort_every" :"0",
    "use_edge" :"1",
    "edge_n_forces" :"1",
    "no_stdout_energy" : "true", 
    #"max_density_multiplier":""
}


class SimulationRunningError(Exception):
    pass #used to indicate we have a running simulation


class PathContext(contextlib.ContextDecorator):
    """Class handling the temporary change of path required for the simulations"""
    def __init__(self, out_dir):
        self.out_dir = out_dir
        self.old_path = None

    def __enter__(self):
        self.old_path = getcwd()
        chdir(self.out_dir)

    def __exit__(self, exc_type, exc_value, traceback):
        chdir(self.old_path)

def path_decorator(func):
    """Helper to wrap the PathContext around the Simulation Class Methods"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        with PathContext(self.out_dir):
            return func(self, *args, **kwargs)
    return wrapper


def get_default_input():
    """
        returns a deepcopy of the default input file
    """
    return deepcopy(default_input_file)

def dump_json(obj:dict[str,str], path:str):
    """
       helper function dumps the given dictionary to the given path as a json file
    """
    with open(path, "w+") as file:
        file.write(dumps(obj))
        
def setup_simulation(top_path:str, dat_path:str, out_dir:str, parameters:Dict[str,str], force_dict:Union[Dict[str,str],None]=None, kill_out_dir=False):
    """
        sets up a simulation in the given output directory
        
        top_path (str): path to the topology file
        dat_path (str): path to the initial configuration
        out_dir (str): path to the output directory
        parameters (dict[str,str]): dictionary of parameters to set
        force_dict (dict): dictionary of forces to set (if empty no forces are set)
        kill_out_dir (bool): if true, the output directory will be deleted if present
    """
    if(exists(out_dir) and not kill_out_dir):
        raise Exception("the output dir is already present, use kill_out_dir to override")
    
    # construct ourseves an output folder
    if(kill_out_dir and out_dir != "./"): 
        # we preemtively kill the dir if it is present
        if(exists(out_dir)):
            rmtree(out_dir)
    mkdir(out_dir)
    #set up input file
    input_file = oxpy.InputFile()
    
    # copy the topology over if not present
    out_top = join(out_dir,basename(top_path))
    if(not exists(out_top)):
        copyfile(top_path, out_top)
    input_file["topology"] = basename(top_path)
    
    #copy the conf over if not present
    out_dat = join(out_dir, basename(dat_path))
    last_conf_path = join(out_dir,"last_conf.dat")
    
    if(not exists(out_dat)):
        copyfile(dat_path, out_dat)
        copyfile(dat_path, last_conf_path)
        
    input_file["conf_file"] = basename(dat_path)
        
    #set the outputs to desired location
    input_file["trajectory_file"] = "trajectory.dat"
    input_file["energy_file"] = "energy.dat"
    input_file["lastconf_file"] ="last_conf.dat"
    
    # if we have defined forces we need to 
    # 1) build the file
    # 2) add it to the input_file
    if force_dict:
        force_path = join(out_dir, "forces.json")
        dump_json(force_dict, force_path)
        input_file["external_forces_as_JSON"] = "true"
        input_file["external_forces"] = "1"
        input_file["external_forces_file"] = "forces.json"
    
    
    # stuff wich was already set triggers logging, 
    # so we need to provide it 
    with oxpy.Context():
        # set all base simulation parameters, overriding previous
        for k,v in parameters.items():
            input_file[k]=str(v)
        
    
    #prevents output
    input_file["log_file"] = "logfile"
    
    input_file_path = join(out_dir,"input")
    with open(input_file_path,"w") as file:
        file.write(str(input_file))
    return input_file_path

# our typical run function for a provided dictionary input file 
def _prun(input_file_path:str):
    with oxpy.Context():
        path = dirname(abspath(input_file_path))
        # change to the directory of the input file path
        with PathContext(path):
            input_file = oxpy.InputFile()
            input_file.init_from_filename(basename(input_file_path))
            manager = oxpy.OxpyManager(input_file)
            manager.run_complete() #run complete run's it till the number steps specified are reached 
            

class Simulation(MutableMapping):
    def __init__(self, input_file_path:str):
        """
            creates a simulation object from the given input file path

            input_file_path (str): path to the input file
        """
        self._input_file = oxpy.InputFile()
        self._input_file.init_from_filename(input_file_path)
        self.p = None # the process refference
        self.out_dir = dirname(abspath(input_file_path))
        self.__input_file_path = input_file_path


    # We support full dictionary like behavior on the Simulation class hiding the InputFile object (MutableMapping inheretance)
    def __setitem__(self, key, value):
        """
            provides a way to modify input file object parameters,
            s = Simulation("./input")
            s["conf_file"] = "./new_name" 
        """
        with oxpy.Context():
            self._input_file[key] = value

    def __getitem__(self, key):
        """ handles the retrival of input file parameters """
        return self._input_file[key]
    
    def __delitem__(self, key):
        with oxpy.Context():
            del self._data[key]

    def __iter__(self):
        return (key for key in self._input_file.keys()) 

    def __len__(self):
        return len(self._input_file.keys())
    #************************************************************************************************
     
    @path_decorator
    def get_init_conf(self):
        """
            returns the initial configuration of the simulation as a rye reader object
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["conf_file"]))
        return (ti, di), get_confs(ti, di, 0, 1)[0]
    
    @path_decorator
    def get_last_conf(self):
        """
            returns the last configuration of the simulation as a rye reader object
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["lastconf_file"]))
        return (ti,di), get_confs(ti, di, 0,1)[0]
                             
    @path_decorator
    def get_conf_count(self):
        """
            returns the number of configurations in the trajectory
        """
        ti,di = describe(abspath(self["topology"]),
                         abspath(self["trajectory_file"]))
        return len(di.idxs)
    
    @path_decorator
    def get_conf(self, id:int):
        """
            returns the configuration at the given index in the trajectory
            as a rye reader object
        """
        ti,di = self.get_trajectory_handle()

        l = len(di.idxs)
        if(id < l):
            return (ti,di), get_confs(ti,di, id, 1)[0]
        else:
            raise Exception("You requested a conf out of bounds.")
    
    @path_decorator
    def get_trajectory_handle(self):
        "usefull for processing stuff"
        ti,di = describe(abspath(self["topology"]),
                         abspath(self["trajectory_file"]))
        return ti,di

    def view_init(self, overlay:Dict[str,List] = None, script_file_path:str = None, inbox_settings:List[str] = ["Monomer", "Origin"], height:int = 500):
        """
            opens the initial configuration in an embeded oxDNA viewer window

            inbox_settings (List[str]): a list of strings, the inbox settings to use
            height (int): height of the view
        """
        (ti,di), conf = self.get_init_conf()        
        oxdna_conf(ti, conf, overlay=overlay, script_file_path = script_file_path, inbox_settings=inbox_settings, height=height)

    def view_last(self, overlay:Dict[str,List] = None, script_file_path:str = None,inbox_settings:List[str] = ["Monomer", "Origin"], height:int = 500):
        """
            opens the last configuration in an embeded oxDNA viewer window

            inbox_settings (List[str]): a list of strings, the inbox settings to use
            height (int): height of the view
        """
        (ti,di), conf = self.get_last_conf()
        oxdna_conf(ti, conf, inbox_settings=inbox_settings, height=height, overlay=overlay,script_file_path = script_file_path)
 
    def view_conf(self, id:int, overlay:Dict[str,List] = None, script_file_path:str = None, inbox_settings:List[str] =  ["Monomer", "Origin"], height:int = 500):
        """ 
            opens the configuration at the given index in the trajectory as an embeded oxDNA viewer window

            inbox_settings (List[str]): a list of strings, the inbox settings to use
            height (int): height of the view
        """
        (ti,di), conf = self.get_conf(id)
        oxdna_conf(ti, conf, inbox_settings=inbox_settings, height=height,overlay=overlay,script_file_path=script_file_path)
    
    def view_traj(self, init:int = 0, overlay:Dict[str,List] = None, script_file_path:str = None, op=None, inbox_settings:List[str] = ["Monomer", "Origin"], height:int = 500):
        """
            opens the trajectory in an embeded oxDNA viewer window

            init (int): the initial configuration to start the trajectory from
            op: an optional observable to plot along side the trajectory
            inbox_settings (List[str]): a list of strings, the inbox settings to use
            height (int): height of the view
        """
        # get the initial conf and the reference to the trajectory 
        (ti,di), cur_conf = self.get_conf(init)
        
        slider = widgets.IntSlider(
            min = 0,
            max = len(di.idxs)-1,
            step=1,
            description="Select:",
            value=init, 
            continuous_update=False,
        )
        
        output = widgets.Output()
        if op:
            min_v,max_v = np.min(op), np.max(op)
        
        def handle(obj=None):
            conf= get_confs(ti,di,slider.value,1)[0]
            with output:
                output.clear_output()
                if op:
                    # make sure our figure is bigger
                    plt.figure(figsize=(15,3)) 
                    plt.plot(op)
                    print(init)
                    plt.plot([slider.value,slider.value],[min_v, max_v], color="r")
                    plt.show()
                oxdna_conf(ti,conf, inbox_settings=inbox_settings, height=height, overlay=overlay, script_file_path=script_file_path)
                
        slider.observe(handle)
        display(slider,output)
        handle(None)       
        
    def is_alive(self):
        """
            returns true if the simulation is still running
        """
        if(self.p):
            print(self.p.is_alive())
        else:
            print(False)

    def terminate(self):
        """
            terminates a running simulation
        """
        if(self.p):
            self.p.terminate()
            print("you evil...")

    def run(self):
        """
            runs the simulation
        """

        # if we have any reference to a simulation
        if(self.p):
            self.terminate()


        #make sure we pick up new settings if anything was changed
        with open(self.__input_file_path,"w") as file:
            file.write(str(self._input_file))
        
        #spawn the process
        self.p = Process(target=_prun, args = (self.__input_file_path,))
        self.p.start()
        return self
    
    @path_decorator 
    def compute_mean(self,ref_conf = None, indixes = None, ncpus=1):
        """returns ti,di of computed mean structure"""
        ti,di = self.get_trajectory_handle()
        # Compute the mean structure and RMSFs
        mean_conf = mean(di, ti, ref_conf, indixes, ncpus)
        write_conf("./mean.dat", mean_conf, False, False)
        #now we get it's handle 
        ti,di = describe(self["topology"], "./mean.dat")
        return (ti,di), mean_conf
    
    @path_decorator
    def compute_align(self, ref_conf = None, indixes = None, center=True, ncpus=1):
        align(self["trajectory_file"], "./aligned.dat", ncpus, indixes, ref_conf, center)
        # handle to the aligned traj
        ti,di = describe(self["topology"], "./aligned.dat")
        return (ti, di), "./aligned.dat"
    
    @path_decorator
    def compute_deviations(self, ref_conf = None,  indixes = None, ncpus=1):
        ti,di = self.get_trajectory_handle()
        RMSDs, RMSFs = deviations(di, ti, ref_conf, indixes,  ncpus)
        # #They come out as numpy arrays, need to be a dict with a list for visualization
        RMSFs = {"RMSF": RMSFs.tolist()}
        return RMSDs, RMSFs

    @path_decorator
    def plot_energy(self, ylim = [-2,0]):
        """
            plots the energy graph of the running simulation
        """
        if not exists(self["energy_file"]):
            print("No energy file yet.")
            return
        
        data = np.loadtxt(self["energy_file"])

        # make sure our figure is bigger and make it pretty
        plt.figure(figsize=(15,3)) 
        plt.ylabel("Energy")
        plt.xlabel("Steps")
        plt.ylim(ylim)
        
        if self._input_file["sim_type"] == "MD":
            dt = float(self["dt"])
            steps = float(self["steps"])
            data = np.loadtxt(self["energy_file"])
            #used to be a df, now we use only numpy
            df = {'time': data[:, 0], 
                  'U': data[:, 1], 
                  #'P': data[:, 2], 
                  #'K': data[:, 3]
            }

            # plot the energy
            plt.plot(df["time"]/dt,
                     df["U"])
        else:
            # assume we have MC
            df = {
                'time': data[:, 0],
                'U': data[:, 1],
                #'accept_translations': data[:, 2],
                #'accept_rotations': data[:, 3],
                #'smth': data[:, 4]
            }
            steps = float(self["steps"])
            # plot the energy
            plt.plot(df["time"],df["U"])

        # and the line indicating the complete run
        plt.plot([steps,steps],ylim, color="r")

    @path_decorator
    def get_log(self):
        """default boilerplate sim output is stored in the log file"""
        with open(self["log_file"]) as file:
            return file.readlines()