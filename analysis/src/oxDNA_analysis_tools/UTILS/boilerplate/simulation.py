"""
The Simulation class: a MutableMapping wrapper around an oxDNA input file.

Provides dictionary-like access to simulation parameters, process management
(run/terminate), configuration inspection, trajectory visualization, and
OAT analysis methods (mean, align, deviations).
"""

import oxpy
from os.path import dirname, abspath, exists
from multiprocessing import Process
import matplotlib.pyplot as plt
import numpy as np
import ipywidgets as widgets
from IPython.display import display
from typing import List, Dict
from collections.abc import MutableMapping

from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, write_conf
from oxDNA_analysis_tools.mean import mean
from oxDNA_analysis_tools.deviations import deviations
from oxDNA_analysis_tools.align import align

from .utils import path_decorator, _prun


class Simulation(MutableMapping):
    """
        MutableMapping wrapper around an oxDNA input file.

        Provides dictionary-like access to simulation parameters, process management
        (run/terminate), configuration inspection, trajectory visualization, and
        OAT analysis methods (mean, align, deviations).

        Parameters:
            input_file_path (str) : Path to an existing oxDNA input file
    """

    def __init__(self, input_file_path: str):
        """
            Creates a simulation object from the given input file path

            Parameters:
                input_file_path (str) : Path to the input file
        """
        self._input_file = oxpy.InputFile()
        self._input_file.init_from_filename(input_file_path)
        self.p = None  # the process reference
        self.out_dir = dirname(abspath(input_file_path))
        self.__input_file_path = input_file_path

    # We support full dictionary like behavior on the Simulation class hiding the InputFile object (MutableMapping inheritance)
    def __setitem__(self, key, value):
        """
            Provides a way to modify input file object parameters.

            Parameters:
                key (str) : The parameter name
                value (str) : The parameter value
        """
        with oxpy.Context():
            self._input_file[key] = value

    def __getitem__(self, key):
        """
            Handles the retrieval of input file parameters.

            Parameters:
                key (str) : The parameter name

            Returns:
                (str) : The parameter value
        """
        return self._input_file[key]

    def __delitem__(self, key):
        with oxpy.Context():
            del self._input_file[key]

    def __iter__(self):
        return (key for key in self._input_file.keys())

    def __len__(self):
        return len(self._input_file.keys())
    #************************************************************************************************

    @path_decorator
    def get_init_conf(self):
        """
            Returns the initial configuration of the simulation as a RyeReader object.

            Returns:
                ((TopInfo, TrajInfo), Configuration) : The topology/trajectory handles and the initial configuration
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["conf_file"]))
        return (ti, di), get_confs(ti, di, 0, 1)[0]

    @path_decorator
    def get_last_conf(self):
        """
            Returns the last configuration of the simulation as a RyeReader object.

            Returns:
                ((TopInfo, TrajInfo), Configuration) : The topology/trajectory handles and the last configuration
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["lastconf_file"]))
        return (ti, di), get_confs(ti, di, 0, 1)[0]

    @path_decorator
    def get_conf_count(self):
        """
            Returns the number of configurations in the trajectory.

            Returns:
                (int) : The number of configurations in the trajectory
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["trajectory_file"]))
        return len(di.idxs)

    @path_decorator
    def get_conf(self, id: int):
        """
            Returns the configuration at the given index in the trajectory
            as a RyeReader object.

            Parameters:
                id (int) : The zero-based index of the desired configuration

            Returns:
                ((TopInfo, TrajInfo), Configuration) : The topology/trajectory handles and the requested configuration
        """
        ti, di = self.get_trajectory_handle()

        l = len(di.idxs)
        if id < l:
            return (ti, di), get_confs(ti, di, id, 1)[0]
        else:
            raise Exception("You requested a conf out of bounds.")

    @path_decorator
    def get_trajectory_handle(self):
        """
            Returns topology and trajectory info handles, useful for processing.

            Returns:
                (TopInfo, TrajInfo) : The topology and trajectory info
        """
        ti, di = describe(abspath(self["topology"]),
                          abspath(self["trajectory_file"]))
        return ti, di

    def view_init(self, overlay: Dict[str, List] = None, script_file_path: str = None, inbox_settings: List[str] = ["Monomer", "Origin"], height: int = 500):
        """
            Opens the initial configuration in an embedded oxDNA viewer window.

            Parameters:
                overlay (Dict[str, List]) : (optional) Dictionary with the color overlay
                script_file_path (str) : (optional) Path to a viewer script file (js)
                inbox_settings (List[str]) : (optional) A list of strings, the inbox settings to use
                height (int) : (optional) Height of the view in pixels. default=500
        """
        (ti, di), conf = self.get_init_conf()
        oxdna_conf(ti, conf, overlay=overlay, script_file_path=script_file_path, inbox_settings=inbox_settings, height=height)

    def view_last(self, overlay: Dict[str, List] = None, script_file_path: str = None, inbox_settings: List[str] = ["Monomer", "Origin"], height: int = 500):
        """
            Opens the last configuration in an embedded oxDNA viewer window.

            Parameters:
                overlay (Dict[str, List]) : (optional) Dictionary with the color overlay
                script_file_path (str) : (optional) Path to a viewer script file (js)
                inbox_settings (List[str]) : (optional) A list of strings, the inbox settings to use
                height (int) : (optional) Height of the view in pixels. default=500
        """
        (ti, di), conf = self.get_last_conf()
        oxdna_conf(ti, conf, inbox_settings=inbox_settings, height=height, overlay=overlay, script_file_path=script_file_path)

    def view_conf(self, id: int, overlay: Dict[str, List] = None, script_file_path: str = None, inbox_settings: List[str] = ["Monomer", "Origin"], height: int = 500):
        """
            Opens the configuration at the given index in the trajectory in an embedded oxDNA viewer window.

            Parameters:
                id (int) : The zero-based index of the configuration to display
                overlay (Dict[str, List]) : (optional) Dictionary with the color overlay
                script_file_path (str) : (optional) Path to a viewer script file (js)
                inbox_settings (List[str]) : (optional) A list of strings, the inbox settings to use
                height (int) : (optional) Height of the view in pixels. default=500
        """
        (ti, di), conf = self.get_conf(id)
        oxdna_conf(ti, conf, inbox_settings=inbox_settings, height=height, overlay=overlay, script_file_path=script_file_path)

    def view_traj(self, init: int = 0, overlay: Dict[str, List] = None, script_file_path: str = None, op=None, inbox_settings: List[str] = ["Monomer", "Origin"], height: int = 500):
        """
            Opens the trajectory in an embedded oxDNA viewer window with an interactive slider.

            Parameters:
                init (int) : (optional) The initial configuration index to display. default=0
                overlay (Dict[str, List]) : (optional) Dictionary with the color overlay
                script_file_path (str) : (optional) Path to a viewer script file (js)
                op (array-like) : (optional) An observable to plot alongside the trajectory
                inbox_settings (List[str]) : (optional) A list of strings, the inbox settings to use
                height (int) : (optional) Height of the view in pixels. default=500
        """
        # get the initial conf and the reference to the trajectory
        (ti, di), cur_conf = self.get_conf(init)

        slider = widgets.IntSlider(
            min=0,
            max=len(di.idxs) - 1,
            step=1,
            description="Select:",
            value=init,
            continuous_update=False,
        )

        output = widgets.Output()
        if op:
            min_v, max_v = np.min(op), np.max(op)

        def handle(obj=None):
            conf = get_confs(ti, di, slider.value, 1)[0]
            with output:
                output.clear_output()
                if op:
                    # make sure our figure is bigger
                    plt.figure(figsize=(15, 3))
                    plt.plot(op)
                    plt.plot([slider.value, slider.value], [min_v, max_v], color="r")
                    plt.show()
                oxdna_conf(ti, conf, inbox_settings=inbox_settings, height=height, overlay=overlay, script_file_path=script_file_path)

        slider.observe(handle)
        display(slider, output)
        handle(None)

    def is_alive(self):
        """
            Prints whether the simulation is still running.
        """
        if self.p:
            print(self.p.is_alive())
        else:
            print(False)

    def terminate(self):
        """
            Terminates a running simulation.
        """
        if self.p:
            self.p.terminate()
            print("you evil...")

    def run(self):
        """
            Runs the simulation as a background process.

            If a previous simulation process exists, it is terminated first.
            The current input file state is written to disk before launching.

            Returns:
                (Simulation) : self, for method chaining
        """

        # if we have any reference to a simulation
        if self.p:
            self.terminate()

        # make sure we pick up new settings if anything was changed
        with open(self.__input_file_path, "w") as file:
            file.write(str(self._input_file))

        # spawn the process
        self.p = Process(target=_prun, args=(self.__input_file_path,))
        self.p.start()
        return self

    @path_decorator
    def compute_mean(self, ref_conf=None, indixes=None, ncpus=1):
        """
            Compute the mean structure from the trajectory and write it to mean.dat.

            Parameters:
                ref_conf (Configuration) : (optional) The reference configuration for alignment
                indixes (List[int]) : (optional) List of nucleotide indices to include
                ncpus (int) : (optional) Number of CPUs for parallel computation. default=1

            Returns:
                ((TopInfo, TrajInfo), Configuration) : The topology/trajectory handles and the mean configuration
        """
        ti, di = self.get_trajectory_handle()
        # Compute the mean structure and RMSFs
        mean_conf = mean(di, ti, ref_conf, indixes, ncpus)
        write_conf("./mean.dat", mean_conf, False, False)
        # now we get its handle
        ti, di = describe(self["topology"], "./mean.dat")
        return (ti, di), mean_conf

    @path_decorator
    def compute_align(self, ref_conf=None, indixes=None, center=True, ncpus=1):
        """
            Align all trajectory frames and write the aligned trajectory to aligned.dat.

            Parameters:
                ref_conf (Configuration) : (optional) The reference configuration for alignment
                indixes (List[int]) : (optional) List of nucleotide indices to align on
                center (bool) : (optional) Whether to center the aligned structures. default=True
                ncpus (int) : (optional) Number of CPUs for parallel computation. default=1

            Returns:
                ((TopInfo, TrajInfo), str) : The topology/trajectory handles and the path to the aligned trajectory
        """
        align(self["trajectory_file"], "./aligned.dat", ncpus, indixes, ref_conf, center)
        # handle to the aligned traj
        ti, di = describe(self["topology"], "./aligned.dat")
        return (ti, di), "./aligned.dat"

    @path_decorator
    def compute_deviations(self, ref_conf=None, indixes=None, ncpus=1):
        """
            Compute RMSD and per-nucleotide RMSF from the trajectory.

            Parameters:
                ref_conf (Configuration) : (optional) The reference configuration
                indixes (List[int]) : (optional) List of nucleotide indices to include
                ncpus (int) : (optional) Number of CPUs for parallel computation. default=1

            Returns:
                (np.ndarray, dict[str, list]) :
                | Root mean squared deviation for each configuration in the trajectory
                | Dictionary with "RMSF" key containing per-nucleotide RMSF values
        """
        ti, di = self.get_trajectory_handle()
        RMSDs, RMSFs = deviations(di, ti, ref_conf, indixes, ncpus)
        # They come out as numpy arrays, need to be a dict with a list for visualization
        RMSFs = {"RMSF": RMSFs.tolist()}
        return RMSDs, RMSFs

    @path_decorator
    def plot_energy(self, ylim=[-2, 0]):
        """
            Plots the energy graph of the running simulation.

            Parameters:
                ylim (list[float]) : (optional) Y-axis limits for the plot as [min, max]. default=[-2, 0]
        """
        if not exists(self["energy_file"]):
            print("No energy file yet.")
            return

        data = np.loadtxt(self["energy_file"])

        # make sure our figure is bigger and make it pretty
        plt.figure(figsize=(15, 3))
        plt.ylabel("Energy")
        plt.xlabel("Steps")
        plt.ylim(ylim)

        if self._input_file["sim_type"] == "MD":
            dt = float(self["dt"])
            steps = float(self["steps"])
            data = np.loadtxt(self["energy_file"])
            df = {'time': data[:, 0],
                  'U': data[:, 1],
            }

            # plot the energy
            plt.plot(df["time"] / dt,
                     df["U"])
        else:
            # assume we have MC
            df = {
                'time': data[:, 0],
                'U': data[:, 1],
            }
            steps = float(self["steps"])
            # plot the energy
            plt.plot(df["time"], df["U"])

        # and the line indicating the complete run
        plt.plot([steps, steps], ylim, color="r")

    @path_decorator
    def get_log(self):
        """
            Read and return the simulation log file contents.

            Returns:
                (list[str]) : Lines from the log file
        """
        with open(self["log_file"]) as file:
            return file.readlines()
