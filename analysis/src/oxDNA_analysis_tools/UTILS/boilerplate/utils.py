"""
Internal utilities for the boilerplate package.

Contains:
    - SimulationRunningError: Exception for running simulations
    - PathContext: Context manager for temporary directory changes
    - path_decorator: Decorator wrapping methods with PathContext
    - _prun: Process runner for launching oxDNA via multiprocessing
"""

import oxpy
import contextlib
from os.path import dirname, abspath, basename
from os import chdir, getcwd
from functools import wraps


class SimulationRunningError(Exception):
    """
        Exception used to indicate we have a running simulation
    """
    pass


class PathContext(contextlib.ContextDecorator):
    """
        Context manager handling the temporary change of working directory required for the simulations

        Parameters:
            out_dir (str) : The directory to change to
    """
    def __init__(self, out_dir):
        self.out_dir = out_dir
        self.old_path = None

    def __enter__(self):
        self.old_path = getcwd()
        chdir(self.out_dir)

    def __exit__(self, exc_type, exc_value, traceback):
        chdir(self.old_path)


def path_decorator(func):
    """
        Decorator that wraps Simulation class methods with PathContext,
        ensuring the working directory is set to the simulation output
        directory for the duration of the call.

        Parameters:
            func (callable) : The method to wrap

        Returns:
            (callable) : The wrapped method
    """
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        with PathContext(self.out_dir):
            return func(self, *args, **kwargs)
    return wrapper


def _prun(input_file_path: str):
    """
        Run an oxDNA simulation to completion (used as multiprocessing target).

        Parameters:
            input_file_path (str) : Path to the oxDNA input file
    """
    with oxpy.Context():
        path = dirname(abspath(input_file_path))
        with PathContext(path):
            input_file = oxpy.InputFile()
            input_file.init_from_filename(basename(input_file_path))
            manager = oxpy.OxpyManager(input_file)
            manager.run_complete()
