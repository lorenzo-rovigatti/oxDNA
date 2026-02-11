"""
Simulation setup utilities.

Contains:
    - dump_json: Write a dictionary to a JSON file
    - setup_simulation: Prepare a complete simulation directory with all
      required input files (topology, configuration, input, forces)
"""

import oxpy
from os.path import join, exists, basename
from os import mkdir
from shutil import copyfile, rmtree
from json import dumps


def dump_json(obj: dict[str, str], path: str):
    """
        Write a dictionary to a file as JSON.

        Parameters:
            obj (dict[str, str]) : The dictionary to serialize
            path (str) : The output file path
    """
    with open(path, "w") as file:
        file.write(dumps(obj))


def setup_simulation(
    top_path: str,
    dat_path: str,
    out_dir: str,
    parameters: dict[str, str],
    force_dict: dict[str, str] | None = None,
    kill_out_dir: bool = False,
):
    """
        Set up an oxDNA simulation in the given output directory.

        Creates the output directory, copies topology and configuration files
        into it, writes the oxDNA input file from the provided parameters,
        and optionally writes an external forces JSON file.

        Parameters:
            top_path (str) : Path to the topology file
            dat_path (str) : Path to the initial configuration (.dat) file
            out_dir (str) : Path to the output directory to create
            parameters (dict[str, str]) : Dictionary of oxDNA input file parameters
            force_dict (dict[str, str]) : (optional) Dictionary of external forces. If provided,
                a forces.json file is written and the input file is configured to use it.
            kill_out_dir (bool) : (optional) If True, delete and recreate out_dir if it already
                exists. default=False

        Returns:
            (str) : Path to the generated input file
    """
    if exists(out_dir) and not kill_out_dir:
        raise Exception("the output dir is already present, use kill_out_dir to override")

    # Construct the output folder, removing it first if requested
    if kill_out_dir and out_dir != "./":
        if exists(out_dir):
            rmtree(out_dir)
    mkdir(out_dir)

    # Set up the input file object
    input_file = oxpy.InputFile()

    # Copy topology if not already present
    out_top = join(out_dir, basename(top_path))
    if not exists(out_top):
        copyfile(top_path, out_top)
    input_file["topology"] = basename(top_path)

    # Copy configuration if not already present
    out_dat = join(out_dir, basename(dat_path))
    last_conf_path = join(out_dir, "last_conf.dat")

    if not exists(out_dat):
        copyfile(dat_path, out_dat)
        copyfile(dat_path, last_conf_path)

    input_file["conf_file"] = basename(dat_path)

    # Set default output file names
    input_file["trajectory_file"] = "trajectory.dat"
    input_file["energy_file"] = "energy.dat"
    input_file["lastconf_file"] = "last_conf.dat"

    # Write external forces file if provided
    if force_dict:
        force_path = join(out_dir, "forces.json")
        dump_json(force_dict, force_path)
        input_file["external_forces_as_JSON"] = "true"
        input_file["external_forces"] = "1"
        input_file["external_forces_file"] = "forces.json"

    # Apply all user-supplied parameters (overrides defaults above)
    with oxpy.Context():
        for k, v in parameters.items():
            input_file[k] = str(v)

    # Redirect oxDNA stdout to a log file
    input_file["log_file"] = "logfile"

    # Write the input file
    input_file_path = join(out_dir, "input")
    with open(input_file_path, "w") as file:
        file.write(str(input_file))
    return input_file_path
