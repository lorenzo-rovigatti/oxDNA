#!/usr/bin/env python3

import os, shutil
import toml
import oxpy
import random

BASE_RUNNER_DIR = "run-meta_0"
RUNNER_INPUT_FILE = "input-meta"
METAD_CONFIG_FILENAME = "original_metad_config.toml"

FILES = {
    "last_conf.dat" : "init_conf.dat",
    "ext_meta.txt" : "ext_meta.txt"
}

LINES = {
    "conf_file" : ["last_conf.dat", "init_conf.dat"],
    "log_file" : ["oxDNA_log.txt", ""],
    "no_stdout_energy" : ["true", ""],
    "print_conf_interval" : ["1e11", "1e7"],
    "restart_step_counter" : ["", "true"],
    "seed" : ["*", ""]
}

def build_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description="Helper script for the metadynamics interface for oxDNA. It creates a folder that can be used to run simulations biased with the final metadynamics bias potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("sampling_dir", help="The path to the folder where the sampling files will be created. This folder will be created if it does not exist, and, if --force is given, overwritten if it already exists.")

    parser.add_argument("--config", type=str, required=True, help="Path to the TOML configuration file printed by the metadynamics interface.")
    parser.add_argument("--force", "-f", action="store_true", help="Whether to overwrite the sampling folder if it already exists.")
    parser.add_argument("--seed", "-s", type=int, default=random.randint(0, 2**32 - 1), help="The seed to use for the oxDNA simulation. If not given, the seed will be chosen randomly.")
    parser.add_argument("--rescale_bias", action="store_true", help="Multiply the final bias by (1 + dT) / dT, so that the sampling simulation is biased by (minus) the free-energy estimate. Without this option the well-tempered bias compensates only a fraction dT / (1 + dT) of the free energy, leaving a residual landscape F / (1 + dT).")

    return parser


def rescale_potential_grids(filename: str, factor: float) -> int:
    """Multiply all the values of the potential grids contained in the given external forces file by factor.
    Handles both 1D grids (comma-separated values) and 2D grids (rows separated by '|'). Returns the number
    of grids that have been rescaled."""
    with open(filename) as f:
        lines = f.readlines()

    N_grids = 0
    for i, line in enumerate(lines):
        key, sep, value = line.partition("=")
        if sep and key.strip() == "potential_grid":
            # empty tokens (e.g. due to trailing commas or '|') are kept as they are
            rows = [",".join(f"{float(v) * factor}" if v.strip() else v for v in row.split(",")) for row in value.strip().split("|")]
            indent = key[:len(key) - len(key.lstrip())]
            lines[i] = f"{indent}potential_grid = {'|'.join(rows)}\n"
            N_grids += 1

    with open(filename, "w") as f:
        f.writelines(lines)

    return N_grids


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    with open(args.config, "r") as f:
        config_data = toml.load(f)

    if args.rescale_bias:
        if "dT" not in config_data:
            print(f"The configuration file {args.config} does not contain the dT key, which is required by --rescale_bias. Use the metad_PID.toml file printed by the metadynamics interface.")
            exit(1)
        bias_factor = (1 + config_data["dT"]) / config_data["dT"]

    # Create the sampling folder
    if os.path.exists(args.sampling_dir):
        if not args.force:
            print(f"Sampling folder {args.sampling_dir} already exists. Use --force to overwrite it.")
            exit(1)
        shutil.rmtree(args.sampling_dir)
    os.makedirs(args.sampling_dir)

    # If the OP is a coordination OP, copy the order-parameter file as well
    if config_data["coordination"] == True:
        FILES["op_coordination.dat"] = "op_coordination.dat"

    # Copy also the metadynamics configuration file
    shutil.copy(args.config, os.path.join(args.sampling_dir, METAD_CONFIG_FILENAME))

    with oxpy.Context(print_coda=False):
        input_file = oxpy.InputFile()
        input_file.init_from_filename(os.path.join(BASE_RUNNER_DIR, RUNNER_INPUT_FILE))
        # Read the name of the topology file from the input file and add it to FILES so that it is copied over to the sampling dir
        topology_file = input_file["topology"]
        FILES[topology_file] = topology_file

        # Replace the lines of the input file that need to be changed to run the sampling
        LINES["seed"][1] = f"{args.seed}"
        for key in LINES:
            if key in input_file:
                option = input_file[key]
                if LINES[key][0] != "" and LINES[key][0] != "*" and option != LINES[key][0]:
                    print(f"Warning: the value of {key} in the input file is not {LINES[key][0]} as expected, but {option}. The line will be replaced anyway, but make sure that we are working with the correct input file.")

            if LINES[key][1] == "" and key in input_file:
                del input_file[key]
            else:
                input_file[key] = LINES[key][1]

        # Print the modified input file to the sampling directory
        input_filename = os.path.join(args.sampling_dir, "input")
        with open(input_filename, "w") as f:
            print(input_file, file=f)

    # Copy the base simulation files
    for f_old, f_new in FILES.items():
        src = os.path.join(BASE_RUNNER_DIR, f_old)
        dst = os.path.join(args.sampling_dir, f_new)
        shutil.copy(src, dst)

    if args.rescale_bias:
        N_grids = rescale_potential_grids(os.path.join(args.sampling_dir, FILES["ext_meta.txt"]), bias_factor)
        print(f"The bias has been multiplied by (1 + dT) / dT = {bias_factor:.6f} (dT = {config_data['dT']}, {N_grids} potential grids rescaled)")
