#!/usr/bin/env python3

import os
import toml
import numpy as np
import oxpy

FORCE_FILE = "ext_meta.txt"
OP_FILE = "op.dat"

def build_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description="Helper script for the metadynamics interface for oxDNA. It creates a folder that can be used to run simulations biased with the final metadynamics bias potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--config", type=str, required=True, help="Path to the TOML configuration file printed by the metadynamics interface.")
    parser.add_argument("--force", "-f", action="store_true", help="Whether to overwrite the sampling folder if it already exists.")

    return parser

def get_beta_FE(values, op_range=None, weights=None, bins=100):
    P, bins = np.histogram(values, range=op_range, bins=bins, weights=weights)
    bins = bins[:-1]
    nonzero_mask = P > 0
    beta_FE = -np.log(P[nonzero_mask])
    beta_FE -= np.min(beta_FE)
    return bins[nonzero_mask], beta_FE


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    with open(args.config, "r") as f:
        config_data = toml.load(f)

    # Read the bias potential from the force file
    with open(FORCE_FILE, "r") as f:
        op_min = None
        op_max = None
        N_grid = None
        bias = []
        for line in f.readlines():
            line = line.strip()
            if line == "}":
                break
            spl = [x.strip() for x in line.split("=")]
            if len(spl) == 2:
                if spl[0] == "potential_grid":
                    spl[1] = spl[1].strip(",")
                    bias = np.array([float(x) for x in spl[1].split(",")])
                elif spl[0] == "coord_min" or spl[0] == "xmin":
                    op_min = float(spl[1])
                elif spl[0] == "coord_max" or spl[0] == "xmax":
                    op_max = float(spl[1])
                elif spl[0] == "N_grid":
                    N_grid = int(spl[1])

        if op_min is None or op_max is None or N_grid is None:
            raise ValueError("Could not find potential_grid, coord_min, coord_max, or N_grid in the force file.")

        dx = (op_max - op_min) / (N_grid - 1)
        op_grid = np.arange(op_min, op_max + dx, dx)
        np.savetxt("bias_from_force_file.dat", np.column_stack((op_grid, bias)))

    # Read the order parameter values from the OP file
    op_data = np.loadtxt(OP_FILE)
    biased_values = op_data[:, 0]
    bias_values = op_data[:, 1]

    op_range = (op_min, op_max)
    biased_beta_FE = get_beta_FE(biased_values, op_range=op_range)
    np.savetxt("biased_beta_FE.dat", np.column_stack(biased_beta_FE))

    # Get the temperature from the input file
    with oxpy.Context(print_coda=False):
        input_file = oxpy.InputFile()
        input_file.init_from_filename("input")
        k_BT = oxpy.get_temperature(input_file["T"])

    weights = np.exp(bias_values / k_BT)

    unbiased_beta_FE = get_beta_FE(biased_values, op_range=op_range, weights=weights)
    np.savetxt("unbiased_beta_FE.dat", np.column_stack(unbiased_beta_FE))

    # If the OP is a coordination OP, we can also compute the unbiased distribution of the coordination number 
    # defined as the proper number of hydrogen bonds, which is stored in the third column of the OP file
    if config_data["coordination"] == True:
        hb_biased_values = op_data[:, 2]
        bins = list(range(int(op_max - op_min) + 1))
        hb_unbiased_beta_FE = get_beta_FE(hb_biased_values, weights=weights, bins=bins)
        np.savetxt("unbiased_hb_beta_FE.dat", np.column_stack(hb_unbiased_beta_FE))
