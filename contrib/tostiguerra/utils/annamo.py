#!/usr/bin/env python3
"""ANNaMo simulation preparation and run tool."""

import argparse
import json
import os
import random
import shutil
import subprocess
import sys

# Ensure sibling modules (generate_annamo, DNA_SL, …) are importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_annamo


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _find_executable(name):
    """Return the path to *name*, searching PATH then ../../build/bin/."""
    found = shutil.which(name)
    if found:
        return found
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.normpath(
        os.path.join(script_dir, "..", "..", "build", "bin", name)
    )
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return None


def _plugin_search_path():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(script_dir, ".."))


# ---------------------------------------------------------------------------
# input.dat generation
# ---------------------------------------------------------------------------


def _write_input_dat(cfg):
    """Write input.dat to the current directory from the JSON config *cfg*."""
    swap = cfg.get("swap", True)
    lambda_val = 1 if swap else 10
    seed = cfg.get("seed", random.randint(0, 2**31 - 1))
    steps = int(cfg.get("steps", 2_000_000_000))
    temperature = cfg.get("temperature", 30)
    print_conf_interval = int(cfg.get("print_conf_interval", 1_000_000))
    print_energy_every = int(cfg.get("print_energy_every", 1_000))
    overrides = cfg.get("oxdna_overrides", {})

    lines = [
        "backend = CPU",
        "sim_type = MD",
        "thermostat = brownian",
        "newtonian_steps = 53",
        "diff_coeff = 0.5",
        f"steps = {steps}",
        "T = 1.0",
        "dt = 0.002",
        "verlet_skin = 0.2",
        "generate_consider_bonded_interactions = true",
        "energy_threshold = 30",
        "reset_com_momentum = true",
        "reset_initial_com_momentum = true",
        "back_in_box = true",
        "restart_step_counter = true",
        "time_scale = linear",
        "refresh_vel = true",
        "",
        "interaction_type = ANNaMoInteraction",
        f"plugin_search_path = {_plugin_search_path()}",
        "",
        "ANNAMO_annamo_version = 1",
        f"ANNAMO_tC = {temperature}",
        "ANNAMO_interaction_matrix_file = dHdS_matrix.dat",
        f"DPS_3b_lambda = {lambda_val}",
        "",
        "topology = topology.top",
        "conf_file = init_conf.dat",
        "trajectory_file = trajectory.dat",
        "energy_file = energy.dat",
        f"print_conf_interval = {print_conf_interval}",
        f"print_energy_every = {print_energy_every}",
        f"seed = {seed}",
        "",
        "data_output_1 = {",
        "name = bonds.dat",
        f"print_every = {print_conf_interval}",
        "col_1 = {",
        "type = pair_energy",
        "}",
        "}",
        "",
        "data_output_2 = {",
        "name = last_backup.dat",
        "only_last = true",
        f"print_every = {print_conf_interval}",
        "col_1 = {",
        "type = configuration",
        "}",
        "}",
    ]

    if overrides:
        lines.append("")
        lines.append("# oxdna_overrides")
        for k, v in overrides.items():
            lines.append(f"{k} = {v}")

    with open("input.dat", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Input file written to: input.dat")


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------


def _validate_config(cfg, json_path):
    errors = []

    material = cfg.get("material")
    if material is None:
        errors.append('missing required field: "material" ("DNA" or "RNA")')
    elif material not in ("DNA", "RNA"):
        errors.append(f'"material" must be "DNA" or "RNA", got: {material!r}')

    strands = cfg.get("strands")
    if strands is None:
        errors.append('missing required field: "strands"')
    elif not isinstance(strands, list) or len(strands) == 0:
        errors.append('"strands" must be a non-empty list')
    else:
        for i, s in enumerate(strands):
            if not isinstance(s, list) or len(s) == 0:
                errors.append(f'"strands[{i}]" must be a non-empty list of beads')

    # salt_concentration is validated but not yet used (salt correction not implemented)
    for key in (
        "temperature",
        "salt_concentration",
        "box_size",
        "print_conf_interval",
        "print_energy_every",
        "steps",
    ):
        val = cfg.get(key)
        if val is not None and not isinstance(val, (int, float)):
            errors.append(f'"{key}" must be a number, got: {val!r}')

    if cfg.get("swap") is not None and not isinstance(cfg["swap"], bool):
        errors.append(f'"swap" must be true or false, got: {cfg["swap"]!r}')

    if cfg.get("oxdna_overrides") is not None and not isinstance(
        cfg["oxdna_overrides"], dict
    ):
        errors.append('"oxdna_overrides" must be a JSON object')

    if errors:
        sys.exit(
            f"Error: invalid configuration in {json_path}:\n"
            + "\n".join(f"  - {e}" for e in errors)
        )


# ---------------------------------------------------------------------------
# prepare subcommand
# ---------------------------------------------------------------------------


def cmd_prepare(args):
    json_path = os.path.abspath(args.system_json)
    if not os.path.isfile(json_path):
        sys.exit(f"Error: file not found: {args.system_json}")

    with open(json_path) as f:
        try:
            cfg = json.load(f)
        except json.JSONDecodeError as e:
            sys.exit(f"Error: could not parse {args.system_json}: {e}")

    _validate_config(cfg, json_path)

    # topology.top and dHdS_matrix.dat (written to cwd by generate_annamo)
    generate_annamo.main(json_path)

    # input.dat (written to cwd)
    _write_input_dat(cfg)

    # init_conf.dat via confGenerator
    conf_gen = _find_executable("confGenerator")
    if conf_gen is None:
        sys.exit(
            "Error: 'confGenerator' not found.\n"
            "Add the oxDNA build/bin directory to your PATH, e.g.:\n"
            "  export PATH=/path/to/oxDNA/build/bin:$PATH"
        )
    box_size = cfg.get("box_size", 30)
    subprocess.run([conf_gen, "input.dat", str(box_size)], check=True)
    print("Initial configuration written to: init_conf.dat")

    print("\nDone. You can now edit input.dat and run:\n  oxDNA input.dat")


# ---------------------------------------------------------------------------
# run subcommand
# ---------------------------------------------------------------------------


def cmd_run(args):
    cmd_prepare(args)

    oxdna = _find_executable("oxDNA")
    if oxdna is None:
        sys.exit(
            "Error: 'oxDNA' not found.\n"
            "Add the oxDNA build/bin directory to your PATH, e.g.:\n"
            "  export PATH=/path/to/oxDNA/build/bin:$PATH"
        )

    print(f"\nLaunching: {oxdna} input.dat")
    subprocess.run([oxdna, "input.dat"], check=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_JSON_EXAMPLE = """\
System JSON fields (required fields marked with *):

  {
    "material":             "DNA",   (* "DNA" or "RNA")
    "strands":              [...],   (* list of strands, each a list of bead sequences)
    "temperature":          30,      (degrees C, default: 30)
    "steps":                2e9,     (default: 2000000000)
    "box_size":             30,      (internal units, default: 30)
    "swap":                 true,    (true → λ=1 bond-swapping; false → λ=10, default: true)
    "seed":                 104123,  (RNG seed, default: random)
    "print_conf_interval":  1e6,     (default: 1000000)
    "print_energy_every":   1e3,     (default: 1000)
    "oxdna_overrides":      {}       (raw key=value pairs appended to input.dat)
  }

Bead design: each bead represents a nucleotide sequence of ideally 3 nt (range 2-4).
Strand division should follow native contacts as described in the ANNaMo paper
(Tosti Guerra et al., J. Chem. Phys. 2024): beads are chosen so that native base pairs
fall at bead boundaries rather than within a single bead.

Example:
  "strands": [
    ["GAA", "GTG", "ACA", "TGG"],
    ["CCA", "TGT", "CAC", "TTC"]
  ]
"""


def main():
    parser = argparse.ArgumentParser(
        prog="annamo",
        description="ANNaMo simulation preparation and run tool.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_JSON_EXAMPLE,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_prepare = sub.add_parser(
        "prepare",
        help="Generate all input files without launching the simulation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_JSON_EXAMPLE,
    )
    p_prepare.add_argument("system_json", help="Path to the system JSON file.")
    p_prepare.set_defaults(func=cmd_prepare)

    p_run = sub.add_parser(
        "run",
        help="Generate all input files and launch the oxDNA simulation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_JSON_EXAMPLE,
    )
    p_run.add_argument("system_json", help="Path to the system JSON file.")
    p_run.set_defaults(func=cmd_run)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
