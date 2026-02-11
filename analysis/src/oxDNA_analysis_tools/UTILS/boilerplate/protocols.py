"""
Simulation protocol parameter dictionaries.

Contains the original default_input_file and three new protocol
parameter sets derived from the example input files:

    - MC relaxation (resolve excluded volume clashes)
    - MD relaxation (relax strained structures under external forces)
    - Production MD (equilibrium sampling)

Each getter function returns a deep copy so callers can modify freely.
"""

from copy import deepcopy


# ---------------------------------------------------------------------------
# Legacy default parameters (preserved from the original boilerplate module)
# ---------------------------------------------------------------------------

default_input_file = {
    "T":                      "20C",
    "steps":                  "1e9",
    "salt_concentration":     "1",
    "backend":                "CUDA",
    "interaction_type":       "DNA2",
    "print_conf_interval":    "1e5",
    "print_energy_every":     "1e4",
    "dt":                     "0.005",
    "sim_type":               "MD",
    "max_density_multiplier": "10",
    "verlet_skin":            "0.5",
    "time_scale":             "linear",
    "ensemble":               "NVT",
    "thermostat":             "john",
    "diff_coeff":             "2.5",
    "backend_precision":      "mixed",
    "refresh_vel":            "1",
    "restart_step_counter":   "1",
    "newtonian_steps":        "103",
    "CUDA_list":              "verlet",
    "CUDA_sort_every":        "0",
    "use_edge":               "1",
    "edge_n_forces":          "1",
    "no_stdout_energy":       "true",
}


# ---------------------------------------------------------------------------
# Monte-Carlo relaxation protocol
# ---------------------------------------------------------------------------
# Derived from example_input_files/input_relax_MC
#
# Purpose: Initial coarse relaxation on CPU to resolve excluded-volume
# clashes present in structures freshly exported from design tools
# (cadnano, Tiamat, etc.). Uses large MC moves at slightly elevated
# temperature with backbone force clamping.

_relax_MC_input = {
    # Program parameters
    "sim_type":               "MC",
    "backend":                "CPU",
    "backend_precision":      "double",
    "verlet_skin":            "0.5",
    # Simulation parameters
    "interaction_type":       "DNA2",
    "steps":                  "5e3",
    "dt":                     "0.05",
    "T":                      "30C",
    "salt_concentration":     "1.0",
    "ensemble":               "nvt",
    "delta_translation":      "0.22",
    "delta_rotation":         "0.22",
    "diff_coeff":             "2.5",
    "max_backbone_force":     "5",
    "max_backbone_force_far": "10",
    # I/O parameters
    "lastconf_file":          "last_conf_MC.dat",
    "trajectory_file":        "MC_traj.dat",
    "energy_file":            "energy_MC.dat",
    "print_conf_interval":    "500",      # steps/10 for 5e3
    "print_energy_every":     "100",      # steps/50 for 5e3
    "time_scale":             "linear",
    "refresh_vel":            "0",
    "restart_step_counter":   "1",
    "no_stdout_energy":       "0",
}


# ---------------------------------------------------------------------------
# MD relaxation protocol (with external forces)
# ---------------------------------------------------------------------------
# Derived from example_input_files/input_relax_MD
#
# Purpose: Fine-grained relaxation under external forces on GPU.
# Typically used after an initial MC relaxation step. The Langevin
# thermostat with a small timestep keeps the system stable while
# backbone force clamping remains active. External forces file must
# be supplied separately.

_relax_MD_input = {
    # Program parameters
    "sim_type":                  "MD",
    "backend":                   "CUDA",
    "backend_precision":         "mixed",
    "use_edge":                  "1",
    "edge_n_forces":             "1",
    "CUDA_list":                 "verlet",
    "CUDA_sort_every":           "0",
    "cells_auto_optimisation":   "true",
    "verlet_skin":               "0.5",
    # Simulation parameters
    "interaction_type":          "DNA2",
    "steps":                     "1e8",
    "dt":                        "0.002",
    "ensemble":                  "nvt",
    "T":                         "20C",
    "salt_concentration":        "1.0",
    "thermostat":                "langevin",
    "newtonian_steps":           "103",
    "diff_coeff":                "2.5",
    "max_backbone_force":        "5",
    "max_backbone_force_far":    "10",
    "external_forces":           "1",
    "external_forces_file":      "external_forces.txt",
    "use_average_seq":           "1",
    # I/O parameters
    "lastconf_file":             "last_conf_MD.dat",
    "trajectory_file":           "trajectory_MD.dat",
    "energy_file":               "energy_MD.dat",
    "print_conf_interval":       "5e5",    # steps/200 for 1e8
    "print_energy_every":        "5e4",    # steps/2000 for 1e8
    "time_scale":                "linear",
    "reset_com_momentum":        "true",
    "refresh_vel":               "1",
    "restart_step_counter":      "1",
    "no_stdout_energy":          "0",
}


# ---------------------------------------------------------------------------
# Production MD protocol
# ---------------------------------------------------------------------------
# Derived from example_input_files/input_run
#
# Purpose: Long production run on GPU for data collection. Uses the
# Brownian thermostat (diffusive dynamics, appropriate for coarse-grained
# DNA) with the standard dt=0.003 timestep. External forces are disabled
# by default; enable and supply a forces file if needed.

_production_MD_input = {
    # Program parameters
    "sim_type":               "MD",
    "backend":                "CUDA",
    "backend_precision":      "mixed",
    "use_edge":               "1",
    "edge_n_forces":          "1",
    "CUDA_list":              "verlet",
    "CUDA_sort_every":        "0",
    "verlet_skin":            "0.5",
    # Simulation parameters
    "interaction_type":       "DNA2",
    "steps":                  "1e9",
    "dt":                     "0.003",
    "ensemble":               "nvt",
    "T":                      "20C",
    "salt_concentration":     "1.0",
    "thermostat":             "brownian",
    "newtonian_steps":        "103",
    "diff_coeff":             "2.5",
    "external_forces":        "0",
    "use_average_seq":        "1",
    # I/O parameters
    "lastconf_file":          "last_conf.dat",
    "trajectory_file":        "trajectory.dat",
    "energy_file":            "energy.dat",
    "print_conf_interval":    "5e5",    # steps/2000 for 1e9
    "print_energy_every":     "1e5",    # steps/10000 for 1e9
    "time_scale":             "linear",
    "refresh_vel":            "1",
    "restart_step_counter":   "1",
    "no_stdout_energy":       "0",
}


# ---------------------------------------------------------------------------
# Public getter functions (return deep copies so callers can modify freely)
# ---------------------------------------------------------------------------

def get_default_input():
    """
        Return a deep copy of the default input parameters.

        This is the legacy general-purpose parameter set from the original
        boilerplate module. Suitable as a starting point for custom simulations.

        Returns:
            (dict[str, str]) : A mutable copy of default_input_file
    """
    return deepcopy(default_input_file)


def get_relax_MC_input():
    """
        Return parameters for a Monte Carlo relaxation protocol.

        This protocol is designed for initial coarse relaxation of structures
        that may contain steric clashes (e.g. freshly converted from cadnano).

        Key characteristics:
            - sim_type: MC on CPU (no GPU required)
            - 5000 steps with large MC moves (delta = 0.22)
            - Slightly elevated temperature (30C) to aid relaxation
            - Backbone force clamping (max_backbone_force = 5)

        Returns:
            (dict[str, str]) : A mutable copy of the MC relaxation parameters
    """
    return deepcopy(_relax_MC_input)


def get_relax_MD_input():
    """
        Return parameters for an MD relaxation protocol with external forces.

        This protocol performs fine-grained relaxation under applied forces,
        typically used after an initial MC relaxation step.

        Key characteristics:
            - sim_type: MD on CUDA (GPU required)
            - 1e8 steps with small timestep (dt = 0.002) for stability
            - Langevin thermostat at 20C
            - External forces enabled (supply forces file separately)
            - Backbone force clamping active

        Returns:
            (dict[str, str]) : A mutable copy of the MD relaxation parameters
    """
    return deepcopy(_relax_MD_input)


def get_production_MD_input():
    """
        Return parameters for a production MD simulation.

        This protocol is for long data-collection runs after the system
        has been properly relaxed.

        Key characteristics:
            - sim_type: MD on CUDA (GPU required)
            - 1e9 steps with standard timestep (dt = 0.003)
            - Brownian thermostat at 20C (diffusive dynamics)
            - External forces disabled by default
            - No backbone force clamping (system should be relaxed)

        Returns:
            (dict[str, str]) : A mutable copy of the production MD parameters
    """
    return deepcopy(_production_MD_input)
