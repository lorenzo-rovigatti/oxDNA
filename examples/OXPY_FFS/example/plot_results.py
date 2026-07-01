#!/usr/bin/env python3
"""
Analyse the OXPY_FFS example output and generate reference results.

This script reads the ffs.log and progress.json files produced by the FLUX
and SHOOT_* stages, estimates the FFS rate and its uncertainty, and writes:

    ffs_reference_results.dat
    ffs_cumulative_rate.png

The plot is intended as a compact diagnostic/reference plot showing how the
final FFS rate estimate is built from the flux and successive conditional
crossing probabilities.
"""

import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


STAGES = [
    ("FLUX", "Flux", "flux"),
    ("SHOOT_01_dist_m3_m2", "P(-3->-2)", "shoot"),
    ("SHOOT_02_dist_m2_m1", "P(-2->-1)", "shoot"),
    ("SHOOT_03_dist_m1_0", "P(-1->0)", "shoot"),
    ("SHOOT_04_bond_0_1", "P(0->1)", "shoot"),
]


def parse_flux_log(path):
    text = path.read_text()

    mean_match = re.search(r"Average number of timesteps.*:\s*([0-9.eE+-]+)", text)
    std_match = re.search(
        r"Average number of timesteps standard deviation:\s*([0-9.eE+-]+)",
        text,
    )
    flux_match = re.search(r"Initial flux:\s*([0-9.eE+-]+)", text)

    if mean_match is None or std_match is None or flux_match is None:
        raise RuntimeError(f"Could not parse flux summary from {path}")

    return float(mean_match.group(1)), float(std_match.group(1)), float(flux_match.group(1))


def parse_shoot_log(path):
    text = path.read_text()

    match = re.search(
        r"nsuccesses:\s*(\d+)\s+nattempts:\s*(\d+)\s+success_prob:\s*([0-9.eE+-]+)",
        text,
    )

    if match is None:
        raise RuntimeError(f"Could not parse shooting summary from {path}")

    return int(match.group(1)), int(match.group(2)), float(match.group(3))


def read_progress(path):
    with open(path) as f:
        return json.load(f)


def main():
    base = Path(__file__).resolve().parent

    rows = []
    cumulative_rate = 1.0
    cumulative_relative_variance = 0.0

    for dirname, label, kind in STAGES:
        stage_dir = base / dirname
        log_path = stage_dir / "ffs.log"
        progress_path = stage_dir / "progress.json"

        if not log_path.exists():
            raise FileNotFoundError(f"Missing {log_path}")
        if not progress_path.exists():
            raise FileNotFoundError(f"Missing {progress_path}")

        progress = read_progress(progress_path)

        if kind == "flux":
            mean_t, std_t, value = parse_flux_log(log_path)
            successes = int(progress["success_count"])
            attempts = None

            # The flux is 1 / <tau>. The uncertainty is obtained by propagating
            # the standard error on the mean first-crossing time.
            sem_t = std_t / math.sqrt(successes)
            error = sem_t / (mean_t * mean_t)

        else:
            successes, attempts, value = parse_shoot_log(log_path)
            mean_t = None
            std_t = None

            # Binomial standard error for a conditional crossing probability.
            error = math.sqrt(value * (1.0 - value) / attempts)

        relative_error = error / value

        cumulative_rate *= value
        cumulative_relative_variance += relative_error**2
        cumulative_error = cumulative_rate * math.sqrt(cumulative_relative_variance)

        rows.append(
            {
                "stage": dirname,
                "label": label,
                "kind": kind,
                "value": value,
                "error": error,
                "relative_error": relative_error,
                "successes": successes,
                "attempts": attempts,
                "mean_timesteps": mean_t,
                "std_timesteps": std_t,
                "cumulative_rate": cumulative_rate,
                "cumulative_rate_error": cumulative_error,
            }
        )

    final_rate = rows[-1]["cumulative_rate"]
    final_rate_error = rows[-1]["cumulative_rate_error"]

    average_waiting_time = 1.0 / final_rate
    average_waiting_time_error = final_rate_error / (final_rate * final_rate)

    with open("ffs_reference_results.dat", "w") as f:
        f.write("# Reference FFS results\n")
        f.write("# Values are estimated from the example run.\n")
        f.write("# Errors correspond to one standard error.\n\n")

        for row in rows:
            f.write(
                f"{row['label']:12s} "
                f"value = {row['value']:.8e} "
                f"error = {row['error']:.8e} "
                f"rel_error = {row['relative_error']:.6f}"
            )

            if row["kind"] == "flux":
                f.write(f" crossings = {row['successes']}")
                f.write(f" mean_timesteps = {row['mean_timesteps']:.8e}")
                f.write(f" std_timesteps = {row['std_timesteps']:.8e}")
            else:
                f.write(f" successes = {row['successes']}")
                f.write(f" attempts = {row['attempts']}")

            f.write("\n")

        f.write("\n")
        f.write(
            f"Final rate = ({final_rate:.8e} +/- "
            f"{final_rate_error:.8e}) timestep^-1\n"
        )
        f.write(
            f"Average waiting time = ({average_waiting_time:.8e} +/- "
            f"{average_waiting_time_error:.8e}) timesteps\n"
        )

    x = list(range(len(rows)))
    labels = [row["label"] for row in rows]
    y = [row["cumulative_rate"] for row in rows]
    yerr = [row["cumulative_rate_error"] for row in rows]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))

    ax.errorbar(x, y, yerr=yerr, marker="o", linestyle="-", capsize=4)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel(r"Cumulative rate estimate [$\mathrm{timestep}^{-1}$]")
    ax.set_title("Construction of the FFS rate estimate")

    ax.text(
        0.03,
        0.05,
        (
            f"Final rate = ({final_rate:.2e} +/- "
            f"{final_rate_error:.1e}) timestep$^{{-1}}$"
        ),
        transform=ax.transAxes,
        va="bottom",
    )

    fig.tight_layout()
    fig.savefig("ffs_cumulative_rate.png", dpi=300)

    print("Wrote ffs_reference_results.dat")
    print("Wrote ffs_cumulative_rate.png")
    print(
        f"Final rate = ({final_rate:.8e} +/- "
        f"{final_rate_error:.8e}) timestep^-1"
    )
    print(
        f"Average waiting time = ({average_waiting_time:.8e} +/- "
        f"{average_waiting_time_error:.8e}) timesteps"
    )


if __name__ == "__main__":
    main()
