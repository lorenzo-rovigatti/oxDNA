#!/usr/bin/env python3
"""
plot_results.py

Turns benchmarks/results/summary.csv (produced by run_benchmark.sh) into a
comparison graph: throughput and peak VRAM of the new CUDAEdgeList
(`CUDA_list = edge`) against the existing Verlet list (with and without
`use_edge`), plus a CPU Verlet-vs-Cells reference panel.

Usage:
    python3 plot_results.py [summary.csv] [output.png]
"""
import csv
import re
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("matplotlib is required: pip install matplotlib")

CASE_RE = re.compile(r"^N(\d+)_(cpu|gpu)_(verlet|cells|verlet_edge|edge)$")

LABELS = {
    "cpu_verlet": "CPU · verlet",
    "cpu_cells": "CPU · cells",
    "gpu_verlet": "GPU · verlet",
    "gpu_verlet_edge": "GPU · verlet + use_edge",
    "gpu_edge": "GPU · edge (this PR)",
}
COLORS = {
    "cpu_verlet": "#999999",
    "cpu_cells": "#444444",
    "gpu_verlet": "#4c72b0",
    "gpu_verlet_edge": "#dd8452",
    "gpu_edge": "#55a868",
}


def load(csv_path):
    data = {}  # variant -> {N: (steps_per_sec, vram_mb)}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            m = CASE_RE.match(row["case"])
            if not m:
                continue
            n, variant = int(m.group(1)), f"{m.group(2)}_{m.group(3)}"
            try:
                sps = float(row["steps_per_sec"])
            except ValueError:
                continue
            vram = row["vram_mb"]
            vram = float(vram) if vram not in ("N/A", "") else None
            data.setdefault(variant, {})[n] = (sps, vram)
    return data


def plot(data, out_path):
    fig, (ax_speed, ax_vram) = plt.subplots(1, 2, figsize=(11, 4.5))

    for variant in ("cpu_verlet", "cpu_cells", "gpu_verlet", "gpu_verlet_edge", "gpu_edge"):
        if variant not in data:
            continue
        sizes = sorted(data[variant])
        sps = [data[variant][n][0] for n in sizes]
        ax_speed.plot(sizes, sps, marker="o", label=LABELS[variant], color=COLORS[variant])

    ax_speed.set_xscale("log")
    ax_speed.set_yscale("log")
    ax_speed.set_xlabel("nucleotides (N)")
    ax_speed.set_ylabel("steps / s (higher is better)")
    ax_speed.set_title("Throughput: edge list vs verlet list")
    ax_speed.legend(fontsize=8)
    ax_speed.grid(True, which="both", alpha=0.3)

    for variant in ("gpu_verlet", "gpu_verlet_edge", "gpu_edge"):
        if variant not in data:
            continue
        sizes = sorted(n for n in data[variant] if data[variant][n][1] is not None)
        vram = [data[variant][n][1] for n in sizes]
        if sizes:
            ax_vram.plot(sizes, vram, marker="o", label=LABELS[variant], color=COLORS[variant])

    ax_vram.set_xscale("log")
    ax_vram.set_xlabel("nucleotides (N)")
    ax_vram.set_ylabel("peak extra VRAM (MB)")
    ax_vram.set_title("GPU memory: edge list vs verlet list")
    ax_vram.legend(fontsize=8)
    ax_vram.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results" / "summary.csv"
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_path.parent / "benchmark.png"
    plot(load(csv_path), out_path)
