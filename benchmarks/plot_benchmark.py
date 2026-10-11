#!/usr/bin/env python3
"""Plot the oxDNA Metal benchmark results (ErikPoppleton/oxDNA_performance systems)."""
import csv, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "results.csv"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(OUTDIR, exist_ok=True)

rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        try:
            r["nuc"] = int(r["nuc"])
            r["ms_per_step"] = float(r["ms_per_step"])
            r["total_1e6_steps_h"] = float(r["total_1e6_steps_h"])
        except ValueError:
            continue
        if r["ms_per_step"] <= 0:
            continue
        rows.append(r)

configs = [
    ("cpu",          "CPU (double)",            "#444444", "o", "-"),
    ("float_verlet", "Metal float · verlet",    "#1f77b4", "s", "-"),
    ("float_edge",   "Metal float · edge",      "#1f77b4", "D", "--"),
    ("mixed_verlet", "Metal mixed · verlet",    "#d62728", "s", "-"),
    ("mixed_edge",   "Metal mixed · edge",      "#d62728", "D", "--"),
    ("hardmixed",    "Metal hardmixed (6 thr)", "#2ca02c", "^", "-"),
]

def series(name):
    d = sorted([r for r in rows if r["config"] == name], key=lambda r: r["nuc"])
    return np.array([r["nuc"] for r in d]), np.array([r["ms_per_step"] for r in d])

# ---- 1. ms/step vs system size ----
fig, ax = plt.subplots(figsize=(8, 5.5))
for name, label, color, marker, ls in configs:
    x, y = series(name)
    if len(x):
        ax.loglog(x, y, marker=marker, color=color, ls=ls, label=label, ms=7)
ax.set_xlabel("system size  (nucleotides)")
ax.set_ylabel("wall time per MD step  (ms)")
ax.set_title("oxDNA DNA2 — time per step, Apple M4\n(ErikPoppleton/oxDNA_performance duplex boxes, dt=0.003)")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "ms_per_step.png"), dpi=130)
plt.close(fig)

# ---- 2. speed-up vs CPU ----
cpu = {r["nuc"]: r["ms_per_step"] for r in rows if r["config"] == "cpu"}
fig, ax = plt.subplots(figsize=(8, 5.5))
for name, label, color, marker, ls in configs:
    if name == "cpu":
        continue
    x, y = series(name)
    xs, sp = [], []
    for xi, yi in zip(x, y):
        if xi in cpu:
            xs.append(xi); sp.append(cpu[xi] / yi)
    if xs:
        ax.semilogx(xs, sp, marker=marker, color=color, ls=ls, label=label, ms=7)
ax.axhline(1, color="#444", lw=0.8)
ax.set_xlabel("system size  (nucleotides)")
ax.set_ylabel("speed-up vs CPU backend  (×)")
ax.set_title("oxDNA DNA2 — Metal speed-up over the double-precision CPU backend")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "speedup_vs_cpu.png"), dpi=130)
plt.close(fig)

# ---- 3. edge vs verlet ----
fig, ax = plt.subplots(figsize=(8, 5.5))
for prec, color in [("float", "#1f77b4"), ("mixed", "#d62728")]:
    xv, yv = series(f"{prec}_verlet")
    xe, ye = series(f"{prec}_edge")
    m = {xi: yi for xi, yi in zip(xe, ye)}
    xs, gain = [], []
    for xi, yi in zip(xv, yv):
        if xi in m:
            xs.append(xi); gain.append(yi / m[xi])
    if xs:
        ax.semilogx(xs, gain, marker="D", color=color, label=f"{prec}: verlet / edge", ms=7)
ax.axhline(1, color="#444", lw=0.8)
ax.set_xlabel("system size  (nucleotides)")
ax.set_ylabel("edge-list speed-up over verlet  (×)")
ax.set_title("oxDNA Metal — edge list vs verlet list")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "edge_vs_verlet.png"), dpi=130)
plt.close(fig)

# ---- 4. projected total wall time for 1e6 steps ----
fig, ax = plt.subplots(figsize=(8, 5.5))
for name, label, color, marker, ls in configs:
    d = sorted([r for r in rows if r["config"] == name], key=lambda r: r["nuc"])
    if d:
        ax.loglog([r["nuc"] for r in d], [r["total_1e6_steps_h"] for r in d],
                  marker=marker, color=color, ls=ls, label=label, ms=7)
ax.set_xlabel("system size  (nucleotides)")
ax.set_ylabel("projected wall time for 1e6 steps  (hours)")
ax.set_title("oxDNA DNA2 — projected cost of a 1e6-step run")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "total_1e6_steps.png"), dpi=130)
plt.close(fig)

# ---- text table ----
print(f"{'system':>8} {'nuc':>8} " + " ".join(f"{c[0]:>14}" for c in configs))
sizes = sorted(set(r["nuc"] for r in rows))
for s in sizes:
    line = f"{'':>8} {s:>8} "
    for name, *_ in configs:
        v = [r["ms_per_step"] for r in rows if r["nuc"] == s and r["config"] == name]
        line += f"{(v[0] if v else float('nan')):>14.4f} "
    print(line)
print("\nwrote: ms_per_step.png  speedup_vs_cpu.png  edge_vs_verlet.png  total_1e6_steps.png")
