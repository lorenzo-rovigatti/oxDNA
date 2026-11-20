#!/usr/bin/env python3
"""Thermostatted (john) run of the NA interaction on CPU or CUDA, and the averages of observables.

The observables are the ones that an oxNA user looks at: the kinetic energy per nucleotide (i.e. the temperature), the potential energy per nucleotide,
the energy of each interaction term per nucleotide and the number of hydrogen bonds (hb_list), with block-averaged standard errors.

Usage: thermo_run.py --system hybrid_duplex --T 50C --steps 2000000 --seed 1 --backend CPU|CUDA [--precision ..] --bin <oxDNA> --tag name --out out.json
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

RUNS = Path(os.environ.get("NA_VALIDATION_RUNS", "runs")).resolve()
DATA = Path(os.environ.get("NA_VALIDATION_DATA", "data")).resolve()
TERMS = ["backbone", "bonded_excluded_volume", "stacking", "nonbonded_excluded_volume", "hydrogen_bonding", "cross_stacking", "coaxial_stacking", "debye_huckel"]


def block_stats(x, nblocks=10):
    x = np.asarray(x, float)
    n = len(x) // nblocks * nblocks
    b = x[:n].reshape(nblocks, -1).mean(axis=1)
    return float(x.mean()), float(b.std(ddof=1) / np.sqrt(nblocks))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--T", default="25C")
    ap.add_argument("--steps", type=int, default=2000000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--backend", required=True)
    ap.add_argument("--precision", default="mixed")
    ap.add_argument("--bin", required=True)
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--dt", type=float, default=0.003)
    ap.add_argument("--salt", type=float, default=0.5)
    ap.add_argument("--every", type=int, default=1000)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", action="append")
    args = ap.parse_args()

    work = RUNS / "thermo" / ("%s_s%d" % (args.tag, args.seed))
    work.mkdir(parents=True, exist_ok=True)
    start = work / "start.dat"
    lines = (DATA / args.system / "init.dat").read_text().splitlines()
    # start from the relaxed configuration of the CPU run (without its velocities being used: refresh_vel)
    traj = DATA / args.system / "traj.dat"
    if traj.exists():
        tl = traj.read_text().splitlines()
        idx = [i for i, l in enumerate(tl) if l.startswith("t =")]
        lines = tl[idx[-1]:]
    start.write_text("\n".join(lines) + "\n")
    inp = [
        "backend = %s" % args.backend, "sim_type = MD", "interaction_type = NA", "use_average_seq = true",
        "salt_concentration = %g" % args.salt, "T = %s" % args.T,
        "topology = %s" % (DATA / args.system / "init.top"), "conf_file = %s" % start,
        "thermostat = john", "newtonian_steps = 103", "diff_coeff = 2.5", "dt = %g" % args.dt, "steps = %d" % args.steps,
        "verlet_skin = 0.5", "refresh_vel = 1", "fix_diffusion = false", "seed = %d" % args.seed,
        "no_stdout_energy = 1", "print_energy_every = %d" % args.every, "print_conf_interval = 1000000000",
        "trajectory_file = /dev/null", "lastconf_file = %s" % (work / "last.dat"),
        "energy_file = %s" % (work / "energy.dat"), "log_file = %s" % (work / "log.dat"),
        "restart_step_counter = 1", "time_scale = linear",
        "data_output_1 = {\n name = %s\n print_every = %d\n col_1 = {\n type = potential_energy\n split = true\n }\n}" % (work / "split.dat", args.every),
        "data_output_2 = {\n name = %s\n print_every = %d\n col_1 = {\n type = hb_list\n only_count = true\n }\n}" % (work / "hb.dat", args.every),
    ]
    if args.backend == "CUDA":
        inp += ["backend_precision = %s" % args.precision, "CUDA_device = %d" % args.gpu, "CUDA_sort_every = 0"]
    inp += args.set or []
    f = work / "input"
    f.write_text("\n".join(inp) + "\n")
    t0 = time.time()
    r = subprocess.run([args.bin, str(f)], capture_output=True, text=True, cwd=str(work))
    wall = time.time() - t0
    res = {"args": vars(args), "returncode": r.returncode, "wall": wall}
    if r.returncode != 0:
        res["error"] = (r.stderr or r.stdout)[-400:]
        Path(args.out).write_text(json.dumps(res, indent=1))
        print("FAILED", res["error"])
        return
    e = np.loadtxt(work / "energy.dat")
    s = np.loadtxt(work / "split.dat")
    h = np.loadtxt(work / "hb.dat")
    h = h.reshape(-1) if h.ndim == 1 else h[:, -1]
    n0 = len(e) // 5  # discard the first 20 %
    res["samples"] = int(len(e) - n0)
    res["U"] = block_stats(e[n0:, 1])
    res["K"] = block_stats(e[n0:, 2])
    res["E"] = block_stats(e[n0:, 3])
    res["terms"] = {TERMS[i]: block_stats(s[n0:, i]) for i in range(len(TERMS))}
    res["n_hb"] = block_stats(h[n0:])
    Path(args.out).write_text(json.dumps(res, indent=1))
    print("%s s%d wall=%.0fs U=%.4f+-%.4f K=%.4f+-%.4f n_hb=%.3f+-%.3f" % (args.tag, args.seed, wall, *res["U"], *res["K"], *res["n_hb"]), flush=True)


if __name__ == "__main__":
    main()
