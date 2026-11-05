#!/usr/bin/env python3
"""CPU versus CUDA parity of the NA interaction on one or more configurations.

CPU reference: oxpy (cpu_ref.py), forces, body-frame torques and energies per nucleotide from a loop over all pairs.
GPU: oxDNA built with CUDA, with CUDA_dump_file set, which dumps the same per-nucleotide quantities after the first
force evaluation (the energy of a nucleotide is the sum over its pairs, each pair is counted twice, as in the CPU loop).

Usage: parity.py --bin <CUDA oxDNA> --precision float|mixed|double --oxpy <dir with oxpy> --systems a,b,c [--frames N] [--salt x]
       [--seqdep] [--set key=value ...] [--tag name] [--json out]
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RUNS = Path(os.environ.get("NA_VALIDATION_RUNS", "runs")).resolve()
DATA = Path(os.environ.get("NA_VALIDATION_DATA", "data")).resolve()


def split_traj(traj, outdir, nmax):
    lines = Path(traj).read_text().splitlines()
    # each configuration has 3 header lines and N lines
    frames, i = [], 0
    n = None
    top = None
    while i < len(lines):
        j = i + 3
        while j < len(lines) and not lines[j].startswith("t ="):
            j += 1
        frames.append(lines[i:j])
        i = j
    sel = frames[-nmax:] if nmax else frames
    out = []
    for k, f in enumerate(sel):
        p = Path(outdir) / ("frame_%02d.dat" % k)
        p.write_text("\n".join(f) + "\n")
        out.append(p)
    return out


def make_input(system, conf, backend, args, workdir, dump=None):
    d = DATA / system
    lines = [
        "backend = %s" % backend,
        "sim_type = MD",
        "interaction_type = NA",
        "salt_concentration = %g" % args.salt,
        "T = %s" % args.T,
        "topology = %s" % (d / "init.top"),
        "conf_file = %s" % conf,
        "thermostat = no",
        "dt = 0.001",
        "steps = 1",
        "verlet_skin = 0.5",
        "refresh_vel = 1",
        "fix_diffusion = false",
        "no_stdout_energy = 1",
        "print_energy_every = 1000000",
        "print_conf_interval = 1000000000",
        "trajectory_file = /dev/null",
        "lastconf_file = %s" % (Path(workdir) / ("last_%s.dat" % backend)),
        "energy_file = %s" % (Path(workdir) / ("energy_%s.dat" % backend)),
        "log_file = %s" % (Path(workdir) / ("log_%s.dat" % backend)),
        "restart_step_counter = 1",
        "time_scale = linear",
        "diff_coeff = 2.5",
        "newtonian_steps = 103",
        "seed = 1",
    ]
    if args.seqdep:
        P = REPO
        lines += ["use_average_seq = false", "seq_dep_file_DNA = %s" % (P / "oxDNA2_sequence_dependent_parameters.txt"),
                  "seq_dep_file_RNA = %s" % (P / "rna_sequence_dependent_parameters.txt"), "seq_dep_file_NA = %s" % (P / "NA_sequence_dependent_parameters.txt")]
    else:
        lines += ["use_average_seq = true"]
    if backend == "CUDA":
        if args.hybrid_only:
            lines += ["CUDA_NA_dump_hybrid_only = true"]
        lines += ["backend_precision = %s" % args.precision, "CUDA_device = %d" % args.gpu, "CUDA_sort_every = 0", "CUDA_dump_file = %s" % dump, "threads_per_block = 64"]
    for kv in args.set or []:
        lines.append(kv)
    return "\n".join(lines) + "\n"


def load_dump(path, N):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]
    idx = np.array([int(r[0]) for r in rows])
    vals = np.array([[float(x) for x in r[1:]] for r in rows])
    o = np.argsort(idx)
    vals = vals[o]
    assert len(vals) == N
    return vals[:, 0:3], vals[:, 4:7], vals[:, 3], vals[:, 7]


def dev(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.abs(a - b)
    mx = np.abs(b).max()
    r = {"max_abs": float(d.max()), "max_ref": float(mx), "norm_max": float(d.max() / (mx + 1e-300))}
    if a.ndim == 2:
        # per-nucleotide relative deviation, with a floor of 1e-3 times the rms of the reference so that tiny vectors do not dominate
        nd = np.linalg.norm(a - b, axis=1)
        nb = np.linalg.norm(b, axis=1)
        floor = 1e-3 * np.sqrt(np.mean(nb ** 2))
        r["rel_max"] = float((nd / np.maximum(nb, floor)).max())
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--precision", default="float")
    ap.add_argument("--oxpy", required=True)
    ap.add_argument("--systems", required=True)
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--salt", type=float, default=0.5)
    ap.add_argument("--T", default="25C")
    ap.add_argument("--seqdep", action="store_true")
    ap.add_argument("--terms", action="store_true")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--set", action="append")
    ap.add_argument("--tag", default="par")
    ap.add_argument("--json")
    ap.add_argument("--hybrid-only", action="store_true")
    ap.add_argument("--cpu-set", action="append", help="key=value only for the CPU input")
    args = ap.parse_args()

    results = []
    for system in args.systems.split(","):
        work = RUNS / args.tag / system
        work.mkdir(parents=True, exist_ok=True)
        confs = [DATA / system / "init.dat"] + (split_traj(DATA / system / "traj.dat", work, args.frames) if (DATA / system / "traj.dat").exists() else [])
        for ci, conf in enumerate(confs):
            tag = "%s_c%d" % (system, ci)
            cpu_in = work / ("cpu_%d.in" % ci)
            cpu_in.write_text(make_input(system, conf, "CPU", args, work) + "".join(kv + "\n" for kv in (args.cpu_set or [])))
            gpu_in = work / ("gpu_%d.in" % ci)
            dump = work / ("dump_%d.txt" % ci)
            gpu_in.write_text(make_input(system, conf, "CUDA", args, work, dump))
            npz = work / ("cpu_%d.npz" % ci)
            if not npz.exists() or args.cpu_set or args.hybrid_only:
                r = subprocess.run([sys.executable, str(HERE / "cpu_ref.py"), str(cpu_in), str(npz), "--oxpy", args.oxpy] + (["--terms"] if args.terms else []) + (["--hybrid-only"] if args.hybrid_only else []), capture_output=True, text=True)
                if r.returncode != 0:
                    print(r.stdout[-1500:], r.stderr[-2500:])
                    raise SystemExit("CPU reference failed for " + tag)
            ref = np.load(npz)
            N = int(ref["N"])
            g = subprocess.run([args.bin, str(gpu_in)], capture_output=True, text=True, cwd=str(work))
            if g.returncode != 0:
                print(g.stdout[-1500:], g.stderr[-2500:])
                raise SystemExit("GPU run failed for " + tag)
            F, T, E, HB = load_dump(dump, N)
            res = {"tag": tag, "N": N, "U_cpu": float(ref["U"]), "U_gpu": float(E.sum() / 2)}
            res["U_rel"] = abs(res["U_cpu"] - res["U_gpu"]) / max(abs(res["U_cpu"]), 1e-300)
            res["force"] = dev(F, ref["F"])
            res["torque"] = dev(T, ref["T"])
            res["energy"] = dev(E, ref["E"])
            if args.terms:
                res["hb_cpu"] = float(ref["terms"][4])
                res["hb_gpu"] = float(HB.sum() / 2)
            results.append(res)
            print("%-22s N=%4d U_cpu=% .10g U_rel=%.2e F=%.2e (rel %.1e) T=%.2e (rel %.1e) E=%.2e" % (tag, N, res["U_cpu"], res["U_rel"], res["force"]["norm_max"], res["force"]["rel_max"], res["torque"]["norm_max"], res["torque"]["rel_max"], res["energy"]["norm_max"]), flush=True)
    summary = {k: max(r[k]["norm_max"] for r in results) for k in ("force", "torque", "energy")}
    summary["U_rel"] = max(r["U_rel"] for r in results)
    summary["force_rel"] = max(r["force"]["rel_max"] for r in results)
    summary["torque_rel"] = max(r["torque"]["rel_max"] for r in results)
    print("WORST", json.dumps(summary))
    if args.json:
        Path(args.json).write_text(json.dumps({"args": vars(args), "results": results, "worst": summary}, indent=1))


if __name__ == "__main__":
    main()
