#!/usr/bin/env python3
"""NVE energy conservation versus time step for the NA interaction, on CPU or CUDA.

The system starts from the last frame of the trajectory of a thermostatted CPU run (positions and velocities), thermostat = no. The
energy file is written by oxDNA itself: the potential energy is computed on the CPU, in double precision, from the (GPU) positions,
and the kinetic energy from the velocities, so a force that is not the gradient of the energy shows up as a drift of the total.

Usage: nve_run.py --system large --backend CPU|CUDA [--precision float|mixed|double] --dts 0.001,0.002 --time 20 --bin <oxDNA> --out result.json
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


def last_frame(traj, out):
    lines = Path(traj).read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("t =")]
    Path(out).write_text("\n".join(lines[starts[-1]:]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--backend", required=True)
    ap.add_argument("--precision", default="double")
    ap.add_argument("--dts", required=True)
    ap.add_argument("--time", type=float, default=20.0)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--salt", type=float, default=0.5)
    ap.add_argument("--T", default="25C")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--timeout", type=int, default=3000)
    ap.add_argument("--nsamp", type=int, default=100)
    ap.add_argument("--set", action="append")
    args = ap.parse_args()

    work = RUNS / "nve" / args.tag
    work.mkdir(parents=True, exist_ok=True)
    start = work / "start.dat"
    traj = DATA / args.system / "traj.dat"
    last_frame(traj if traj.exists() else DATA / args.system / "init.dat", start)
    N = int(Path(DATA / args.system / "init.top").read_text().split()[0])

    results = []
    for dt in [float(x) for x in args.dts.split(",")]:
        steps = int(round(args.time / dt))
        every = max(1, steps // args.nsamp)
        tag = "dt%g" % dt
        lines = [
            "backend = %s" % args.backend, "sim_type = MD", "interaction_type = NA", "use_average_seq = true",
            "salt_concentration = %g" % args.salt, "T = %s" % args.T,
            "topology = %s" % (DATA / args.system / "init.top"), "conf_file = %s" % start,
            "thermostat = no", "dt = %g" % dt, "steps = %d" % steps, "verlet_skin = 0.5", "refresh_vel = 0", "fix_diffusion = false",
            "no_stdout_energy = 1", "print_energy_every = %d" % every, "print_conf_interval = 1000000000",
            "trajectory_file = /dev/null", "lastconf_file = %s" % (work / ("last_%s.dat" % tag)),
            "energy_file = %s" % (work / ("energy_%s.dat" % tag)), "log_file = %s" % (work / ("log_%s.dat" % tag)),
            "restart_step_counter = 1", "time_scale = linear", "seed = 1",
        ]
        if args.backend == "CUDA":
            lines += ["backend_precision = %s" % args.precision, "CUDA_device = %d" % args.gpu, "CUDA_sort_every = 0"]
        lines += args.set or []
        inp = work / ("in_%s" % tag)
        inp.write_text("\n".join(lines) + "\n")
        efile = work / ("energy_%s.dat" % tag)
        if efile.exists():
            efile.unlink()
        t0 = time.time()
        status = "ok"
        try:
            r = subprocess.run([args.bin, str(inp)], capture_output=True, text=True, timeout=args.timeout, cwd=str(work))
            if r.returncode != 0:
                status = "crash: " + (r.stderr or r.stdout).strip().splitlines()[-1][:160] if (r.stderr or r.stdout).strip() else "crash"
        except subprocess.TimeoutExpired:
            status = "timeout"
        wall = time.time() - t0
        res = {"dt": dt, "steps": steps, "status": status, "wall": wall}
        try:
            d = np.loadtxt(efile)
            d = np.atleast_2d(d)
            t, U, K, E = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
            res["samples"] = int(len(t))
            res["complete"] = bool(abs(t[-1] - args.time) < 2 * args.time / args.nsamp and np.isfinite(E).all())
            res["K_mean"] = float(K.mean())
            res["U_mean"] = float(U.mean())
            if len(t) > 3 and np.isfinite(E).all():
                A = np.vstack([t, np.ones_like(t)]).T
                slope, icpt = np.linalg.lstsq(A, E, rcond=None)[0]
                resid = E - (slope * t + icpt)
                res["drift_per_N_per_time"] = float(slope)
                res["drift_over_K"] = float(slope * args.time / K.mean())
                res["fluctuation"] = float(resid.std())
                res["E_first"] = float(E[0])
                res["E_last"] = float(E[-1])
                res["max_dev"] = float(np.abs(E - E[0]).max())
        except Exception as e:  # noqa
            res["parse_error"] = str(e)
        results.append(res)
        print("%s dt=%g %s wall=%.0fs %s" % (args.tag, dt, status, wall, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in res.items() if k in ("complete", "K_mean", "drift_per_N_per_time", "drift_over_K", "fluctuation", "max_dev")}), flush=True)
        Path(args.out).write_text(json.dumps({"args": vars(args), "N": N, "results": results}, indent=1))


if __name__ == "__main__":
    main()
