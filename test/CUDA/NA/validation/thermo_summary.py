#!/usr/bin/env python3
"""Summarises the thermostatted runs: per system and temperature, the mean over seeds (with the standard error of the seed means, or the
largest block error if smaller) of each observable on each backend, and the deviation of each CUDA precision from the CPU in units of the
combined standard error."""
import glob
import json
import sys
from collections import defaultdict

import numpy as np

runs = defaultdict(list)
for f in sorted(glob.glob("res/thermo/*.json")):
    d = json.load(open(f))
    if d.get("returncode") != 0:
        continue
    a = d["args"]
    runs[(a["system"], a["T"], a["backend"] if a["backend"] == "CPU" else a["precision"])].append(d)

OBS = ["U", "K", "n_hb"] + ["terms/" + t for t in ("stacking", "hydrogen_bonding", "cross_stacking", "coaxial_stacking", "debye_huckel", "nonbonded_excluded_volume")]


def get(d, o):
    if o.startswith("terms/"):
        return d["terms"][o[6:]]
    return d[o]


def agg(lst, o):
    m = np.array([get(d, o)[0] for d in lst])
    e = np.array([get(d, o)[1] for d in lst])
    # mean over seeds; error: combination of the block errors (each seed mean has error e_i) and the scatter between seeds
    sem_blocks = np.sqrt((e ** 2).sum()) / len(m)
    sem_seeds = m.std(ddof=1) / np.sqrt(len(m)) if len(m) > 1 else 0.0
    return m.mean(), max(sem_blocks, sem_seeds)


cases = sorted({(s, T) for (s, T, b) in runs})
for (s, T) in cases:
    print("=== %s, T = %s  (seeds: %s)" % (s, T, {b: len(runs[(s, T, b)]) for b in ("CPU", "float", "mixed", "double") if (s, T, b) in runs}))
    base = runs.get((s, T, "CPU"))
    for o in OBS:
        row = []
        for b in ("CPU", "float", "mixed", "double"):
            l = runs.get((s, T, b))
            if not l:
                row.append("%24s" % "-")
                continue
            m, e = agg(l, o)
            z = ""
            if b != "CPU" and base:
                mb, eb = agg(base, o)
                z = " (%+.1f s)" % ((m - mb) / np.sqrt(e ** 2 + eb ** 2 + 1e-30))
            row.append("%10.4f+-%.4f%-10s" % (m, e, z))
        print("%-26s" % o + " ".join(row))
