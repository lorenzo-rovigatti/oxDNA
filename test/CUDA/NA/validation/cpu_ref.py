#!/usr/bin/env python3
"""CPU reference for the NA interaction on one configuration, computed with oxpy.

Writes an .npz with the per-nucleotide forces (lab frame), torques (body frame) and energies (each pair adds its energy to both
nucleotides, which is what the CUDA code does), the potential energy and the energy of each interaction term (summed over pairs).
Usage: cpu_ref.py <input file> <out.npz> [--terms]
"""
import os
import sys

import numpy as np

NTERMS = 8
TERM_NAMES = ["backbone", "bonded_excluded_volume", "stacking", "nonbonded_excluded_volume", "hydrogen_bonding", "cross_stacking", "coaxial_stacking", "debye_huckel"]


def main():
    if "--oxpy" in sys.argv:
        sys.path.insert(0, sys.argv[sys.argv.index("--oxpy") + 1])
    import oxpy

    inp = os.path.abspath(sys.argv[1])
    out = sys.argv[2]
    want_terms = "--terms" in sys.argv
    hybrid_only = "--hybrid-only" in sys.argv
    os.chdir(os.path.dirname(inp))
    with oxpy.Context():
        manager = oxpy.OxpyManager(os.path.basename(inp))
        ci = manager.config_info()
        interaction = ci.interaction
        parts = ci.particles()
        N = len(parts)
        for p in parts:
            p.force = np.zeros(3)
            p.torque = np.zeros(3)
        E = np.zeros(N)
        U = 0.0
        topfile = [l.split("=")[1].strip() for l in open(inp) if l.split("=")[0].strip() == "topology"][0]
        species = []
        for l in open(topfile).read().splitlines()[1:]:
            sp = l.split()
            species += [("type=RNA" in sp)] * len(sp[0])
        assert len(species) == N
        for i in range(N):
            for j in range(i + 1, N):
                if hybrid_only and (parts[i].is_bonded(parts[j]) or species[i] == species[j]):
                    continue
                e = interaction.pair_interaction(parts[i], parts[j], True, True)
                if e != 0.0:
                    E[i] += e
                    E[j] += e
                    U += e
        F = np.array([np.array(p.force) for p in parts])
        T = np.array([np.array(p.torque) for p in parts])
        terms = np.zeros(NTERMS)
        if want_terms:
            for i in range(N):
                for j in range(i + 1, N):
                    bonded = parts[i].is_bonded(parts[j])
                    for t in range(NTERMS):
                        if (t < 3) != bonded:
                            continue
                        terms[t] += interaction.pair_interaction_term(t, parts[i], parts[j], True, False)
        np.savez(out, F=F, T=T, E=E, U=U, terms=terms, N=N, btype=np.array([p.btype for p in parts]))
    print("U = %.15g (N = %d)" % (U, N))


if __name__ == "__main__":
    main()
