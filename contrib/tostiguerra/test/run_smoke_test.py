#!/usr/bin/env python3
"""Smoke test of the ANNaMo plugin: generate a configuration, run a short simulation, check the energies.

Run through `make test_annamo` from the build directory, which passes the paths below. In CUDA builds (--cuda) the
simulation is also run on the GPU, if there is one, and the energy computed by the GPU kernels is compared with the
CPU one.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, 'RNA_pseudoknot')
N_BEADS = 11
GPU_TOLERANCE = 1e-4  # per bead, mixed precision


def fail(msg):
    print(f'FAIL: {msg}')
    sys.exit(1)


def run(cmd, cwd):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    except OSError as e:
        fail(f'cannot run {cmd[0]}: {e.strerror}')
    if p.returncode != 0:
        log = os.path.join(cwd, 'log.dat')
        errors = [l.strip() for l in open(log) if 'ERROR' in l] if os.path.exists(log) else []
        errors += [l.strip() for l in (p.stdout + p.stderr).splitlines() if 'ERROR' in l]
        fail(f'{" ".join(cmd)} exited with code {p.returncode}' + (f': {errors[-1]}' if errors else ''))


def read_energies(path):
    """(time, potential, kinetic, total, GPU energy or None) for every line of an energy file"""
    rows = []
    for line in open(path):
        gpu = None
        if 'CUDA_energy:' in line:
            line, rest = line.split('CUDA_energy:', 1)
            gpu = float(rest.split()[0])
        s = line.split()
        if s and not s[0].startswith('#'):
            rows.append([float(x) for x in s[:4]] + [gpu])
    return rows


def simulate(workdir, args, extra):
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for f in ('topology.dat', 'dHdS_matrix.dat'):
        shutil.copy(os.path.join(EXAMPLE, f), workdir)
    lines = open(os.path.join(EXAMPLE, 'input')).read().splitlines()
    lines.append(f'plugin_search_path = {os.path.abspath(args.plugins)}')
    lines += [f'{k} = {v}' for k, v in extra.items()]
    open(os.path.join(workdir, 'input'), 'w').write('\n'.join(lines) + '\n')

    run([args.confgen, 'input', '100'], workdir)
    conf = os.path.join(workdir, 'init_conf.dat')
    if not os.path.exists(conf) or len(open(conf).read().splitlines()) != 3 + N_BEADS:
        fail(f'confGenerator did not write a configuration of {N_BEADS} beads')

    run([args.oxdna, 'input'], workdir)
    energies = read_energies(os.path.join(workdir, 'energy.dat'))
    if len(energies) < 2:
        fail('the energy file has fewer than two lines')
    if not all(math.isfinite(x) for row in energies for x in row[:4]):
        fail('non-finite energies')
    if not os.path.exists(os.path.join(workdir, 'last_conf.dat')):
        fail('no last_conf.dat')
    return energies


def gpu_available():
    try:
        return subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True).returncode == 0
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--oxdna', required=True)
    ap.add_argument('--confgen', required=True)
    ap.add_argument('--plugins', required=True, help='folder with ANNaMoInteraction.so (and CUDAANNaMoInteraction.so)')
    ap.add_argument('--work', required=True, help='scratch folder, emptied by the test')
    ap.add_argument('--cuda', choices=('mixed', 'double'), help='also run on CUDA, with this precision')
    args = ap.parse_args()
    # the simulations run in their own folders
    for k in ('oxdna', 'confgen', 'plugins', 'work'):
        setattr(args, k, os.path.abspath(getattr(args, k)))

    cpu = simulate(os.path.join(args.work, 'cpu'), args, {'backend': 'CPU'})
    print(f'CPU: {len(cpu)} energies, potential energy {cpu[0][1]:.6f} -> {cpu[-1][1]:.6f} per bead')

    if args.cuda:
        if not gpu_available():
            print('CUDA: no GPU available, skipped (the CUDA plugin was compiled)')
        else:
            cuda = simulate(os.path.join(args.work, 'cuda'), args, {
                'backend': 'CUDA', 'backend_precision': args.cuda, 'CUDA_list': 'verlet', 'CUDA_sort_every': 0,
                'use_edge': 0, 'CUDA_print_energy': 'true'})
            pairs = [(row[1], row[4]) for row in cuda if row[4] is not None]
            if not pairs:
                fail('CUDA: no GPU energies in the energy file (CUDA_print_energy)')
            worst = max(abs(gpu - host) for host, gpu in pairs)
            if not all(math.isfinite(gpu) for _, gpu in pairs) or worst > GPU_TOLERANCE:
                fail(f'CUDA: the GPU energy differs from the CPU one by up to {worst:.2e} per bead (tolerance {GPU_TOLERANCE:g})')
            print(f'CUDA: {len(cuda)} energies, max |GPU - CPU| energy {worst:.1e} per bead')

    print('PASS')


if __name__ == '__main__':
    main()
