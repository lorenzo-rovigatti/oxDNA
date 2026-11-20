#!/usr/bin/env python3
"""
Tests of the CUDA implementation of the NA interaction (DNA2 + RNA2 + DNA/RNA hybrid), which compare it with the CPU one.

They need a GPU and an oxDNA binary compiled with -DCUDA=ON, which is used both for the CPU and for the CUDA runs. If the binary cannot
run on a GPU the tests are skipped. Usage:

    python3 test_cuda_na.py /path/to/oxDNA [-v]

or, with pytest, OXDNA_BIN=/path/to/oxDNA pytest test_cuda_na.py

The precision of the CUDA backend (float, mixed or double) can be chosen with the OXDNA_CUDA_PRECISION environment variable (defaults to mixed).
For double precision the binary has to be compiled with -DCUDA_DOUBLE=ON.

The CUDA runs are performed with the neighbour-matrix kernel and with the edge-based one (use_edge = true, which is not available in double
precision), and the two are also compared directly through the per-nucleotide dump of forces, torques and energies (CUDA_dump_file).

The configurations are small systems that were relaxed and sampled with the CPU code: a DNA duplex, an RNA duplex, a DNA:RNA hybrid
with a single-stranded overhang, and a 179-nucleotide system made of DNA and RNA duplexes, overhangs, hairpins and single strands.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ["dna_duplex", "rna_duplex", "hybrid_overhang", "mixed"]

OXDNA = os.path.abspath(os.environ["OXDNA_BIN"]) if os.environ.get("OXDNA_BIN") else None
PRECISION = os.environ.get("OXDNA_CUDA_PRECISION", "mixed")
# the CPU code uses interpolated meshes for the angular modulations, the GPU code the analytic functions: the energies of the two codes
# differ by ~1e-5, and the forces by up to ~5e-3 of the largest force (also in double precision) because of the meshes
ENERGY_TOLERANCE = 1e-4
POSITION_TOLERANCE = {"double": 1e-3, "mixed": 5e-3, "float": 5e-3}
# the edge-based and the neighbour-matrix kernels evaluate the same pair interactions but sum them in a different order, so their results
# differ only by round-off, relative to the largest force/torque/energy
EDGE_TOLERANCE = 1e-4

# extra input options of the CUDA runs: neighbour matrix, edge list, edge list with several force slots per nucleotide
CUDA_MODES = [("matrix", "")]
if PRECISION != "double":
    CUDA_MODES += [("edge", "use_edge = true\n"), ("edge_4_forces", "use_edge = true\nedge_n_forces = 4\n")]

BASE_INPUT = """
interaction_type = NA
sim_type = MD
salt_concentration = 0.5
use_average_seq = true
T = 25C
topology = {top}
conf_file = {conf}
verlet_skin = 0.5
dt = {dt}
steps = {steps}
thermostat = no
refresh_vel = 0
fix_diffusion = false
no_stdout_energy = 1
print_energy_every = {every}
print_conf_interval = 1000000
restart_step_counter = 1
time_scale = linear
trajectory_file = /dev/null
energy_file = energy.dat
lastconf_file = last.dat
log_file = log.dat
seed = 1
"""


def run(case, backend, workdir, steps=1, dt=0.001, every=1, extra=""):
    """Runs oxDNA on one of the configurations and returns the completed process"""
    d = HERE / case
    text = BASE_INPUT.format(top=d / "init.top", conf=d / "init.dat", dt=dt, steps=steps, every=every)
    text += "backend = %s\n" % backend
    if backend == "CUDA":
        text += "backend_precision = %s\nCUDA_print_energy = true\n" % PRECISION
    text += extra
    inp = Path(workdir) / "input"
    inp.write_text(text)
    return subprocess.run([OXDNA, "input"], cwd=workdir, capture_output=True, text=True)


def read_energy_file(workdir):
    """Returns the lines of the energy file as lists of floats, and the GPU energy (if present) of each line"""
    rows = []
    for line in (Path(workdir) / "energy.dat").read_text().splitlines():
        m = re.search(r"CUDA_energy:\s*(\S+)", line)
        gpu = float(m.group(1)) if m else None
        rows.append(([float(x) for x in line.split("CUDA_energy:")[0].split()], gpu))
    return rows


def read_positions(path):
    lines = Path(path).read_text().splitlines()[3:]
    return [[float(x) for x in l.split()[:3]] for l in lines if l.strip()]


class CUDANATest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if OXDNA is None or not Path(OXDNA).exists():
            raise unittest.SkipTest("oxDNA binary not found: pass its path as the first argument")
        with tempfile.TemporaryDirectory() as d:
            r = run("dna_duplex", "CUDA", d)
            if r.returncode != 0:
                raise unittest.SkipTest("the CUDA backend cannot be used (%s)" % (r.stderr or r.stdout).strip()[-200:])

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_energy_cpu_vs_gpu(self):
        """The potential energy per nucleotide computed on the GPU is the one computed by the CPU code on the same configuration"""
        for case in CASES:
            for mode, extra in CUDA_MODES:
                with self.subTest(case=case, mode=mode):
                    r = run(case, "CUDA", self.tmp, extra=extra)
                    self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
                    # the energy of the GPU is printed from the first step on: both energies refer to the configuration after one step
                    (values, gpu_energy) = read_energy_file(self.tmp)[-1]
                    cpu_energy = values[1]
                    self.assertIsNotNone(gpu_energy, "the GPU energy was not printed")
                    self.assertLess(abs(gpu_energy - cpu_energy), ENERGY_TOLERANCE * abs(cpu_energy), "CPU: %g, GPU: %g" % (cpu_energy, gpu_energy))

    def test_nve_trajectories(self):
        """Short NVE trajectories on the CPU and on the GPU, started from the same configuration, stay together (this checks forces and torques)"""
        steps, dt = 300, 0.001
        for case in CASES:
            for mode, extra in CUDA_MODES:
                with self.subTest(case=case, mode=mode):
                    cpu_dir = Path(self.tmp) / case / mode / "cpu"
                    gpu_dir = Path(self.tmp) / case / mode / "gpu"
                    cpu_dir.mkdir(parents=True)
                    gpu_dir.mkdir(parents=True)
                    for backend, d in (("CPU", cpu_dir), ("CUDA", gpu_dir)):
                        r = run(case, backend, d, steps=steps, dt=dt, every=steps, extra=extra if backend == "CUDA" else "")
                        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
                    cpu_rows = read_energy_file(cpu_dir)
                    gpu_rows = read_energy_file(gpu_dir)
                    # total energy per nucleotide at the end of the run
                    self.assertLess(abs(cpu_rows[-1][0][3] - gpu_rows[-1][0][3]), 1e-3)
                    cpu_pos = read_positions(cpu_dir / "last.dat")
                    gpu_pos = read_positions(gpu_dir / "last.dat")
                    self.assertEqual(len(cpu_pos), len(gpu_pos))
                    max_dev = max(abs(a - b) for p, q in zip(cpu_pos, gpu_pos) for a, b in zip(p, q))
                    self.assertLess(max_dev, POSITION_TOLERANCE[PRECISION], "largest position difference: %g" % max_dev)

    @unittest.skipIf(PRECISION == "double", "use_edge is not compatible with double precision")
    def test_edge_forces_match_matrix_forces(self):
        """The edge-based kernels give the same per-nucleotide forces, torques and energies as the neighbour-matrix one"""
        dumps = {}
        for mode, extra in CUDA_MODES:
            d = Path(self.tmp) / mode
            d.mkdir()
            dump = d / "dump.txt"
            r = run("mixed", "CUDA", d, extra=extra + "CUDA_dump_file = %s\nCUDA_sort_every = 0\n" % dump)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            rows = [l.split() for l in dump.read_text().splitlines() if l.strip() and not l.startswith("#")]
            rows.sort(key=lambda row: int(row[0]))
            dumps[mode] = [[float(x) for x in row[1:]] for row in rows]
        reference = dumps["matrix"]
        # columns: force (3), energy, torque (3), hydrogen-bonding energy
        scales = [max(abs(row[c]) for row in reference) for c in range(8)]
        for mode, values in dumps.items():
            with self.subTest(mode=mode):
                self.assertEqual(len(values), len(reference))
                for c, name in enumerate(["fx", "fy", "fz", "energy", "tx", "ty", "tz", "hb_energy"]):
                    dev = max(abs(a[c] - b[c]) for a, b in zip(values, reference))
                    self.assertLess(dev, EDGE_TOLERANCE * scales[c], "%s: largest deviation %g (scale %g)" % (name, dev, scales[c]))

    def test_external_force(self):
        """The CUDA external forces act on the NA nucleotides as the CPU ones do"""
        x, y, z = read_positions(HERE / "hybrid_overhang" / "init.dat")[0]
        forces = Path(self.tmp) / "forces.dat"
        forces.write_text("{\ntype = trap\nparticle = 0\npos0 = %g, %g, %g\nstiff = 5.\nrate = 0.\ndir = 1., 0., 0.\n}\n" % (x + 0.8, y, z))
        extra = "external_forces = true\nexternal_forces_file = %s\n" % forces
        steps, dt = 300, 0.001
        final = {}
        for backend in ("CPU", "CUDA"):
            d = Path(self.tmp) / backend
            d.mkdir()
            r = run("hybrid_overhang", backend, d, steps=steps, dt=dt, every=steps, extra=extra)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            final[backend] = read_positions(d / "last.dat")
        max_dev = max(abs(a - b) for p, q in zip(final["CPU"], final["CUDA"]) for a, b in zip(p, q))
        self.assertLess(max_dev, POSITION_TOLERANCE[PRECISION], "largest position difference: %g" % max_dev)

    def test_unsupported_options_fail(self):
        """Options that the CUDA NA interaction does not support make oxDNA exit with an error"""
        for extra, message in (("CUDA_update_stress_tensor_every = 1\n", "stress tensor"),):
            with self.subTest(option=extra.strip()):
                r = run("dna_duplex", "CUDA", self.tmp, extra=extra)
                self.assertNotEqual(r.returncode, 0)
                log = (Path(self.tmp) / "log.dat").read_text() if (Path(self.tmp) / "log.dat").exists() else ""
                self.assertIn(message, r.stdout + r.stderr + log)


if __name__ == "__main__":
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        OXDNA = os.path.abspath(sys.argv.pop(1))
    unittest.main()
