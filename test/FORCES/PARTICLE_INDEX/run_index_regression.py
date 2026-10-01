#!/usr/bin/env python3
"""Regression checks for single-particle indexes in external force files.

The stock TestSuite treats a non-zero exit and any log line starting with
ERROR as a failed simulation, then skips compare checks. Expected rejections
therefore cannot be ordinary quick_input cases. This driver runs oxDNA itself.
`make test_quick` and `make test_run` reach it through the ParticleIndex
entry in quick_compare and run_compare.

A rejection counts only when the process exits non-zero, is not a crash or
timeout, and the combined stdout, stderr, and log contain every required
diagnostic. A crash, timeout, or unrelated initialization error is a failure.

Positive cases use the 15-particle SSDNA fixture, so particle 0 and "last"
(14) are different indexes. Trap checks require the reference index and the
attachment index on the same log line.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOPOLOGY = REPO / "test" / "DNA" / "SSDNA15" / "MD" / "ssdna15.top"
CONFIGURATION = REPO / "test" / "DNA" / "SSDNA15" / "MD" / "init.dat"
N_PARTICLES = 15
LAST = N_PARTICLES - 1
TIMEOUT_SECONDS = 60

# Windows NTSTATUS values, Unix signal exits, and abort codes. oxDNA reports
# a handled oxDNAException as exit code 1, which is not in this set.
CRASH_CODES = {
    3,
    134,
    136,
    139,
    -6,
    -11,
    3221225477,  # STATUS_ACCESS_VIOLATION
    -1073741819,
    3221225725,  # STATUS_STACK_OVERFLOW
    -1073741571,
    3221226505,  # STATUS_STACK_BUFFER_OVERRUN
    -1073740791,
}


def is_crash(returncode):
    if returncode is None or returncode < 0:
        return True
    if returncode in CRASH_CODES:
        return True
    if returncode >= 0xC0000000:
        return True
    return False


def input_text(force_name):
    return f"""backend = CPU
sim_type = MC
ensemble = NVT
steps = 1
seed = 42
delta_translation = 0.05
delta_rotation = 0.10
T = 0.1
verlet_skin = 0.5
topology = ssdna15.top
conf_file = init.dat
trajectory_file = trajectory.dat
energy_file = energy.dat
log_file = run.log
no_stdout_energy = 1
restart_step_counter = 1
print_conf_interval = 100
print_energy_every = 100
time_scale = linear
external_forces = 1
external_forces_file = {force_name}
"""


def adding_line(force_label, ref_particle, particle):
    """Match one log line that records both indexes, when a reference exists."""
    if ref_particle is None:
        return re.compile(
            rf"Adding a {re.escape(force_label)} \(.* on particle {particle}\b"
        )
    return re.compile(
        rf"Adding a {re.escape(force_label)} \(.*ref_particle={ref_particle}\b.* on particle {particle}\b"
    )


CASES = [
    {
        "name": "mutual_numeric",
        "expect": "accept",
        "force": "type = mutual_trap\nparticle = 1\nref_particle = 2\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("MutualTrap", 2, 1),
    },
    {
        "name": "mutual_ref_last",
        "expect": "accept",
        "force": "type = mutual_trap\nparticle = 0\nref_particle = last\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("MutualTrap", LAST, 0),
    },
    {
        "name": "mutual_particle_last",
        "expect": "accept",
        "force": "type = mutual_trap\nparticle = last\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("MutualTrap", 0, LAST),
    },
    {
        "name": "constant_numeric",
        "expect": "accept",
        "force": "type = constant_trap\nparticle = 1\nref_particle = 2\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("ConstantTrap", 2, 1),
    },
    {
        "name": "constant_ref_last",
        "expect": "accept",
        "force": "type = constant_trap\nparticle = 0\nref_particle = last\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("ConstantTrap", LAST, 0),
    },
    {
        "name": "constant_particle_last",
        "expect": "accept",
        "force": "type = constant_trap\nparticle = last\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "pattern": adding_line("ConstantTrap", 0, LAST),
    },
    {
        "name": "alignment_numeric",
        "expect": "accept",
        "force": "type = alignment_field\nparticle = 0\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "pattern": adding_line("AlignmentField", None, 0),
    },
    {
        "name": "alignment_last",
        "expect": "accept",
        "force": "type = alignment_field\nparticle = last\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "pattern": adding_line("AlignmentField", None, LAST),
    },
    {
        "name": "mutual_particle_garbage",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = not_an_index\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap particle", "not_an_index"],
    },
    {
        "name": "mutual_ref_garbage",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = 0\nref_particle = not_an_index\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap ref_particle", "not_an_index"],
    },
    {
        "name": "mutual_particle_out_of_range",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = 99\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap particle", "non-existent particle 99"],
    },
    {
        "name": "mutual_ref_out_of_range",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = 0\nref_particle = 99\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap ref_particle", "non-existent particle 99"],
    },
    {
        "name": "mutual_particle_list",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = 0,1\nref_particle = 2\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap particle: expected exactly one particle"],
    },
    {
        "name": "mutual_ref_all",
        "expect": "reject",
        "force": "type = mutual_trap\nparticle = 0\nref_particle = all\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["MutualTrap ref_particle: expected exactly one particle", '"all"'],
    },
    {
        "name": "constant_particle_garbage",
        "expect": "reject",
        "force": "type = constant_trap\nparticle = not_an_index\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["ConstantTrap particle", "not_an_index"],
    },
    {
        "name": "constant_ref_list",
        "expect": "reject",
        "force": "type = constant_trap\nparticle = 0\nref_particle = 1,2\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["ConstantTrap ref_particle: expected exactly one particle"],
    },
    {
        "name": "constant_particle_out_of_range",
        "expect": "reject",
        "force": "type = constant_trap\nparticle = 99\nref_particle = 0\nstiff = 1.\nr0 = 5\n",
        "diagnostics": ["ConstantTrap particle", "non-existent particle 99"],
    },
    {
        "name": "alignment_garbage",
        "expect": "reject",
        "force": "type = alignment_field\nparticle = not_an_index\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "diagnostics": ["AlignmentField particle", "not_an_index"],
    },
    {
        "name": "alignment_out_of_range",
        "expect": "reject",
        "force": "type = alignment_field\nparticle = 99\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "diagnostics": ["AlignmentField particle", "non-existent particle 99"],
    },
    {
        "name": "alignment_list",
        "expect": "reject",
        "force": "type = alignment_field\nparticle = 0-2\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "diagnostics": ["AlignmentField particle: expected exactly one particle"],
    },
    {
        "name": "alignment_omitted_particle",
        "expect": "reject",
        "force": "type = alignment_field\nv_idx = 0\nF = 1.\ndir = 0,0,1\n",
        "diagnostics": ["Mandatory key `particle' not found"],
    },
]


def combined_output(stdout, stderr, log_path):
    parts = [stdout or "", stderr or ""]
    if log_path.is_file():
        parts.append(log_path.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


def run_oxdna(executable, workdir, timeout):
    log_path = workdir / "run.log"
    try:
        completed = subprocess.run(
            [str(executable), "input", "log_file=run.log"],
            cwd=workdir,
            timeout=timeout,
            capture_output=True,
            text=True,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = exc.stderr if isinstance(exc.stderr, str) else ""
        return {
            "status": "timeout",
            "returncode": None,
            "output": combined_output(stdout, stderr, log_path),
        }
    output = combined_output(completed.stdout, completed.stderr, log_path)
    if is_crash(completed.returncode):
        status = "crash"
    else:
        status = "exited"
    return {"status": status, "returncode": completed.returncode, "output": output}


def evaluate(case, run):
    output = run["output"]
    if case["expect"] == "accept":
        if run["status"] != "exited" or run["returncode"] != 0:
            return False, f"{run['status']} returncode={run['returncode']}"
        if case["pattern"].search(output) is None:
            return False, "accepted run did not log the expected particle and reference indexes"
        return True, "accepted with the expected indexes"
    if run["status"] == "timeout":
        return False, "timed out; not counted as a rejection"
    if run["status"] == "crash":
        return False, f"crash returncode={run['returncode']}; not counted as a rejection"
    if run["returncode"] == 0:
        return False, "exited 0 instead of rejecting the force file"
    missing = [text for text in case["diagnostics"] if text not in output]
    if missing:
        return False, "rejected without the expected diagnostic; missing " + ", ".join(missing)
    return True, f"rejected with returncode={run['returncode']} and the expected diagnostic"


def prepare_case(case, root):
    workdir = root / case["name"]
    workdir.mkdir()
    shutil.copy(TOPOLOGY, workdir / "ssdna15.top")
    shutil.copy(CONFIGURATION, workdir / "init.dat")
    (workdir / "input").write_text(input_text("forces.txt"), encoding="utf-8")
    (workdir / "forces.txt").write_text("{\n" + case["force"] + "}\n", encoding="utf-8")
    return workdir


def run_suite(executable, timeout):
    executable = Path(executable).resolve()
    if not executable.is_file():
        raise SystemExit(f"oxDNA executable not found: {executable}")
    if not TOPOLOGY.is_file() or not CONFIGURATION.is_file():
        raise SystemExit("SSDNA15 topology or configuration is missing")
    results = []
    with tempfile.TemporaryDirectory(prefix="oxdna-index-") as temporary:
        root = Path(temporary)
        for case in CASES:
            workdir = prepare_case(case, root)
            run = run_oxdna(executable, workdir, timeout)
            passed, detail = evaluate(case, run)
            results.append({
                "name": case["name"],
                "expect": case["expect"],
                "passed": passed,
                "detail": detail,
            })
    return results


def format_results(label, results):
    lines = [f"[{label}]"]
    for result in results:
        state = "PASS" if result["passed"] else "FAIL"
        lines.append(f"  {state}  {result['name']} ({result['expect']}): {result['detail']}")
    passed = sum(1 for result in results if result["passed"])
    lines.append(f"  {passed}/{len(results)} passed")
    return "\n".join(lines)


def compare_results(parent, patched):
    lines = ["[parent -> patched]"]
    by_parent = {result["name"]: result for result in parent}
    for result in patched:
        before = by_parent[result["name"]]
        if (not before["passed"]) and result["passed"]:
            change = "fail -> pass"
        elif before["passed"] and result["passed"]:
            change = "pass -> pass"
        elif (not before["passed"]) and (not result["passed"]):
            change = "fail -> fail"
        else:
            change = "pass -> fail"
        lines.append(f"  {change}  {result['name']}")
        if not before["passed"]:
            lines.append(f"    parent: {before['detail']}")
        if not result["passed"]:
            lines.append(f"    patched: {result['detail']}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", nargs="?", help="oxDNA binary to test")
    parser.add_argument("--parent", help="unpatched oxDNA binary")
    parser.add_argument("--patched", help="patched oxDNA binary")
    parser.add_argument("--label", default="oxDNA", help="label used when testing one binary")
    parser.add_argument("--timeout", type=int, default=TIMEOUT_SECONDS)
    args = parser.parse_args()

    if args.parent or args.patched:
        if not args.parent or not args.patched:
            parser.error("--parent and --patched are required together")
        parent = run_suite(args.parent, args.timeout)
        patched = run_suite(args.patched, args.timeout)
        print(format_results("parent", parent))
        print(format_results("patched", patched))
        print(compare_results(parent, patched))
        return 0 if all(result["passed"] for result in patched) else 1

    if not args.executable:
        parser.error("pass an oxDNA executable, or --parent and --patched")
    results = run_suite(args.executable, args.timeout)
    print(format_results(args.label, results))
    return 0 if all(result["passed"] for result in results) else 1


if __name__ == "__main__":
    sys.exit(main())
