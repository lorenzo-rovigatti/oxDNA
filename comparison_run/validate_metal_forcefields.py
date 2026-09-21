#!/usr/bin/env python3
"""
Validate CPU-vs-Metal parity for core forcefield families.

This runner executes short deterministic MD jobs for:
  - DNA
  - DNA2
  - LJ
  - RNA
  - RNA2
  - patchy
  - TEP

and checks that relative potential-energy error stays within tolerance.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CPU_BIN = ROOT / "build_cpu" / "bin" / "oxDNA"
DEFAULT_METAL_BIN = ROOT / "build_metal" / "bin" / "oxDNA"
DEFAULT_CONFGEN_BIN = ROOT / "build_metal" / "bin" / "confGenerator"
DEFAULT_SHADER_LIB = ROOT / "build_metal" / "bin" / "shaders.metallib"


class ValidationError(RuntimeError):
    pass


@dataclass
class Scenario:
    name: str
    interaction_type: str
    temperature: str
    dt: float
    verlet_skin: float
    steps: int
    extra_input: dict[str, str] = field(default_factory=dict)
    setup_fn: Callable[[Path, Path], None] | None = None


@dataclass
class ScenarioResult:
    name: str
    status: str
    max_rel_error: float | None = None
    step_index: int | None = None
    message: str = ""


def run_cmd(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True)


def write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def parse_energy_file(path: Path) -> list[tuple[float, float]]:
    if not path.exists():
        raise ValidationError(f"Missing energy file: {path}")

    rows: list[tuple[float, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            t = float(parts[0])
            pot = float(parts[1])
        except ValueError:
            continue
        if not math.isfinite(t) or not math.isfinite(pot):
            raise ValidationError(f"NaN/Inf in energy file: {path}")
        rows.append((t, pot))

    if not rows:
        raise ValidationError(f"Energy file has no parsable data: {path}")
    return rows


def ensure_no_error_markers(log_path: Path, stdouterr: str) -> None:
    chunks: list[str] = [stdouterr]
    if log_path.exists():
        chunks.append(log_path.read_text(encoding="utf-8", errors="replace"))
    full = "\n".join(chunks)
    if "ERROR:" in full:
        raise ValidationError(f"Detected ERROR marker in run output/log: {log_path}")


def write_input_file(
    path: Path,
    backend: str,
    scenario: Scenario,
    seed: int,
    trajectory_name: str,
    energy_name: str,
    log_name: str,
    metal_avoid_cpu_calculations: int,
) -> None:
    lines: list[str] = [
        f"backend = {backend}",
        "sim_type = MD",
        f"interaction_type = {scenario.interaction_type}",
        f"steps = {scenario.steps}",
        "newtonian_steps = 0",
        "thermostat = no",
        f"T = {scenario.temperature}",
        f"dt = {scenario.dt}",
        f"verlet_skin = {scenario.verlet_skin}",
        "topology = topology.top",
        "conf_file = init.dat",
        f"trajectory_file = {trajectory_name}",
        f"energy_file = {energy_name}",
        f"log_file = {log_name}",
        "print_conf_interval = 1000000",
        "print_energy_every = 1",
        "refresh_vel = 1",
        f"seed = {seed}",
        "restart_step_counter = 1",
        "external_forces = 0",
        "no_stdout_energy = 1",
        "time_scale = linear",
    ]

    if backend == "Metal":
        lines.append(f"Metal_avoid_cpu_calculations = {metal_avoid_cpu_calculations}")

    for key, value in scenario.extra_input.items():
        lines.append(f"{key} = {value}")

    write_text(path, "\n".join(lines) + "\n")


def copy_fixture(src_top: Path, src_conf: Path, dst_dir: Path) -> None:
    shutil.copy2(src_top, dst_dir / "topology.top")
    shutil.copy2(src_conf, dst_dir / "init.dat")


def setup_dna(dst_dir: Path, _confgen_bin: Path) -> None:
    copy_fixture(
        ROOT / "test" / "METAL" / "SIMPLE_MD" / "ssdna15.top",
        ROOT / "test" / "METAL" / "SIMPLE_MD" / "init.dat",
        dst_dir,
    )


def setup_dna2(dst_dir: Path, _confgen_bin: Path) -> None:
    copy_fixture(
        ROOT / "test" / "DNA" / "FORCE_FIELD" / "init.top",
        ROOT / "test" / "DNA" / "FORCE_FIELD" / "init.dat",
        dst_dir,
    )


def setup_lj(dst_dir: Path, _confgen_bin: Path) -> None:
    copy_fixture(
        ROOT / "test" / "LJ" / "topology.dat",
        ROOT / "test" / "LJ" / "init_conf.dat",
        dst_dir,
    )


def setup_rna(dst_dir: Path, _confgen_bin: Path) -> None:
    copy_fixture(
        ROOT / "test" / "RNA" / "FORCE_FIELD" / "init.top",
        ROOT / "test" / "RNA" / "FORCE_FIELD" / "init.dat",
        dst_dir,
    )


def setup_patchy(dst_dir: Path, confgen_bin: Path) -> None:
    write_text(dst_dir / "topology.top", "16 0\n")
    write_text(
        dst_dir / "input_gen",
        "\n".join(
            [
                "backend = CPU",
                "interaction_type = patchy",
                "PATCHY_N = 4",
                "T = 1.0",
                "topology = topology.top",
                "trajectory_file = /dev/null",
                "conf_file = init.dat",
                "seed = 12345",
            ]
        )
        + "\n",
    )
    result = run_cmd([str(confgen_bin), "input_gen", "4"], dst_dir)
    if result.returncode != 0:
        raise ValidationError(f"patchy confGenerator failed:\n{result.stdout}\n{result.stderr}")
    if not (dst_dir / "init.dat").exists():
        raise ValidationError("patchy confGenerator did not produce init.dat")


def setup_tep(dst_dir: Path, confgen_bin: Path) -> None:
    write_text(dst_dir / "topology.top", "12 1\n12\n")
    write_text(
        dst_dir / "input_gen",
        "\n".join(
            [
                "backend = CPU",
                "interaction_type = TEP",
                "T = 1.0",
                "topology = topology.top",
                "trajectory_file = /dev/null",
                "conf_file = init.dat",
                "seed = 12345",
            ]
        )
        + "\n",
    )
    result = run_cmd([str(confgen_bin), "input_gen", "20"], dst_dir)
    if result.returncode != 0:
        raise ValidationError(f"TEP confGenerator failed:\n{result.stdout}\n{result.stderr}")
    if not (dst_dir / "init.dat").exists():
        raise ValidationError("TEP confGenerator did not produce init.dat")


def max_relative_error(
    cpu_rows: list[tuple[float, float]],
    metal_rows: list[tuple[float, float]],
) -> tuple[float, int]:
    n = min(len(cpu_rows), len(metal_rows))
    if n == 0:
        raise ValidationError("No overlap between CPU and Metal energy samples")

    max_rel = 0.0
    max_idx = 0
    for i in range(n):
        cpu_pot = cpu_rows[i][1]
        metal_pot = metal_rows[i][1]
        diff = abs(cpu_pot - metal_pot)

        if abs(cpu_pot) < 1e-12:
            rel = 0.0 if diff < 1e-12 else diff
        else:
            rel = diff / abs(cpu_pot)

        if rel > max_rel:
            max_rel = rel
            max_idx = i

    return max_rel, max_idx


def ensure_no_explosion(rows: list[tuple[float, float]], label: str) -> None:
    for _, pot in rows:
        if not math.isfinite(pot):
            raise ValidationError(f"{label}: potential contains NaN/Inf")
        if abs(pot) > 1e8:
            raise ValidationError(f"{label}: potential appears to explode (|U| > 1e8)")


def run_single_scenario(
    scenario: Scenario,
    base_tmp: Path,
    cpu_bin: Path,
    metal_bin: Path,
    confgen_bin: Path,
    shader_lib: Path,
    tol: float,
    seed: int,
    metal_avoid_cpu_calculations: int,
) -> ScenarioResult:
    scenario_dir = base_tmp / scenario.name
    scenario_dir.mkdir(parents=True, exist_ok=True)

    if scenario.setup_fn is None:
        return ScenarioResult(scenario.name, "FAIL", message="Internal error: missing setup_fn")

    scenario.setup_fn(scenario_dir, confgen_bin)

    cpu_dir = scenario_dir / "cpu"
    metal_dir = scenario_dir / "metal"
    cpu_dir.mkdir()
    metal_dir.mkdir()

    shutil.copy2(scenario_dir / "topology.top", cpu_dir / "topology.top")
    shutil.copy2(scenario_dir / "init.dat", cpu_dir / "init.dat")
    shutil.copy2(scenario_dir / "topology.top", metal_dir / "topology.top")
    shutil.copy2(scenario_dir / "init.dat", metal_dir / "init.dat")

    write_input_file(
        cpu_dir / "input.dat",
        "CPU",
        scenario,
        seed,
        "trajectory_cpu.dat",
        "energy_cpu.dat",
        "log_cpu.dat",
        metal_avoid_cpu_calculations=metal_avoid_cpu_calculations,
    )

    write_input_file(
        metal_dir / "input.dat",
        "Metal",
        scenario,
        seed,
        "trajectory_metal.dat",
        "energy_metal.dat",
        "log_metal.dat",
        metal_avoid_cpu_calculations=metal_avoid_cpu_calculations,
    )

    cpu_run = run_cmd([str(cpu_bin), "input.dat"], cpu_dir)
    if cpu_run.returncode != 0:
        return ScenarioResult(
            scenario.name,
            "FAIL",
            message=f"CPU run failed (exit {cpu_run.returncode})",
        )
    try:
        ensure_no_error_markers(cpu_dir / "log_cpu.dat", cpu_run.stdout + cpu_run.stderr)
    except ValidationError as exc:
        return ScenarioResult(scenario.name, "FAIL", message=str(exc))

    if shader_lib.exists():
        shutil.copy2(shader_lib, metal_dir / "shaders.metallib")

    metal_run = run_cmd([str(metal_bin), "input.dat"], metal_dir)
    if metal_run.returncode != 0:
        return ScenarioResult(
            scenario.name,
            "FAIL",
            message=f"Metal run failed (exit {metal_run.returncode})",
        )
    try:
        ensure_no_error_markers(metal_dir / "log_metal.dat", metal_run.stdout + metal_run.stderr)
    except ValidationError as exc:
        return ScenarioResult(scenario.name, "FAIL", message=str(exc))

    try:
        cpu_rows = parse_energy_file(cpu_dir / "energy_cpu.dat")
        metal_rows = parse_energy_file(metal_dir / "energy_metal.dat")
        ensure_no_explosion(cpu_rows, f"{scenario.name} CPU")
        ensure_no_explosion(metal_rows, f"{scenario.name} Metal")
        max_rel, idx = max_relative_error(cpu_rows, metal_rows)
    except ValidationError as exc:
        return ScenarioResult(scenario.name, "FAIL", message=str(exc))

    if max_rel > tol:
        return ScenarioResult(
            scenario.name,
            "FAIL",
            max_rel_error=max_rel,
            step_index=idx,
            message=f"relative error {max_rel:.6g} exceeds tolerance {tol:.6g}",
        )

    return ScenarioResult(
        scenario.name,
        "PASS",
        max_rel_error=max_rel,
        step_index=idx,
        message="ok",
    )


def build_scenarios(default_steps: int) -> list[Scenario]:
    return [
        Scenario(
            name="dna",
            interaction_type="DNA",
            temperature="300K",
            dt=0.001,
            verlet_skin=0.05,
            steps=default_steps,
            setup_fn=setup_dna,
        ),
        Scenario(
            name="dna2",
            interaction_type="DNA2",
            temperature="300K",
            dt=0.001,
            verlet_skin=0.05,
            steps=default_steps,
            extra_input={"salt_concentration": "1.0"},
            setup_fn=setup_dna2,
        ),
        Scenario(
            name="lj",
            interaction_type="LJ",
            temperature="1.5",
            dt=0.001,
            verlet_skin=0.2,
            steps=default_steps,
            setup_fn=setup_lj,
        ),
        Scenario(
            name="rna",
            interaction_type="RNA",
            temperature="300K",
            dt=0.001,
            verlet_skin=0.05,
            steps=default_steps,
            setup_fn=setup_rna,
        ),
        Scenario(
            name="rna2",
            interaction_type="RNA2",
            temperature="300K",
            dt=0.001,
            verlet_skin=0.05,
            steps=default_steps,
            extra_input={"salt_concentration": "1.0"},
            setup_fn=setup_rna,
        ),
        Scenario(
            name="patchy",
            interaction_type="patchy",
            temperature="1.0",
            dt=0.001,
            verlet_skin=0.2,
            steps=min(default_steps, 20),
            extra_input={"PATCHY_N": "4"},
            setup_fn=setup_patchy,
        ),
        Scenario(
            name="tep",
            interaction_type="TEP",
            temperature="1.0",
            dt=0.0005,
            verlet_skin=0.2,
            steps=default_steps,
            setup_fn=setup_tep,
        ),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate CPU-vs-Metal forcefield parity.")
    parser.add_argument("--cpu-bin", type=Path, default=DEFAULT_CPU_BIN)
    parser.add_argument("--metal-bin", type=Path, default=DEFAULT_METAL_BIN)
    parser.add_argument("--confgen-bin", type=Path, default=DEFAULT_CONFGEN_BIN)
    parser.add_argument("--shader-lib", type=Path, default=DEFAULT_SHADER_LIB)
    parser.add_argument("--tol", type=float, default=1e-3)
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--seed", type=int, default=129382)
    parser.add_argument("--scenarios", nargs="*", default=[])
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--tmp-dir", type=Path, default=None)
    parser.add_argument("--metal-avoid-cpu-calculations", type=int, choices=[0, 1], default=0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    for binary in [args.cpu_bin, args.metal_bin, args.confgen_bin]:
        if not binary.exists():
            print(f"ERROR: missing executable: {binary}", file=sys.stderr)
            return 2

    if args.metal_avoid_cpu_calculations != 0:
        print(
            "WARNING: Metal_avoid_cpu_calculations=1 is experimental. "
            "Parity failures are expected.",
            file=sys.stderr,
        )

    scenarios = build_scenarios(args.steps)
    if args.scenarios:
        wanted = {name.strip().lower() for name in args.scenarios}
        scenarios = [s for s in scenarios if s.name.lower() in wanted]
        if not scenarios:
            print("ERROR: no matching scenarios after --scenarios filter", file=sys.stderr)
            return 2

    if args.tmp_dir is None:
        tmp_root_obj = tempfile.TemporaryDirectory(prefix="validate_metal_forcefields_")
        tmp_root = Path(tmp_root_obj.name)
    else:
        tmp_root_obj = None
        tmp_root = args.tmp_dir
        tmp_root.mkdir(parents=True, exist_ok=True)

    print(f"Working directory: {tmp_root}")
    print(f"Tolerance: {args.tol:.3e}")
    print(f"Seed: {args.seed}")
    print("")

    results: list[ScenarioResult] = []
    for scenario in scenarios:
        print(f"[RUN] {scenario.name}")
        try:
            result = run_single_scenario(
                scenario=scenario,
                base_tmp=tmp_root,
                cpu_bin=args.cpu_bin,
                metal_bin=args.metal_bin,
                confgen_bin=args.confgen_bin,
                shader_lib=args.shader_lib,
                tol=args.tol,
                seed=args.seed,
                metal_avoid_cpu_calculations=args.metal_avoid_cpu_calculations,
            )
        except Exception as exc:  # noqa: BLE001
            result = ScenarioResult(scenario.name, "FAIL", message=f"Unhandled exception: {exc}")
        results.append(result)

        if result.status == "PASS":
            print(
                f"  PASS  max_rel={result.max_rel_error:.6g} "
                f"(sample_index={result.step_index})"
            )
        else:
            print(f"  FAIL  {result.message}")
        print("")

    n_pass = sum(1 for r in results if r.status == "PASS")
    n_fail = len(results) - n_pass

    print("Summary:")
    for r in results:
        rel = "-" if r.max_rel_error is None else f"{r.max_rel_error:.6g}"
        idx = "-" if r.step_index is None else str(r.step_index)
        print(f"  {r.name:8s} {r.status:4s} max_rel={rel:>10s} sample={idx:>4s}  {r.message}")

    print("")
    print(f"Passed: {n_pass}/{len(results)}")
    print(f"Failed: {n_fail}/{len(results)}")

    if args.keep_temp:
        print(f"Kept temporary data in: {tmp_root}")
    elif args.tmp_dir is not None:
        print(f"Temporary data available in --tmp-dir: {tmp_root}")

    if tmp_root_obj is not None and not args.keep_temp:
        tmp_root_obj.cleanup()

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
