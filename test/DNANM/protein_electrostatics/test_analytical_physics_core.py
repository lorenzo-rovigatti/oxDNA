#!/usr/bin/env python3

import math
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


OXDNA = Path(
    os.environ["OXDNA_BIN"]
).resolve()

if not OXDNA.is_file():
    raise RuntimeError(
        f"Invalid OXDNA_BIN: {OXDNA}"
    )


# -----------------------------------------------------------------
# oxDNA2 Debye-Huckel reference parameters.
#
# These are the upstream defaults at:
#   T = 300 K = 0.1 simulation units
#   salt = 1 M
# -----------------------------------------------------------------

T_SIM = 0.1
SALT = 1.0

DH_LAMBDA_FACTOR = 0.3616455
DH_STRENGTH = 0.0543

LAMBDA = (
    DH_LAMBDA_FACTOR
    * math.sqrt(T_SIM / 0.1)
    / math.sqrt(SALT)
)

KAPPA = 1.0 / LAMBDA

RHIGH = 3.0 * LAMBDA

# Same expression used by upstream DNA2Interaction.
B_SMOOTH = -(
    math.exp(-RHIGH / LAMBDA)
    * DH_STRENGTH
    * DH_STRENGTH
    * (RHIGH + LAMBDA)
    * (RHIGH + LAMBDA)
) / (
    -4.0
    * RHIGH
    * RHIGH
    * RHIGH
    * LAMBDA
    * LAMBDA
    * DH_STRENGTH
)

RCUT = (
    RHIGH
    * (
        DH_STRENGTH * RHIGH
        + 3.0 * DH_STRENGTH * LAMBDA
    )
    / (
        DH_STRENGTH
        * (RHIGH + LAMBDA)
    )
)

DNA_EFFECTIVE_CHARGE_E = 0.815

# DNA2 BACK interaction centre at identity orientation.
BACK = (-0.3400, 0.3408, 0.0)

DNA_POS = (5.0, 5.0, 5.0)

# This places protein local index 0 exactly 1.02 units
# from the DNA BACK site along +z.
P1 = (
    4.6600,
    5.3408,
    6.0200,
)

# Spectator.
P2 = (
    12.0,
    12.0,
    12.0,
)

# Protein local index 2 is 0.8 units from P1 along +x.
P3 = (
    5.4600,
    5.3408,
    6.0200,
)


TOP = """4 4 1 3 1
1 A -1 -1
-1 K -1 -1
-2 K -1 -1
-3 K -1 -1
"""


def vec_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vec_sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def vec_scale(a, s):
    return tuple(x * s for x in a)


def vec_norm(a):
    return math.sqrt(sum(x * x for x in a))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dh_energy_force_on_second(
    charge_product,
    rvec,
):
    """
    Return:
      U
      force acting on the second site/particle

    rvec points from the first interaction site to the second.
    """

    r = vec_norm(rvec)

    if r <= 0.0 or r >= RCUT:
        return 0.0, (0.0, 0.0, 0.0)

    rhat = vec_scale(
        rvec,
        1.0 / r,
    )

    if r < RHIGH:

        exponential = math.exp(
            -KAPPA * r
        )

        energy = (
            charge_product
            * exponential
            * DH_STRENGTH
            / r
        )

        force_scalar = (
            charge_product
            * DH_STRENGTH
            * exponential
            * (
                KAPPA / r
                + 1.0 / (r * r)
            )
        )

        force = vec_scale(
            rhat,
            force_scalar,
        )

    else:

        energy = (
            charge_product
            * B_SMOOTH
            * (r - RCUT) ** 2
        )

        force_scalar = (
            -charge_product
            * 2.0
            * B_SMOOTH
            * (r - RCUT)
        )

        force = vec_scale(
            rhat,
            force_scalar,
        )

    return energy, force


def charges_text(
    charges,
    groups=(0, 1, 2),
):
    n_groups = len(set(groups))

    lines = [
        "# DNANM_ES_CHARGES 1",
        "# index_space = protein_local",
        "# index_base = 0",
        "# charge_model = analytical-unit-test",
        "# charge_source = analytical-unit-test",
        "# n_protein_beads = 3",
        f"# n_groups = {n_groups}",
        (
            "# columns = protein_index "
            "topology_strand residue charge_e group_id"
        ),
    ]

    strands = (-1, -2, -3)

    for i, (
        strand,
        q,
        group,
    ) in enumerate(
        zip(
            strands,
            charges,
            groups,
        )
    ):
        lines.append(
            f"{i} {strand} K "
            f"{q:.12f} {group}"
        )

    return "\n".join(lines) + "\n"


def conf_text(
    p1=P1,
    p2=P2,
    p3=P3,
    dna_theta_x=0.0,
):
    """
    oxDNA configuration format stores a1 and a3.

    Rotation of DNA around +x:
       a1 = x
       a3 = R_x(theta) z
          = (0, -sin(theta), cos(theta))
    """

    s = math.sin(dna_theta_x)
    c = math.cos(dna_theta_x)

    dna_a1 = (
        1.0,
        0.0,
        0.0,
    )

    dna_a3 = (
        0.0,
        -s,
        c,
    )

    identity_a1 = (
        1.0,
        0.0,
        0.0,
    )

    identity_a3 = (
        0.0,
        0.0,
        1.0,
    )

    # Moment(a) are intentionally nonzero in the input format.
    # refresh_vel=true then generates a deterministic set from seed.
    momenta = (
        1.0e-12,
        -2.0e-12,
        3.0e-12,
        1.0e-12,
        -2.0e-12,
        3.0e-12,
    )

    def particle_line(
        pos,
        a1,
        a3,
    ):
        vals = (
            *pos,
            *a1,
            *a3,
            *momenta,
        )

        return " ".join(
            f"{x:.15g}"
            for x in vals
        )

    return (
        "t = 0\n"
        "b = 20 20 20\n"
        "E = 0 0 0\n"
        + particle_line(
            DNA_POS,
            dna_a1,
            dna_a3,
        )
        + "\n"
        + particle_line(
            p1,
            identity_a1,
            identity_a3,
        )
        + "\n"
        + particle_line(
            p2,
            identity_a1,
            identity_a3,
        )
        + "\n"
        + particle_line(
            p3,
            identity_a1,
            identity_a3,
        )
        + "\n"
    )


def input_text(
    interaction_type,
    pd,
    pp,
    pair,
    *,
    half_charged_ends=False,
    protein_charge_scale=1.0,
    pd_scale=1.0,
    pp_scale=1.0,
    exclude_groups=True,
    parfile="none",
):
    extra_angular = ""

    if interaction_type == "DNACT":
        extra_angular = (
            "bending_k = 0.0\n"
            "torsion_k = 0.0\n"
        )

    return f"""backend = CPU
interaction_type = {interaction_type}

sim_type = MD
T = 300K
dt = 0.002

thermostat = no
refresh_vel = 1
seed = 314159

restart_step_counter = 1
time_scale = linear

verlet_skin = 0.25
ensemble = NVT

external_forces = 0
max_io = 100

salt_concentration = 1.0
dh_half_charged_ends = {"true" if half_charged_ends else "false"}

{extra_angular}parfile = {parfile}

topology = system.top
conf_file = init.dat

trajectory_file = trajectory.dat
lastconf_file = last_conf.dat
energy_file = energy.dat
log_file = log.dat

no_stdout_energy = 1

steps = 1
print_conf_interval = 1
print_energy_every = 1

protein_dna_electrostatics = {"true" if pd else "false"}
protein_protein_electrostatics = {"true" if pp else "false"}

protein_charge_file = protein_charges.dat

protein_charge_scale = {protein_charge_scale:.15g}
protein_dna_dh_scale = {pd_scale:.15g}
protein_protein_dh_scale = {pp_scale:.15g}

dna_effective_charge_e = 0.815

exclude_intragroup_electrostatics = {"true" if exclude_groups else "false"}

data_output_1 = {{
 name = pair_energy.dat
 print_every = 1

 col_1 = {{
  type = pair_energy
   print_header = true
 }}
}}

data_output_2 = {{
 name = pair_force.dat
 print_every = 1

 col_1 = {{
  type = pair_force
  particle_id = {pair[0]}
 }}
}}
"""


def first_numeric_row(path):
    for raw in path.read_text().splitlines():
        line = raw.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        try:
            return [
                float(x)
                for x in line.split()
            ]
        except ValueError:
            continue

    raise RuntimeError(
        f"No numeric row found in {path}"
    )


def first_block_pair_energy(path, pair):
    """Read selected pair from all-pairs mode using actual map keys.

    Original upstream PairEnergy prints each interacting pair once, with
    a final total. Stop at the first block boundary to preserve t=0.
    Return zero only when that pair is absent in the first frame.
    """
    selected = frozenset(pair)
    started = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("#id1 id2"):
            if started:
                break
            started = True
            continue
        if not started:
            continue
        if line.startswith("#Total energy"):
            break
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 3:
            raise RuntimeError(f"Malformed pair energy row: {line}")
        try:
            a, b = int(fields[0]), int(fields[1])
            total = float(fields[-1])
        except ValueError as exc:
            raise RuntimeError(f"Invalid pair energy row: {line}") from exc
        if frozenset((a,b)) == selected:
            if not math.isfinite(total):
                raise RuntimeError(f"Non-finite pair energy: {line}")
            return total
    if not started:
        raise RuntimeError(f"No first PairEnergy block header in {path}")
    return 0.0


def parse_pair_force(
    path,
    pair,
):
    """
    Return force/torque on each member of the selected pair.

    Result:
      {
        particle_id: {
          "force": (fx,fy,fz),
          "torque": (tx,ty,tz),
        }
      }

    If PairForce prints no selected interaction, return zeros.
    """

    wanted = set(pair)

    found = {}

    if path.exists():
        for raw in path.read_text().splitlines():

            line = raw.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) != 8:
                continue

            try:
                i = int(fields[0])
                j = int(fields[1])

                vals = [
                    float(x)
                    for x in fields[2:]
                ]

            except ValueError:
                continue

            if {i, j} != wanted:
                continue

            # Keep only the first output block (t = 0).
            if i not in found:
                found[i] = {
                    "force": tuple(vals[0:3]),
                    "torque": tuple(vals[3:6]),
                }

            if all(
                pid in found
                for pid in pair
            ):
                break

    for pid in pair:
        found.setdefault(
            pid,
            {
                "force": (
                    0.0,
                    0.0,
                    0.0,
                ),
                "torque": (
                    0.0,
                    0.0,
                    0.0,
                ),
            },
        )

    return found


def run_case(
    root,
    name,
    *,
    interaction_type="DNANM",
    pd=False,
    pp=False,
    pair=(0, 1),
    charges=(0.0, 0.0, 0.0),
    groups=(0, 1, 2),
    half_charged_ends=False,
    protein_charge_scale=1.0,
    pd_scale=1.0,
    pp_scale=1.0,
    exclude_groups=True,
    p1=P1,
    p2=P2,
    p3=P3,
    dna_theta_x=0.0,
    bond_par=None,
):
    case = root / name
    case.mkdir()

    (case / "system.top").write_text(
        TOP
    )

    (case / "init.dat").write_text(
        conf_text(
            p1=p1,
            p2=p2,
            p3=p3,
            dna_theta_x=dna_theta_x,
        )
    )

    (case / "protein_charges.dat").write_text(
        charges_text(
            charges,
            groups,
        )
    )

    parfile = "none"

    if bond_par is not None:
        parfile = "anm.par"

        (case / "anm.par").write_text(
            bond_par
        )
    elif interaction_type == "DNACT":
        # Explicit no-spring file preserves the former none fixture mechanics.
        parfile = "anm.par"
        n_protein = int(TOP.splitlines()[0].split()[3])
        (case / "anm.par").write_text(f"{n_protein}\n")

    (case / "input").write_text(
        input_text(
            interaction_type,
            pd,
            pp,
            pair,
            half_charged_ends=half_charged_ends,
            protein_charge_scale=protein_charge_scale,
            pd_scale=pd_scale,
            pp_scale=pp_scale,
            exclude_groups=exclude_groups,
            parfile=parfile,
        )
    )

    proc = subprocess.run(
        [
            str(OXDNA),
            "input",
        ],
        cwd=case,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
    )

    log_text = ""

    logfile = case / "log.dat"

    if logfile.exists():
        log_text = logfile.read_text(
            errors="replace"
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"{name} failed\n\n"
            f"STDOUT/STDERR:\n{proc.stdout}\n\n"
            f"LOG:\n{log_text}"
        )

    pair_energy = first_block_pair_energy(
        case / "pair_energy.dat", pair
    )

    force_info = parse_pair_force(
        case / "pair_force.dat",
        pair,
    )

    return {
        "energy": pair_energy,
        "pair": pair,
        "force_info": force_info,
        "case": case,
    }


class AnalyticalPhysicsCore(unittest.TestCase):

    def assertScalarClose(
        self,
        observed,
        expected,
        *,
        atol=2.0e-6,
        rtol=5.0e-5,
        label="scalar",
    ):
        tol = max(
            atol,
            rtol * abs(expected),
        )

        self.assertAlmostEqual(
            observed,
            expected,
            delta=tol,
            msg=(
                f"{label}: "
                f"observed={observed:.12e}, "
                f"expected={expected:.12e}, "
                f"tol={tol:.3e}"
            ),
        )

    def assertVecClose(
        self,
        observed,
        expected,
        *,
        atol=3.0e-6,
        rtol=1.0e-4,
        label="vector",
    ):
        for axis, (
            obs,
            exp,
        ) in enumerate(
            zip(
                observed,
                expected,
            )
        ):
            self.assertScalarClose(
                obs,
                exp,
                atol=atol,
                rtol=rtol,
                label=f"{label}[{axis}]",
            )

    # -------------------------------------------------------------
    # PD analytic energy, force, Newton III, BACK torque
    # -------------------------------------------------------------

    def test_01_pd_analytic_energy_force_torque(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pd_analytic_"
        ) as td:

            root = Path(td)

            result = run_case(
                root,
                "pd_positive",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            # Vector from DNA BACK site to protein centre.
            rback = (
                0.0,
                0.0,
                1.02,
            )

            charge_product = (
                -(+1.0 / DNA_EFFECTIVE_CHARGE_E)
            )

            expected_E, expected_Fprotein = (
                dh_energy_force_on_second(
                    charge_product,
                    rback,
                )
            )

            self.assertLess(
                expected_E,
                0.0,
            )

            self.assertScalarClose(
                result["energy"],
                expected_E,
                label="PD energy",
            )

            Fdna = result[
                "force_info"
            ][0]["force"]

            Fprotein = result[
                "force_info"
            ][1]["force"]

            self.assertVecClose(
                Fprotein,
                expected_Fprotein,
                label="PD protein force",
            )

            self.assertVecClose(
                Fdna,
                vec_scale(
                    expected_Fprotein,
                    -1.0,
                ),
                label="PD DNA force",
            )

            self.assertVecClose(
                vec_add(
                    Fdna,
                    Fprotein,
                ),
                (
                    0.0,
                    0.0,
                    0.0,
                ),
                atol=4.0e-6,
                label="PD Newton III",
            )

            expected_tau_dna = cross(
                BACK,
                Fdna,
            )

            tau_dna = result[
                "force_info"
            ][0]["torque"]

            self.assertVecClose(
                tau_dna,
                expected_tau_dna,
                atol=4.0e-6,
                label="PD BACK torque",
            )

            # Point charge is at protein COM:
            # electrostatic term should not torque the protein.
            self.assertVecClose(
                result[
                    "force_info"
                ][1]["torque"],
                (
                    0.0,
                    0.0,
                    0.0,
                ),
                atol=4.0e-6,
                label="PD protein torque",
            )

            print(
                "PD_ANALYTIC",
                result["energy"],
                expected_E,
                Fprotein,
                tau_dna,
            )

    # -------------------------------------------------------------
    # PD sign reversal
    # -------------------------------------------------------------

    def test_02_pd_charge_sign_controls_attraction_repulsion(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pd_sign_"
        ) as td:

            root = Path(td)

            attractive = run_case(
                root,
                "positive_protein",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            repulsive = run_case(
                root,
                "negative_protein",
                pd=True,
                pair=(0, 1),
                charges=(
                    -1.0,
                    0.0,
                    0.0,
                ),
            )

            self.assertLess(
                attractive["energy"],
                0.0,
            )

            self.assertGreater(
                repulsive["energy"],
                0.0,
            )

            self.assertLess(
                attractive[
                    "force_info"
                ][1]["force"][2],
                0.0,
            )

            self.assertGreater(
                repulsive[
                    "force_info"
                ][1]["force"][2],
                0.0,
            )

            self.assertScalarClose(
                repulsive["energy"],
                -attractive["energy"],
                atol=3.0e-6,
                label="PD charge sign energy",
            )

    # -------------------------------------------------------------
    # PP same-charge analytic repulsion
    # -------------------------------------------------------------

    def test_03_pp_analytic_same_charge_repulsion(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pp_analytic_"
        ) as td:

            root = Path(td)

            result = run_case(
                root,
                "pp_like",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            qred = (
                1.0
                /
                DNA_EFFECTIVE_CHARGE_E
            )

            charge_product = (
                qred * qred
            )

            expected_E, expected_F3 = (
                dh_energy_force_on_second(
                    charge_product,
                    (
                        0.8,
                        0.0,
                        0.0,
                    ),
                )
            )

            self.assertGreater(
                expected_E,
                0.0,
            )

            self.assertScalarClose(
                result["energy"],
                expected_E,
                label="PP energy",
            )

            F1 = result[
                "force_info"
            ][1]["force"]

            F3 = result[
                "force_info"
            ][3]["force"]

            self.assertVecClose(
                F3,
                expected_F3,
                label="PP force second",
            )

            self.assertVecClose(
                F1,
                vec_scale(
                    expected_F3,
                    -1.0,
                ),
                label="PP force first",
            )

            self.assertVecClose(
                vec_add(
                    F1,
                    F3,
                ),
                (
                    0.0,
                    0.0,
                    0.0,
                ),
                atol=4.0e-6,
                label="PP Newton III",
            )

            self.assertGreater(
                F3[0],
                0.0,
            )

            print(
                "PP_ANALYTIC",
                result["energy"],
                expected_E,
                F3,
            )

    # -------------------------------------------------------------
    # PP opposite-charge attraction
    # -------------------------------------------------------------

    def test_04_pp_opposite_charge_attraction(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pp_opposite_"
        ) as td:

            root = Path(td)

            result = run_case(
                root,
                "pp_opposite",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    -1.0,
                ),
            )

            self.assertLess(
                result["energy"],
                0.0,
            )

            self.assertLess(
                result[
                    "force_info"
                ][3]["force"][0],
                0.0,
            )

    # -------------------------------------------------------------
    # DNANM and DNACT must share exactly the same ES physics here.
    # -------------------------------------------------------------

    def test_05_dnanm_dnact_shared_pd_physics(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_shared_pd_"
        ) as td:

            root = Path(td)

            a = run_case(
                root,
                "dnanm",
                interaction_type="DNANM",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            b = run_case(
                root,
                "dnact",
                interaction_type="DNACT",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            self.assertScalarClose(
                b["energy"],
                a["energy"],
                atol=2.0e-7,
                rtol=1.0e-7,
                label="DNANM/DNACT PD energy",
            )

            self.assertVecClose(
                b["force_info"][1]["force"],
                a["force_info"][1]["force"],
                atol=2.0e-7,
                rtol=1.0e-7,
                label="DNANM/DNACT PD force",
            )

    def test_06_dnanm_dnact_shared_pp_physics(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_shared_pp_"
        ) as td:

            root = Path(td)

            a = run_case(
                root,
                "dnanm",
                interaction_type="DNANM",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            b = run_case(
                root,
                "dnact",
                interaction_type="DNACT",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            self.assertScalarClose(
                b["energy"],
                a["energy"],
                atol=2.0e-7,
                rtol=1.0e-7,
                label="DNANM/DNACT PP energy",
            )

            self.assertVecClose(
                b["force_info"][3]["force"],
                a["force_info"][3]["force"],
                atol=2.0e-7,
                rtol=1.0e-7,
                label="DNANM/DNACT PP force",
            )

    # -------------------------------------------------------------
    # oxDNA2 terminal half-charge convention
    # -------------------------------------------------------------

    def test_07_pd_half_charged_terminal_nucleotide(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_half_charge_"
        ) as td:

            root = Path(td)

            full = run_case(
                root,
                "full",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                half_charged_ends=False,
            )

            half = run_case(
                root,
                "half",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                half_charged_ends=True,
            )

            self.assertScalarClose(
                half["energy"],
                0.5 * full["energy"],
                atol=2.0e-6,
                label="PD terminal half energy",
            )

            self.assertVecClose(
                half["force_info"][1]["force"],
                vec_scale(
                    full["force_info"][1]["force"],
                    0.5,
                ),
                atol=3.0e-6,
                label="PD terminal half force",
            )

    # -------------------------------------------------------------
    # Explicit scaling semantics
    # -------------------------------------------------------------

    def test_08_explicit_pd_pp_scaling(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_scaling_"
        ) as td:

            root = Path(td)

            pd_base = run_case(
                root,
                "pd_base",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            pd_half = run_case(
                root,
                "pd_half_scale",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                pd_scale=0.5,
            )

            self.assertScalarClose(
                pd_half["energy"],
                0.5 * pd_base["energy"],
                atol=2.0e-6,
                label="PD DH scale",
            )

            pp_base = run_case(
                root,
                "pp_base",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            pp_half = run_case(
                root,
                "pp_half_scale",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                pp_scale=0.5,
            )

            self.assertScalarClose(
                pp_half["energy"],
                0.5 * pp_base["energy"],
                atol=2.0e-6,
                label="PP DH scale",
            )

            pp_charge_half = run_case(
                root,
                "pp_charge_half",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                protein_charge_scale=0.5,
            )

            # Both protein charges are scaled:
            # PP energy therefore scales quadratically.
            self.assertScalarClose(
                pp_charge_half["energy"],
                0.25 * pp_base["energy"],
                atol=2.0e-6,
                label="PP protein charge scale squared",
            )

    # -------------------------------------------------------------
    # Group exclusion switch
    # -------------------------------------------------------------

    def test_09_intragroup_exclusion_toggle(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_groups_"
        ) as td:

            root = Path(td)

            excluded = run_case(
                root,
                "excluded",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                groups=(
                    7,
                    8,
                    7,
                ),
                exclude_groups=True,
            )

            active = run_case(
                root,
                "active",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                groups=(
                    7,
                    8,
                    7,
                ),
                exclude_groups=False,
            )

            self.assertScalarClose(
                excluded["energy"],
                0.0,
                atol=1.0e-10,
                rtol=0.0,
                label="same-group exclusion",
            )

            self.assertGreater(
                active["energy"],
                0.0,
            )

    # -------------------------------------------------------------
    # ANM-connected PP pair must not receive PP-DH.
    # Two springs are used because current upstream DNANM parser
    # treats a single spring specially for N > 2.
    # -------------------------------------------------------------

    def test_10_anm_connected_pair_excluded_from_pp_dh(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_anm_exclusion_"
        ) as td:

            root = Path(td)

            nonbonded = run_case(
                root,
                "nonbonded",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            bonds = """3
1 3 0.8 s 0.0
1 2 12.2065556157 s 0.0
"""

            bonded = run_case(
                root,
                "bonded",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                bond_par=bonds,
            )

            self.assertGreater(
                nonbonded["energy"],
                0.0,
            )

            self.assertScalarClose(
                bonded["energy"],
                0.0,
                atol=2.0e-9,
                rtol=0.0,
                label="ANM bonded PP exclusion",
            )

    # -------------------------------------------------------------
    # Translational finite differences
    # -------------------------------------------------------------

    def test_11_pd_force_matches_energy_finite_difference(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pd_fd_"
        ) as td:

            root = Path(td)

            h = 1.0e-3

            base = run_case(
                root,
                "base",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
            )

            plus_p1 = (
                P1[0],
                P1[1],
                P1[2] + h,
            )

            minus_p1 = (
                P1[0],
                P1[1],
                P1[2] - h,
            )

            plus = run_case(
                root,
                "plus",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                p1=plus_p1,
            )

            minus = run_case(
                root,
                "minus",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                p1=minus_p1,
            )

            numerical_force_z = -(
                plus["energy"]
                -
                minus["energy"]
            ) / (
                2.0 * h
            )

            observed_force_z = (
                base[
                    "force_info"
                ][1]["force"][2]
            )

            self.assertScalarClose(
                observed_force_z,
                numerical_force_z,
                atol=4.0e-5,
                rtol=4.0e-3,
                label="PD translational finite difference",
            )

    def test_12_pp_force_matches_energy_finite_difference(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pp_fd_"
        ) as td:

            root = Path(td)

            h = 1.0e-3

            base = run_case(
                root,
                "base",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
            )

            plus_p3 = (
                P3[0] + h,
                P3[1],
                P3[2],
            )

            minus_p3 = (
                P3[0] - h,
                P3[1],
                P3[2],
            )

            plus = run_case(
                root,
                "plus",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                p3=plus_p3,
            )

            minus = run_case(
                root,
                "minus",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                p3=minus_p3,
            )

            numerical_force_x = -(
                plus["energy"]
                -
                minus["energy"]
            ) / (
                2.0 * h
            )

            observed_force_x = (
                base[
                    "force_info"
                ][3]["force"][0]
            )

            self.assertScalarClose(
                observed_force_x,
                numerical_force_x,
                atol=5.0e-5,
                rtol=4.0e-3,
                label="PP translational finite difference",
            )

    # -------------------------------------------------------------
    # Rotational finite difference of the DNA BACK torque.
    # Rotation axis is DNA a1 == laboratory +x at theta=0.
    # -------------------------------------------------------------

    def test_13_pd_back_torque_matches_rotational_finite_difference(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_pd_torque_fd_"
        ) as td:

            root = Path(td)

            h = 1.0e-3

            base = run_case(
                root,
                "base",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                dna_theta_x=0.0,
            )

            plus = run_case(
                root,
                "plus",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                dna_theta_x=+h,
            )

            minus = run_case(
                root,
                "minus",
                pd=True,
                pair=(0, 1),
                charges=(
                    +1.0,
                    0.0,
                    0.0,
                ),
                dna_theta_x=-h,
            )

            numerical_tau_x = -(
                plus["energy"]
                -
                minus["energy"]
            ) / (
                2.0 * h
            )

            observed_tau_x = (
                base[
                    "force_info"
                ][0]["torque"][0]
            )

            self.assertScalarClose(
                observed_tau_x,
                numerical_tau_x,
                atol=5.0e-5,
                rtol=6.0e-3,
                label="PD rotational finite difference",
            )

    # -------------------------------------------------------------
    # Smoothed branch + cutoff.
    # This directly tests reuse of the oxDNA2 DH tail.
    # -------------------------------------------------------------

    def test_14_smoothed_branch_and_cutoff(self):

        with tempfile.TemporaryDirectory(
            prefix="oxdna_dh_tail_"
        ) as td:

            root = Path(td)

            smooth_r = (
                RHIGH + RCUT
            ) * 0.5

            smooth_p3 = (
                P1[0] + smooth_r,
                P1[1],
                P1[2],
            )

            smooth = run_case(
                root,
                "smooth",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                p3=smooth_p3,
            )

            qred = (
                1.0
                /
                DNA_EFFECTIVE_CHARGE_E
            )

            expected_E, expected_F3 = (
                dh_energy_force_on_second(
                    qred * qred,
                    (
                        smooth_r,
                        0.0,
                        0.0,
                    ),
                )
            )

            self.assertScalarClose(
                smooth["energy"],
                expected_E,
                atol=2.0e-6,
                label="DH smoothed energy",
            )

            self.assertVecClose(
                smooth["force_info"][3]["force"],
                expected_F3,
                atol=4.0e-6,
                label="DH smoothed force",
            )

            outside_r = (
                RCUT + 0.05
            )

            outside_p3 = (
                P1[0] + outside_r,
                P1[1],
                P1[2],
            )

            outside = run_case(
                root,
                "outside",
                pp=True,
                pair=(1, 3),
                charges=(
                    +1.0,
                    0.0,
                    +1.0,
                ),
                p3=outside_p3,
            )

            self.assertScalarClose(
                outside["energy"],
                0.0,
                atol=1.0e-10,
                rtol=0.0,
                label="DH zero beyond cutoff",
            )


if __name__ == "__main__":

    print(
        "ANALYTICAL_DH_CONSTANTS",
        f"lambda={LAMBDA:.12e}",
        f"RHIGH={RHIGH:.12e}",
        f"RCUT={RCUT:.12e}",
        f"B={B_SMOOTH:.12e}",
    )

    unittest.main(
        verbosity=2
    )
