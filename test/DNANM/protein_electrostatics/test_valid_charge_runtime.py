#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve()
REPO = HERE.parents[3]
KDPG = REPO / "examples" / "DNANM" / "KDPG"


def replace_setting(text, key, value):
    pattern = re.compile(
        rf"(?m)^[ \t]*{re.escape(key)}[ \t]*=.*$"
    )

    replacement = f"{key} = {value}"

    if pattern.search(text):
        return pattern.sub(replacement, text)

    return text.rstrip() + "\n" + replacement + "\n"


def make_charge_file(topology: Path, output: Path):
    raw = [
        line.strip()
        for line in topology.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    header = raw[0].split()

    npro = int(header[3])

    proteins = []

    for line in raw[1:]:
        fields = line.split()

        strand = int(fields[0])

        if strand < 0:
            proteins.append(
                (strand, fields[1].upper())
            )

    if len(proteins) != npro:
        raise RuntimeError(
            f"Topology declares {npro} proteins "
            f"but {len(proteins)} protein rows were found"
        )

    lines = [
        "# DNANM_ES_CHARGES 1",
        "# index_space = protein_local",
        "# index_base = 0",
        "# charge_model = tdd-synthetic",
        "# charge_source = upstream-port-test",
        f"# n_protein_beads = {npro}",
        f"# n_groups = {npro}",
        "# columns = protein_index topology_strand residue charge_e group_id",
    ]

    for local, (strand, residue) in enumerate(proteins):

        if local == 0:
            charge = 1.0
        elif local == 1:
            charge = -1.0
        else:
            charge = 0.0

        lines.append(
            f"{local} {strand} {residue} "
            f"{charge:.8f} {local}"
        )

    output.write_text(
        "\n".join(lines) + "\n"
    )


class ValidChargeRuntime(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.oxdna = Path(
            os.environ["OXDNA_BIN"]
        ).resolve()

        if not cls.oxdna.is_file():
            raise RuntimeError(
                f"Invalid OXDNA_BIN: {cls.oxdna}"
            )

    def run_case(self, input_name):
        with tempfile.TemporaryDirectory(
            prefix="oxdna_es_valid_"
        ) as td:

            case = Path(td) / "KDPG"

            shutil.copytree(
                KDPG,
                case
            )

            charge_file = case / "protein_charges.dat"

            make_charge_file(
                case / "generated.top",
                charge_file
            )

            inp = case / input_name
            text = inp.read_text()

            text = replace_setting(
                text,
                "steps",
                "1"
            )

            text = replace_setting(
                text,
                "seed",
                "27"
            )

            text += (
                "\n"
                "protein_dna_electrostatics = true\n"
                "protein_protein_electrostatics = true\n"
                "protein_charge_file = protein_charges.dat\n"
                "protein_charge_scale = 0.01\n"
                "exclude_intragroup_electrostatics = true\n"
            )

            inp.write_text(text)

            proc = subprocess.run(
                [str(self.oxdna), inp.name],
                cwd=case,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=120,
            )

            log_text = ""

            logfile = case / "log.dat"

            if logfile.exists():
                log_text = logfile.read_text(
                    errors="replace"
                )

            combined = (
                proc.stdout
                + "\n"
                + log_text
            )

            self.assertEqual(
                proc.returncode,
                0,
                msg=combined
            )

            self.assertIn(
                "Protein electrostatics enabled",
                combined
            )

            self.assertIn(
                "Loaded",
                combined
            )

            self.assertIn(
                "protein charges",
                combined
            )

    def test_dnanm_valid_v1_charge_file(self):
        self.run_case("input_dnanm")

    def test_dnact_valid_v1_charge_file(self):
        self.run_case("input_dnact")


if __name__ == "__main__":
    unittest.main(verbosity=2)
