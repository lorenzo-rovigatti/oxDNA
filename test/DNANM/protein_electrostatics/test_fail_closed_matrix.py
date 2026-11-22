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

OXDNA = Path(
    os.environ["OXDNA_BIN"]
).resolve()


def replace_setting(text, key, value):
    pattern = re.compile(
        rf"(?m)^[ \t]*{re.escape(key)}[ \t]*=.*$"
    )

    replacement = f"{key} = {value}"

    if pattern.search(text):
        return pattern.sub(
            replacement,
            text,
        )

    return (
        text.rstrip()
        + "\n"
        + replacement
        + "\n"
    )


def parse_proteins(topology):
    raw = [
        line.strip()
        for line in topology.read_text().splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    header = raw[0].split()

    npro = int(header[3])

    proteins = []

    for line in raw[1:]:

        fields = line.split()

        strand = int(fields[0])

        if strand < 0:
            proteins.append(
                {
                    "strand": strand,
                    "residue": fields[1].upper(),
                }
            )

    if len(proteins) != npro:
        raise RuntimeError(
            f"Topology declares npro={npro}, "
            f"but {len(proteins)} protein rows found"
        )

    return proteins


def base_records(proteins):
    records = []

    for i, p in enumerate(proteins):

        # Non-zero values on a few residues are enough to
        # exercise runtime loading while remaining benign.
        if i % 17 == 0:
            charge = 1.0
        elif i % 19 == 0:
            charge = -1.0
        else:
            charge = 0.0

        records.append(
            {
                "index": i,
                "strand": p["strand"],
                "residue": p["residue"],
                "charge": charge,
                "group": i,
            }
        )

    return records


def render_v1(
    records,
    *,
    version=1,
    index_space="protein_local",
    index_base="0",
    n_protein_beads=None,
    n_groups=None,
    include_index_space=True,
    include_index_base=True,
    include_n_protein=True,
    include_n_groups=True,
):
    if n_protein_beads is None:
        n_protein_beads = len(records)

    if n_groups is None:
        n_groups = len(
            set(
                r["group"]
                for r in records
            )
        )

    lines = [
        f"# DNANM_ES_CHARGES {version}",
    ]

    if include_index_space:
        lines.append(
            f"# index_space = {index_space}"
        )

    if include_index_base:
        lines.append(
            f"# index_base = {index_base}"
        )

    lines.extend(
        [
            "# charge_model = fail-closed-test",
            "# charge_source = fail-closed-test",
        ]
    )

    if include_n_protein:
        lines.append(
            f"# n_protein_beads = {n_protein_beads}"
        )

    if include_n_groups:
        lines.append(
            f"# n_groups = {n_groups}"
        )

    lines.append(
        "# columns = protein_index "
        "topology_strand residue charge_e group_id"
    )

    for r in records:
        lines.append(
            f"{r['index']} "
            f"{r['strand']} "
            f"{r['residue']} "
            f"{r['charge']} "
            f"{r['group']}"
        )

    return "\n".join(lines) + "\n"


def render_legacy(records):
    return (
        "\n".join(
            f"{r['index']} "
            f"{r['charge']} "
            f"{r['group']}"
            for r in records
        )
        + "\n"
    )


def combined_output(case, proc):
    log = ""

    logfile = case / "log.dat"

    if logfile.exists():
        log = logfile.read_text(
            errors="replace"
        )

    return (
        proc.stdout
        + "\n"
        + log
    )


def prepare_case(
    root,
    name,
    mode,
    *,
    charge_content=None,
    create_charge_file=True,
    es_enabled=True,
    extra_input="",
):
    case = root / name

    shutil.copytree(
        KDPG,
        case,
    )

    input_name = (
        "input_dnanm"
        if mode == "DNANM"
        else "input_dnact"
    )

    inp = case / input_name

    text = inp.read_text()

    text = replace_setting(
        text,
        "steps",
        "1",
    )

    text = replace_setting(
        text,
        "seed",
        "424242",
    )

    text = replace_setting(
        text,
        "print_energy_every",
        "1",
    )

    text = replace_setting(
        text,
        "print_conf_interval",
        "1",
    )

    if es_enabled:
        text += (
            "\n"
            "protein_dna_electrostatics = true\n"
            "protein_protein_electrostatics = true\n"
            "protein_charge_file = protein_charges.dat\n"
            "protein_charge_scale = 1.0\n"
            "protein_dna_dh_scale = 1.0\n"
            "protein_protein_dh_scale = 1.0\n"
            "dna_effective_charge_e = 0.815\n"
            "exclude_intragroup_electrostatics = true\n"
        )
    else:
        # Deliberately invalid path. With ES OFF it must
        # remain irrelevant to upstream behaviour.
        text += (
            "\n"
            "protein_dna_electrostatics = false\n"
            "protein_protein_electrostatics = false\n"
            "protein_charge_file = DOES_NOT_EXIST.dat\n"
        )

    text += extra_input

    inp.write_text(text)

    if create_charge_file:
        if charge_content is None:
            raise RuntimeError(
                "create_charge_file=True but no content supplied"
            )

        (case / "protein_charges.dat").write_text(
            charge_content
        )

    return case, input_name


def run_case(
    root,
    name,
    mode,
    *,
    charge_content=None,
    create_charge_file=True,
    es_enabled=True,
    extra_input="",
):
    case, input_name = prepare_case(
        root,
        name,
        mode,
        charge_content=charge_content,
        create_charge_file=create_charge_file,
        es_enabled=es_enabled,
        extra_input=extra_input,
    )

    proc = subprocess.run(
        [
            str(OXDNA),
            input_name,
        ],
        cwd=case,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=90,
    )

    return (
        proc.returncode,
        combined_output(
            case,
            proc,
        ),
    )


class FailClosedMatrix(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        if not OXDNA.is_file():
            raise RuntimeError(
                f"Invalid OXDNA_BIN: {OXDNA}"
            )

        cls.proteins = parse_proteins(
            KDPG / "generated.top"
        )

        cls.records = base_records(
            cls.proteins
        )

        cls.valid_v1 = render_v1(
            cls.records
        )

        cls.valid_legacy = render_legacy(
            cls.records
        )

    def assert_rejected(
        self,
        mode,
        label,
        *,
        charge_content=None,
        create_charge_file=True,
        extra_input="",
        needles=(),
    ):
        with tempfile.TemporaryDirectory(
            prefix=f"oxdna_fail_{mode.lower()}_"
        ) as td:

            rc, output = run_case(
                Path(td),
                label,
                mode,
                charge_content=charge_content,
                create_charge_file=create_charge_file,
                es_enabled=True,
                extra_input=extra_input,
            )

            self.assertNotEqual(
                rc,
                0,
                msg=(
                    f"{mode}/{label} unexpectedly accepted\n\n"
                    f"{output}"
                ),
            )

            if needles:
                self.assertTrue(
                    any(
                        needle in output
                        for needle in needles
                    ),
                    msg=(
                        f"{mode}/{label} rejected, but not "
                        f"with the expected diagnostic.\n"
                        f"Expected one of: {needles}\n\n"
                        f"{output}"
                    ),
                )

    def assert_accepted(
        self,
        mode,
        label,
        *,
        charge_content=None,
        create_charge_file=True,
        es_enabled=True,
        extra_input="",
        required_text=(),
    ):
        with tempfile.TemporaryDirectory(
            prefix=f"oxdna_accept_{mode.lower()}_"
        ) as td:

            rc, output = run_case(
                Path(td),
                label,
                mode,
                charge_content=charge_content,
                create_charge_file=create_charge_file,
                es_enabled=es_enabled,
                extra_input=extra_input,
            )

            self.assertEqual(
                rc,
                0,
                msg=(
                    f"{mode}/{label} unexpectedly rejected\n\n"
                    f"{output}"
                ),
            )

            for needle in required_text:
                self.assertIn(
                    needle,
                    output,
                )

    # -------------------------------------------------------------
    # Positive controls
    # -------------------------------------------------------------

    def test_01_valid_v1_both_modes(self):
        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_accepted(
                    mode,
                    "valid_v1",
                    charge_content=self.valid_v1,
                    required_text=(
                        "Loaded",
                        "protein charges",
                    ),
                )

    def test_02_legacy_format_both_modes(self):
        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_accepted(
                    mode,
                    "legacy",
                    charge_content=self.valid_legacy,
                    required_text=(),
                )

    def test_03_es_off_ignores_charge_file(self):
        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_accepted(
                    mode,
                    "es_off",
                    charge_content=None,
                    create_charge_file=False,
                    es_enabled=False,
                )

    # -------------------------------------------------------------
    # File presence / structure
    # -------------------------------------------------------------

    def test_04_missing_charge_file(self):
        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "missing_file",
                    charge_content=None,
                    create_charge_file=False,
                    needles=(
                        "Cannot open protein charge file",
                    ),
                )

    def test_05_empty_charge_file(self):
        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "empty",
                    charge_content="",
                    needles=(
                        "Charge file has 0 entries, expected 605",
                    ),
                )

    def test_06_unsupported_version(self):
        bad = render_v1(
            self.records,
            version=2,
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "version2",
                    charge_content=bad,
                    needles=(
                        "Invalid versioned charge header on line 1",
                    ),
                )

    # -------------------------------------------------------------
    # Metadata validation
    # -------------------------------------------------------------

    def test_07_wrong_index_space(self):
        bad = render_v1(
            self.records,
            index_space="global",
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "wrong_index_space",
                    charge_content=bad,
                    needles=(
                        "Charge indices must be zero-based protein_local",
                    ),
                )

    def test_08_wrong_index_base(self):
        bad = render_v1(
            self.records,
            index_base="1",
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "wrong_index_base",
                    charge_content=bad,
                    needles=(
                        "Charge indices must be zero-based protein_local",
                    ),
                )

    def test_09_missing_required_metadata(self):
        bad = render_v1(
            self.records,
            include_index_space=False,
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "missing_metadata",
                    charge_content=bad,
                    needles=(
                        "Charge file missing metadata 'index_space'",
                    ),
                )

    def test_10_wrong_n_protein_beads(self):
        bad = render_v1(
            self.records,
            n_protein_beads=(
                len(self.records) + 1
            ),
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "wrong_npro",
                    charge_content=bad,
                    needles=(
                        "Charge file topology protein count mismatch",
                    ),
                )

    def test_11_wrong_n_groups(self):
        bad = render_v1(
            self.records,
            n_groups=(
                len(self.records) + 1
            ),
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "wrong_ngroups",
                    charge_content=bad,
                    needles=(
                        "Charge metadata n_groups does not match file",
                    ),
                )

    # -------------------------------------------------------------
    # Record-level fail closed
    # -------------------------------------------------------------

    def test_12_duplicate_index(self):
        records = [
            dict(r)
            for r in self.records
        ]

        records[1] = dict(records[0])

        bad = render_v1(
            records,
            n_protein_beads=len(self.records),
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "duplicate",
                    charge_content=bad,
                    needles=(
                        "Duplicate or invalid protein charge local index",
                    ),
                )

    def test_13_missing_index(self):
        records = [
            dict(r)
            for r in self.records[:-1]
        ]

        bad = render_v1(
            records,
            n_protein_beads=len(self.records),
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "missing_index",
                    charge_content=bad,
                    needles=(
                        "contains",
                        "expected",
                    ),
                )

    def test_14_out_of_range_index(self):
        records = [
            dict(r)
            for r in self.records
        ]

        records[0]["index"] = len(
            self.records
        )

        bad = render_v1(
            records
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "out_of_range",
                    charge_content=bad,
                    needles=(
                        "Duplicate or invalid protein charge local index",
                    ),
                )

    def test_15_strand_mismatch(self):
        records = [
            dict(r)
            for r in self.records
        ]

        records[0]["strand"] = -999999

        bad = render_v1(
            records
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "strand_mismatch",
                    charge_content=bad,
                    needles=(
                        "Protein charge/topology strand or residue mismatch",
                    ),
                )

    def test_16_residue_mismatch(self):
        records = [
            dict(r)
            for r in self.records
        ]

        actual = records[0]["residue"]

        records[0]["residue"] = (
            "A"
            if actual != "A"
            else "K"
        )

        bad = render_v1(
            records
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "residue_mismatch",
                    charge_content=bad,
                    needles=(
                        "Protein charge/topology strand or residue mismatch",
                    ),
                )

    def test_17_negative_group_id(self):
        records = [
            dict(r)
            for r in self.records
        ]

        records[0]["group"] = -1

        bad = render_v1(
            records
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "negative_group",
                    charge_content=bad,
                    needles=(
                        "Negative charge group for protein",
                    ),
                )

    def test_18_extra_record_field(self):
        lines = self.valid_v1.splitlines()

        for i, line in enumerate(lines):

            if line and not line.startswith("#"):
                lines[i] = (
                    line
                    + " UNEXPECTED_FIELD"
                )
                break

        bad = "\n".join(lines) + "\n"

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "extra_field",
                    charge_content=bad,
                    needles=(
                        "Malformed version-1 charge record",
                    ),
                )

    def test_19_malformed_record(self):
        lines = self.valid_v1.splitlines()

        for i, line in enumerate(lines):

            if line and not line.startswith("#"):

                fields = line.split()

                # Remove group_id.
                lines[i] = " ".join(
                    fields[:-1]
                )

                break

        bad = "\n".join(lines) + "\n"

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "malformed_record",
                    charge_content=bad,
                    needles=(
                        "Malformed version-1 charge record",
                    ),
                )

    def test_20_nonfinite_charge(self):
        lines = self.valid_v1.splitlines()

        for i, line in enumerate(lines):

            if line and not line.startswith("#"):

                fields = line.split()
                fields[3] = "nan"

                lines[i] = " ".join(
                    fields
                )

                break

        bad = "\n".join(lines) + "\n"

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(mode=mode):
                self.assert_rejected(
                    mode,
                    "nonfinite_charge",
                    charge_content=bad,
                    needles=(
                        "Non-finite charge",
                        "Malformed version-1 charge record",
                        "out of range or non-finite",
                    ),
                )

    # -------------------------------------------------------------
    # Input-level numeric validation
    # -------------------------------------------------------------

    def numeric_rejection(
        self,
        key,
        value,
        diagnostic,
    ):
        extra = (
            f"\n{key} = {value}\n"
        )

        for mode in (
            "DNANM",
            "DNACT",
        ):
            with self.subTest(
                mode=mode,
                key=key,
                value=value,
            ):
                self.assert_rejected(
                    mode,
                    (
                        key
                        .replace("_", "-")
                    ),
                    charge_content=self.valid_v1,
                    extra_input=extra,
                    needles=diagnostic,
                )

    def test_21_negative_protein_charge_scale(self):
        self.numeric_rejection(
            "protein_charge_scale",
            "-1.0",
            (
                "protein_charge_scale must be non-negative",
            ),
        )

    def test_22_negative_pd_scale(self):
        self.numeric_rejection(
            "protein_dna_dh_scale",
            "-0.1",
            (
                "protein_dna_dh_scale must be non-negative",
            ),
        )

    def test_23_negative_pp_scale(self):
        self.numeric_rejection(
            "protein_protein_dh_scale",
            "-0.1",
            (
                "protein_protein_dh_scale must be non-negative",
            ),
        )

    def test_24_nonpositive_dna_effective_charge(self):
        self.numeric_rejection(
            "dna_effective_charge_e",
            "0.0",
            (
                "dna_effective_charge_e must be positive",
            ),
        )

    def test_25_invalid_numeric_string(self):
        self.numeric_rejection(
            "protein_charge_scale",
            "not-a-number",
            (
                "Numeric input key `protein_charge_scale'",
                "invalid value",
            ),
        )

    def test_26_nonfinite_numeric_input(self):
        # All spellings below are accepted by strtod-family
        # parsers as non-finite values and must fail closed.
        for value in (
            "nan",
            "+nan",
            "-nan",
            "nan(payload)",
            "NaN",
            "inf",
            "+inf",
            "-inf",
            "infinity",
            "-infinity",
        ):
            with self.subTest(value=value):
                self.numeric_rejection(
                    "protein_charge_scale",
                    value,
                    (
                        "out of range or non-finite",
                    ),
                )

    def test_27_overflow_numeric_input(self):
        self.numeric_rejection(
            "protein_charge_scale",
            "1e9999",
            (
                "out of range or non-finite",
            ),
        )


if __name__ == "__main__":
    unittest.main(
        verbosity=2
    )
