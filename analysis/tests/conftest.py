import logging
from pathlib import Path
from shutil import copy
from typing import Any
from typing import Dict

import pytest


MINI_TRAJ = "minitraj.dat"
ALIGN_TRAJ = "aligntraj.dat"
MEAN_CONF = "mean.dat"
INDEX = "index.txt"
PAIRS = "pairs.txt"

RNA_TILE_TOP = "rna_tile.top"
INPUT_RNA = "input_rna"
SEQUENCE_DEPS_RNA = "rna_sequence_dependent_parameters.txt"


logger = logging.getLogger(__name__)


def get_test_resource(file_name: str) -> Path:
    base_dir = Path.cwd()
    return base_dir / "tests" / "resources" / file_name


@pytest.fixture(scope="session")
def align_input(tmp_path_factory: pytest.FixtureRequest) -> Dict[str, Any]:

    dest_dir = tmp_path_factory.mktemp("test_cli")
    data = {}

    mini_traj = get_test_resource(f"{MINI_TRAJ}")
    data["traj"] = copy(mini_traj, dest_dir)
    data["outfile"] = dest_dir / ALIGN_TRAJ
    data["indexes"] = []
    data["ref_conf"] = None
    data["center"] = True
    data["ncpus"] = 1

    test_result = get_test_resource(f"{ALIGN_TRAJ}")
    data["test"] = test_result

    return data

@pytest.fixture(scope="session")
def mean_input(tmp_path_factory: pytest.FixtureRequest) -> Dict[str, Any]:
    from oxDNA_analysis_tools.UTILS.RyeReader import describe
    dest_dir = tmp_path_factory.mktemp("test_cli")
    data = {}

    mini_traj = get_test_resource(f"{MINI_TRAJ}")
    data["top_info"], data["traj_info"] = describe(None, mini_traj.as_posix())
    data["indexes"] = []
    data["ref_conf"] = None
    data["ncpus"] = 1

    test_result = get_test_resource(f"{MEAN_CONF}")
    data["test"] = test_result

    return data

@pytest.fixture(scope="session")
def mean_cli_input(tmp_path_factory: pytest.FixtureRequest) -> Dict[str, Any]:

    data = {}

    mini_traj = get_test_resource(f"{MINI_TRAJ}")
    data["traj_file"] = mini_traj

    test_result = get_test_resource(f"{MEAN_CONF}")
    data["test"] = test_result

    return data