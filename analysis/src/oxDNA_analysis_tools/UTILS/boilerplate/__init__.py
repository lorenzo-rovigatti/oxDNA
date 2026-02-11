"""
oxDNA_analysis_tools.UTILS.boilerplate -- Simulation setup and management.

This package provides a high-level interface for setting up, running, and
analyzing oxDNA simulations. It includes pre-defined protocol parameter
sets for common workflow stages (MC relaxation, MD relaxation, production).

Backward compatibility:
    All public names from the original boilerplate module are re-exported
    here, so existing imports continue to work unchanged::

        from oxDNA_analysis_tools.UTILS.boilerplate import Simulation
        from oxDNA_analysis_tools.UTILS.boilerplate import setup_simulation
        from oxDNA_analysis_tools.UTILS.boilerplate import get_default_input
"""

# Utilities
from .utils import (
    SimulationRunningError,
    PathContext,
    path_decorator,
    _prun,
)

# Protocols and parameter dictionaries
from .protocols import (
    default_input_file,
    get_default_input,
    get_relax_MC_input,
    get_relax_MD_input,
    get_production_MD_input,
)

# Setup functions
from .setup import (
    dump_json,
    setup_simulation,
)

# The Simulation class
from .simulation import Simulation

__all__ = [
    # Original public API (backward compatible)
    "default_input_file",
    "SimulationRunningError",
    "PathContext",
    "path_decorator",
    "get_default_input",
    "dump_json",
    "setup_simulation",
    "_prun",
    "Simulation",
    # New protocol getters
    "get_relax_MC_input",
    "get_relax_MD_input",
    "get_production_MD_input",
]
