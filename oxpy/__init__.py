import sys, os
# we need to load the library so that its symbols are visible outside or plugins won't work. See here: https://stackoverflow.com/a/60841073/5140209
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

from .core import *

from . import utils

__title__ = 'oxpy'
__description__ = 'Python bindings for oxDNA'
__url__ = 'https://github.com/lorenzo-rovigatti/oxDNA'
__author__ = 'Lorenzo Rovigatti, Flavio Romano, Petr Sulc and others'
__author_email__ = 'lorenzo.rovigatti@uniroma1.it'
__license__ = 'GNU GPL 3.0'
__copyright__ = 'Copyright 2025 Lorenzo Rovigatti, Flavio Romano, Petr Sulc and others'

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
