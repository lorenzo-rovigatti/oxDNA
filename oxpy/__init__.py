import sys, os
# we need to load the library so that its symbols are visible outside or plugins won't work. see here: https://stackoverflow.com/a/60841073/5140209
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

from .core import *

from . import utils

from .__version__ import __title__, __description__, __url__, __version__
from .__version__ import __author__, __author_email__, __license__, __copyright__
