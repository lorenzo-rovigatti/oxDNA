import sys
try:
    from numpy import *
    print >> sys.stderr, "## using default numpy"
except:
    from mynumpy_core import *
    print >> sys.stderr, "## Warning: using crappy, stripped down, slow version of numpy."

