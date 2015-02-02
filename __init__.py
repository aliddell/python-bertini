"""NAGlib is a Python library for numerical algebraic geometry. NAGlib aims to
be an interface to Bertini 2 and any other software packages using techniques
from NAG. It depends on NumPy, SciPy, SymPy, and MPMath."""

from __future__ import absolute_import, print_function

from naglib.release import __version__

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.7 or above is required for NAGlib.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys

def __naglib_debug():
    # helper function so we don't import os globally
    import os
    debug_str = os.getenv('NAGLIB_DEBUG', 'False')
    if debug_str in ('True', 'False'):
        return eval(debug_str)
    else:
        raise RuntimeError('unrecognized value for NAGLIB_DEBUG: {0}'.format(debug_str))
NAGLIB_DEBUG = __naglib_debug()

from .bertini import *
#from .core import *