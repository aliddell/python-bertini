"""NAGlib is a Python library for numerical algebraic geometry. NAGlib aims to
be an interface to Bertini 2 and any other software packages using techniques
from NAG. It depends on NumPy, SciPy, SymPy, and MPMath."""

from __future__ import absolute_import, print_function

from .startup import *
from .release import __version__, __authors__
from .exceptions import *
from .bertini import *
from .core import *