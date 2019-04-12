"""NAGlib is a Python library for numerical algebraic geometry. NAGlib aims to
be an interface to Bertini and any other software packages using techniques
from numerical algebraic geometry. It depends on NumPy, SciPy, and SymPy."""

from .exceptions import *
from .release import __version__, __authors__ # this is actually everything

from .bertini import *
from .core import *
