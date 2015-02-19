"""Define constants and register exit functions"""
import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.7 or above is required for NAGlib.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed
del sys

# create a temporary directory in which to work
import os
TEMPDIR = '/tmp/naglib/'
if not os.path.exists(TEMPDIR):
    os.makedirs(TEMPDIR)
del os


def __naglib_debug():
    # helper function so we don't import os globally
    import os
    debug_str = os.getenv('NAGLIB_DEBUG', 'False')
    if debug_str in ('True', 'False'):
        return eval(debug_str)
    else:
        raise RuntimeError('unrecognized value for NAGLIB_DEBUG: {0}'.format(debug_str))
NAGLIB_DEBUG = __naglib_debug()

import atexit
@atexit.register
def cleanup():
    if NAGLIB_DEBUG:
        pass
    # bug here
    #else:
        #import os
        #import shutil
        #if os.path.exists(TEMPDIR):
            #shutil.rmtree(TEMPDIR)
del atexit