from __future__ import print_function

from subprocess import check_output, CalledProcessError

class NoBertiniError(Exception):

    def __init__(self):
        self.message = "You don't seem to have Bertini installed anywhere " \
                       "I can find it."
    def __str__(self):
        return self.message

def __os():
    from sys import platform

    if platform.startswith('win'):
        return 'WINDOWS'
    elif platform.startswith('cygwin'):
        return 'CYGWIN'
    elif platform.startswith('linux'):
        return 'LINUX'
    elif platform.startswith('darwin'):
        return 'OSX'

def __has_bertini():
    platform = __os()
    if platform == 'WINDOWS':
        cmd = 'where.exe'
    else:
        cmd = 'which'

    try:
        bertinipath = check_output([cmd, 'bertini'])
    except CalledProcessError:
        bertinipath = ''

    return bertinipath.strip()

HAS_BERTINI = __has_bertini()

def bertini(input_file, start_file='', cmd=HAS_BERTINI, suppress=True):
    if not cmd:
        raise(NoBertiniError)
    if not start_file:
        arg = [cmd, input_file]
    else:
        arg = [cmd, input_file, start_file]
    output = check_output(arg)
    if not suppress:
        print(output)