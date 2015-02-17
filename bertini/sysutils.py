from __future__ import print_function

from os import chdir
from os.path import dirname
from subprocess import check_output, CalledProcessError

from naglib.exceptions import BertiniError, NoBertiniException

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

def __proc_err_output(output):
    lines = output.split('\n')
    dex = -1
    for l in lines:
        if l.startswith('ERROR'):
            dex = lines.index(l)
    if dex == -1:
        return output
    else:
        l = lines[dex]
        # strip 'ERROR: '
        lines[dex] = l[l.index(' ')+1:]
        # strip 'Bertini will now exit due to this error'
        return '\n'.join(lines[dex:-1])

BERTINI = __has_bertini()

def call_bertini(input_file, start_file='', cmd=BERTINI, stdin=None):
    """
    Call Bertini
    
    Keyword arguments:
    input_file -- string, the path of a Bertini input file
    start_file -- optional string, the path of a Bertini start file
    cmd        -- optional string, the path of the Bertini executable
    redirect   -- optional string, path to stdin redirect
    """
    if not cmd:
        raise(NoBertiniError)
    if not start_file:
        arg = [cmd, input_file]
    elif start_file:
        arg = [cmd, input_file, start_file]
        
    wd = dirname(input_file)
    if wd:
        chdir(wd)
    if stdin:
        stdin = open(stdin, 'r')
    try:
        output = check_output(arg, stdin=stdin)
    except CalledProcessError as e:
        raise(BertiniError(__proc_err_output(e.output)))

    if stdin:
        stdin.close()
    return output