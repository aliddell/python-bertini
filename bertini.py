from __future__ import print_function

from fractions import Fraction
from os.path import isfile
import os
import re
import sys

import sympy

# INTERNAL METHODS
def _get_components(filename):
    dirname = os.path.dirname(filename)
    if not dirname:
        dirname = '.'
    main_data = dirname + '/main_data'
    fh = open(main_data, 'r')
    lines = _striplines(fh.readlines())
    fh.close()

    components = []
    dim = -1
    dimcomponents = set()
    for l in lines:
        if l.startswith('----------DIMENSION'):
            dim = int(re.sub(r'(DIMENSION|\-|\s+)', '', l))
            dimcomponents = set()
        elif l.startswith('Component number'):
            # capture number at end of line
            c = int(re.sub(r'[^(\d+$)]', '', l))
            if c not in dimcomponents:
                components.append(dict([('dim',dim),('id',c), ('degree',1),
                                        ('directory', dirname)]))
                dimcomponents.add(c)
            else:
                components[-1]['degree'] += 1

    return components

def _jacobian(F, variables):
    return sympy.Matrix([[f.diff(v) for v in variables] for f in F])

def _read_input(filename):
    fh = open(filename, 'r')
    lines = _striplines(fh.readlines())
    fh.close()

    var_candidates = []
    fun_candidates = []
    for line in lines:
        if 'variable_group' in line:
            var_candidates.append(line)
        if 'function' in line:
            fun_candidates.append(line)

    if len(var_candidates) == 0 or len(fun_candidates) == 0:
        exit('malformed system file, pal')

    if len(var_candidates) == 1:
        variable_group = var_candidates[0]
    else:
        variable_group = None
        for v in var_candidates:
            if v.startswith('variable_group'):
                variable_group = v
                break
        if not variable_group:
            exit("where's your variable_group?")
    if len(fun_candidates) == 1:
        function_group = fun_candidates[0]
    else:
        function_group = None
        for f in fun_candidates:
            if f.startswith('function'):
                function_group = f
                break
        if not function_group:
            exit("where's your function?")

    variable_group = re.sub(string=variable_group, pattern=r'variable_group\s+', repl='')
    variable_group = re.sub(string=variable_group, pattern=r'\s*;\s*$', repl='')
    variable_group = re.split(string=variable_group, pattern=r',\s*')

    function_group = re.sub(string=function_group, pattern=r'function\s+', repl='')
    function_group = re.sub(string=function_group, pattern=r'\s*;\s*$', repl='')
    function_group = re.split(string=function_group, pattern=r',\s*')

    # break out the functions
    functions = []
    for f in function_group:
        for line in lines:
            if line.startswith(f) and '=' in line:
                functions.append(re.sub(string=line, pattern=r'\s*;\s*', repl=''))
    if len(functions) != len(function_group):
        exit("you either haven't defined all your functions or you've defined too many")

    return (functions, variable_group)

def _striplines(lines, nonempty=True):
    if nonempty:
        return [l.strip() for l in lines if l != '\n']
    else:
        return [l.strip() for l in lines]

def _parselines(lines, **kwargs):
    """Returns a list of n-tuples of complex numbers, realized as ordered pairs"""
    lines = _striplines(lines)
    points = []

    d_limit = 10**8

    # parameters
    if 'tol' in kwargs:
        withtol = True
        tol = kwargs['tol']
    else:
        withtol = False
    if 'rational' in kwargs and kwargs['rational']:
        rational = True
    else:
        rational = False
    if 'limit' in kwargs and kwargs['limit']:
        limit = True
    else:
        limit = False

    numpoints = int(lines[0])
    lines = lines[1:]
    length = len(lines)

    numvar = length//numpoints

    for i in range(0,length,numvar):
        point = lines[i:i+numvar]
        point = [re.split(r'\s+', p) for p in point]
        if withtol:
            newpoint = []
            for p in point:
                if rational and limit:
                    p0 = Fraction(p[0]).limit_denominator(d_limit)
                    p1 = Fraction(p[1]).limit_denominator(d_limit)
                elif rational:
                    p0 = Fraction(p[0])
                    p1 = Fraction(p[1])
                else:
                    p0 = float(p[0])
                    p1 = float(p[1])

                if abs(p0) < tol:
                    p[0] = 0
                else:
                    p[0] = p0
                if abs(p1) < tol:
                    p[1] = 0
                else:
                    p[1] = p1
                newpoint.append(tuple(p))
            points.append(tuple(newpoint))
        else:
            if rational and limit:
                point = [(Fraction(p[0]).limit_denominator(d_limit), Fraction(p[1]).limit_denominator(d_limit)) for p in point]
            elif rational:
                point = [(Fraction(p[0]), Fraction(p[1])) for p in point]
            else:
                point = [(float(p[0]), float(p[1])) for p in point]

            points.append(tuple(point))

    return list(set(points))

def fprint(points, **kwargs):
    """Prints a set of points in Bertini output fashion, optionally to a file"""
    if 'filename' in kwargs:
        usestdout = False
        filename = kwargs['filename']
        if isfile(filename) and 'overwrite' not in kwargs:
            res = input('{0} already exists! Overwrite? '.format(filename)).lower()
            if res[0] == 'y':
                fh = open(filename, 'w')
            else:
                print('Bailing out', file=sys.stderr)
                res = None
                return res
        elif isfile(filename):
            if kwargs['overwrite']:
                fh = open(filename, 'w')
            else:
                print('Bailing out', file=sys.stderr)
                res = None
                return res
        else:
            fh = open(filename, 'w')
    else:
        usestdout = True
        fh = sys.stdout

    numpoints = len(points)
    print('{0}'.format(numpoints), file=fh)
    for p in points:
        for q in p:
            print('{0} {1}'.format(q[0], q[1]), file=fh)
        print('', file=fh)

    if not usestdout:
        fh.close()
    return

def readfile(filename, **kwargs):
    """Reads in a file and returns a set of n-tuples of complex numbers"""
    if not isfile(filename):
        return None

    # parameters
    if 'tol' in kwargs:
        withtol = True
        tol = kwargs['tol']
    else:
        withtol = False
    if 'rational' in kwargs and kwargs['rational']:
        rational = True
    else:
        rational = False
    if 'limit' in kwargs and kwargs['limit']:
        limit = True
    else:
        limit = False

    fh = open(filename, 'r')
    lines = fh.readlines()
    fh.close()

    points = _parselines(lines,**kwargs)
    return points

def real(points, **kwargs):
    """Returns the subset of real points from a given complex set"""
    # parameters
    if 'tol' in kwargs:
        tol = kwargs['tol']
    else:
        tol = 1e-10


    real_points = []

    for p in points:
        isreal = True
        for pi in p:
            compart = pi[1]
            if abs(compart) > tol:
                isreal = False
                break

        if isreal:
            real_points.append(p)

    return real_points

def tofloat(points):
    new_points = []

    for p in points:
        q = []
        for pi in p:
            q.append((float(pi[0]),float(pi[1])))
        new_points.append(tuple(q))

    return new_points

def torational(points,**kwargs):
    if 'limit' in kwargs and kwargs['limit']:
        limit = True
        d_limit = 10**8
    else:
        limit = False

    new_points = []

    for p in points:
        q = []
        for pi in p:
            if limit:
                qi = (Fraction(pi[0]).limit_denominator(d_limit),Fraction(pi[1]))
            else:
                qi = (Fraction(pi[0]),Fraction(pi[1]))
            q.append(qi)
        new_points.append(tuple(q))

    return new_points

######################## EXTERNALLY FACING STUFF ##############################

def get_components(input_file):
    """Return a list of (dimension, component_id) ordered pairs"""
    return _get_components(input_file)

def jacobian(F, variables):
    """Return ths Jacobian of F with respect to variables"""
    return _jacobian(F, variables)

def read_input(input_file):
    """Read in an input file and return lists of functions and variables"""
    return _read_input(input_file)

def striplines(lines, nonempty=True):
    """Return a list of (potentially nonempty) strings,
    stripped of newlines"""
    return _striplines(lines, nonempty)