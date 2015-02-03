from __future__ import print_function

import re

from fractions import Fraction as fraction
from os import chdir
from tempfile import mkdtemp

from sympy import Matrix

from naglib.core.datatypes import Component
from naglib.bertini.sysutils import call_bertini
from fileutils import striplines, write_system

def compute_NID(system):
    dirname = mkdtemp()
    input_file = dirname + 'input'

    variables = set()
    for f in system:
        variables = variables.union(f.free_symbols)

    write_system(variables=variables, params=None, functions=system,
                 constants=None, tracktype=1, filename=input_file)
    chdir(dirname)
    call_bertini(input_file)

    components = get_components(dirname)

    return components

def get_components(dirname='.'):
    # this is a quite crude way to get the NID
    # this will change in the near future
    main_data = dirname + '/main_data'
    fh = open(main_data, 'r')
    lines = striplines(fh.readlines())
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
                components.append(Component(dim,1,None))
                dimcomponents.add(c)
            else:
                components[-1].degree += 1

    return components

def jacobian(F, variables):
    return Matrix([[f.diff(v) for v in variables] for f in F])

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
                qi = (fraction(pi[0]).limit_denominator(d_limit),fraction(pi[1]))
            else:
                qi = (fraction(pi[0]),fraction(pi[1]))
            q.append(qi)
        new_points.append(tuple(q))

    return new_points