from __future__ import print_function

import os
import re

from fractions import Fraction as fraction

from sympy import Matrix

from fileutils import striplines

def get_components(filename):
    dirname = os.path.dirname(filename)
    if not dirname:
        dirname = '.'
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
                components.append(dict([('dim',dim),('id',c), ('degree',1),
                                        ('directory', dirname)]))
                dimcomponents.add(c)
            else:
                components[-1]['degree'] += 1

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