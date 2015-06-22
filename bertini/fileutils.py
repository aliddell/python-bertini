from __future__ import print_function

from os.path import isfile
from sys import stdout

from sympy import I, Float, sympify

from naglib.startup import TOL
from naglib.core import AffinePoint, ProjectivePoint
from naglib.core.misc import striplines

def parselines(lines, tol=TOL, projective=False, as_set=False):
    """    
    Keyword arguments:
    lines -- iterable of strings, first entry the number of points;
             the rest, "%s %s" % real, imag
    tol   -- optional float, smallest allowable nonzero value
    """
    from re import split as resplit
    from naglib.core.misc import dps
    
    lines = striplines(lines)
    points = []

    numpoints = int(lines[0])
    if numpoints == 0:
        return points
    lines = lines[1:]
    length = len(lines)
    numvar = length//numpoints

    for i in range(0, length, numvar):
        point = lines[i:i+numvar]
        point = [resplit(r'\s+', p) for p in point]
        newpoint = []
        for p in point:
            real,imag = p
            if imag[-1] == ';':
                imag = imag[:-1]
                
            real = Float(real, dps(real))
            imag = Float(imag, dps(imag))

            if abs(real) < tol:
                real = 0
            if abs(imag) < tol:
                imag = 0
                
            newpoint.append(real + I*imag)
        if projective:
            points.append(ProjectivePoint(newpoint))
        else:
            points.append(AffinePoint(newpoint))

    if as_set:
        points = list(set(points))

    return points

def read_points(filename, tol=TOL, projective=False, as_set=False):
    """
    Reads in a file and return a set of Float numbers
    """
    if not isfile(filename):
        msg = "{0} does not exist".format(filename)
        raise IOError(msg)

    fh = open(filename, 'r')
    lines = striplines(fh.readlines())
    fh.close()

    points = parselines(lines, tol=tol, projective=projective, as_set=as_set)
    return points

# write utils

def fprint(points, filename=''):
    """
    Print a set of points in Bertini output fashion, optionally to a file
    
    Keyword arguments:
    points   -- numeric iterable, the points to print
    filename -- optional string, path to filename
    """
    if filename:
        fh = open(filename, 'w')
    else:
        fh = stdout

    numpoints = len(points)
    print('{0}\n'.format(numpoints), file=fh)
    for p in points:
        for coordinate in p.coordinates:
            coordinate = sympify(coordinate)
            real, imag = coordinate.as_real_imag()
            print('{0} {1}'.format(real, imag), file=fh)
        print('', file=fh)

    if filename:
        fh.close()
    return
    