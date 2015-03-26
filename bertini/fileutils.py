from __future__ import print_function

from fractions import Fraction as fraction
from os.path import isfile
from sys import stderr, stdout

from sympy import I, Integer, Float, Rational, Matrix as spmatrix, sympify

from naglib.exceptions import UnclassifiedException
from naglib.startup import TOL, DPS
from naglib.core import AffinePoint, ProjectivePoint, WitnessPoint
from naglib.core.misc import striplines

def dps(s):
    return max(DPS, len(s))

# read utils

def parse_witness_data(filename):
    """
    Parse witness_data file into usable data
    
    Keyword arguments:
    filename -- string, path to witness_data file
    """
    if not isfile(filename):
        msg = "{0} does not exist".format(filename)
        raise IOError(msg)
    fh = open(filename, 'r')
    lines = striplines(fh.readlines())
    fh.close()

    num_vars, nonempty_codims = int(lines[0]), int(lines[1])
    # previous line includes additional homogenizing variable(s),
    # as appropriate
    lines = lines[2:]

    # following block represents a single codim; repeated for each
    codims = []
    for i in range(nonempty_codims):
        codim = int(lines[0])
        codims.append({'codim':codim})
        num_points = int(lines[1])
        lines = lines[2:]
        # following block represents a single point; repeated for each
        pts = []
        for j in range(num_points):
            prec = int(lines[0])

            lines = lines[1:]
            pt = []
            for k in range(num_vars):
                real,imag = lines[k].split(' ')
                pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))
            pt = spmatrix(pt)
            lines = lines[num_vars:]
            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines[0])
            lines = lines[1:]
            approx_pt = []
            for k in range(num_vars):
                real,imag = lines[k].split(' ')
                approx_pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))

            lines = lines[num_vars:]
            condition_number = float(lines[0])
            corank = float(lines[1]) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines[2])
            largest_zero_singular = float(lines[3])
            pt_type = int(lines[4])
            multiplicity = int(lines[5])
            component_number = int(lines[6])
            if component_number == -1:
                msg = "components in {0} have unclassified points".format(filename)
                raise UnclassifiedException(msg)
            deflations = int(lines[7])
            lines = lines[8:]
            pts.append({'coordinates':pt,
                        'corank':corank,
                        'condition number':condition_number,
                        'sigma_1':smallest_nonzero_singular,
                        'sigma_0':largest_zero_singular,
                        'type':pt_type,
                        'multiplicity':multiplicity,
                        'component number':component_number,
                        'deflations':deflations})
        codims[-1]['points'] = pts

    # -1 designates the end of witness points
    lines = lines[1:]

    INT = 0
    DOUBLE = 1
    RATIONAL = 2

    # remaining data is related to slices, randomization,
    # homogenization, and patches
    num_format = int(lines[0])
    # previous line describes format for remainder of data
    lines = lines[1:]

    # the following block is repeated for each nonempty codim.
    # first, matrix A used for randomization
    # second, matrix W
    for i in range(nonempty_codims):
        num_rows, num_cols = lines[0].split(' ')
        num_rows = int(num_rows)
        num_cols = int(num_cols)
        AW_size = num_rows*num_cols

        lines = lines[1:]

        if AW_size == 0:
            A = None
            W = None
        else:
            A = lines[:AW_size]
            lines = lines[AW_size:]
            W = lines[:AW_size]
            lines = lines[AW_size:]

            A = [a.split(' ') for a in A] # A is complex-valued
            if num_format == INT:
                A = [Integer(a[0]) + I*Integer(a[1]) for a in A]
            elif num_format == DOUBLE:
                for j in range(len(A)):
                    a = A[j]
                    real,imag = a.split(' ')
                    A[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                A = [Rational(a[0]) + I*Rational(a[1]) for a in A]
            A = [A[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            A = spmatrix(A)

            W = [int(w) for w in W] # W is integer-valued
            W = [W[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            W = spmatrix(W)

        # third, a vector H used for homogenization
        # random if projective input
        H_size = int(lines[0])
        lines = lines[1:]
        H = lines[:H_size]
        H = [h.split(' ') for h in H] # H is complex-valued
        if num_format == INT:
            H = [Integer(h[0]) + I*Integer(h[1]) for h in H]
        elif num_format == DOUBLE:
            for j in range(len(H)):
                h = H[j]
                real,imag = h.split(' ')
                H[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            H = [Rational(h[0]) + I*Rational(h[1]) for h in H]

        H = spmatrix(H)
        lines = lines[H_size:]

        # fourth, a number homVarConst
        # 0 for affine, random for projective
        hvc = lines[0].split(' ')
        if num_format == INT:
            hvc = Integer(hvc[0]) + I*Integer(hvc[1])
        elif num_format == DOUBLE:
            real,imag = hvc
            hvc = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            hvc = Rational(hvc[0]) + I*Rational(hvc[1])

        lines = lines[1:]

        # fifth, matrix B for linear slice coefficients
        num_rows, num_cols = lines[0].split(' ')
        num_rows, num_cols = int(num_rows), int(num_cols)
        B_size = num_rows*num_cols
        lines = lines[1:]

        if B_size == 0:
            B = None
        else:
            B = lines[:B_size]
            lines = lines[B_size:]

            B = [b.split(' ') for b in B] # B is complex-valued
            if num_format == INT:
                B = [Integer(b[0]) + I*Integer(b[1]) for b in B]
            elif num_format == DOUBLE:
                for j in range(len(B)):
                    real,imag = B[j]
                    B[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                B = [Rational(b[0]) + I*Rational(b[1]) for b in B]
            B = [B[j:j+num_cols] for j in range(0,B_size,num_cols)]
            B = spmatrix(B)

        # sixth and finally, vector p for patch coefficients
        p_size = int(lines[0])
        lines = lines[1:]

        p = lines[:p_size]
        p = [q.split(' ') for q in p]
        if num_format == INT:
            p = [Integer(q[0]) + I*Integer(q[1]) for q in p]
        elif num_format == DOUBLE:
            for j in range(len(p)):
                real,imag = p[j]
                p[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            p = [Rational(q[0]) + I*Rational(q[1]) for q in p]

        p = spmatrix(p)
        codims[i]['A'] = A
        codims[i]['W'] = W
        codims[i]['H'] = H
        codims[i]['homVarConst'] = hvc
        codims[i]['slice'] = B
        codims[i]['p'] = p

    return codims

def parselines(lines, tol=TOL, projective=False, as_set=False):
    """    
    Keyword arguments:
    lines -- iterable of strings, first entry the number of points;
                the rest, "%s %s" % real, imag
    tol   -- optional numeric, smallest allowable nonzero value
    """
    from re import split as resplit
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
    Reads in a file and return a set of mpc numbers
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

def write_input(system, config):
    """
    Writes a system to a Bertini-format file
    
    Keyword arguments:
    system -- PolynomialSystem, the system we want to solve
    config -- dict, the configurations for Bertini
    """
    # Python exponentiation: x**y; bertini exponentiation: x^y
    import re
    
    polynomials  = [p.factor() for p in system.polynomials]
    variables    = system.variables
    parameters   = system.parameters
    homvar       = system.homvar
    num_polys    = system.shape[0]
    
    options = config.keys()
    
    if 'filename' not in options:
        raise(IOError('no file to write to'))
    else:
        filename = config.pop('filename')
        options = config.keys()
    
    str_poly = [str(p) for p in polynomials]
    str_poly = [re.sub(string=p, pattern=r'\*\*', repl='^') for p in str_poly]
    str_vars = [str(v) for v in variables]
    str_pars = [str(p) for p in parameters]
    
    poly_names = ['f{0}'.format(i+1) for i in range(num_polys)]
    polys_named = zip(poly_names, str_poly)
    
    poly_list = ','.join([f for f in poly_names])
    vars_list = ','.join([v for v in str_vars])
    pars_list = ','.join([p for p in str_pars])

    fh = open(filename, 'w')
    
    # write the CONFIG section
    print('CONFIG', file=fh)
    for option in options:
        print('{0}:{1};'.format(option, config[option]), file=fh)
    print('END', file=fh)
    
    # write the INPUT section
    print('INPUT', file=fh)
    if parameters:
        print('parameter {0};'.format(pars_list), file=fh)
    if homvar:
        print('hom_variable_group {0};'.format(vars_list), file=fh)
    else:
        print('variable_group {0};'.format(vars_list), file=fh)
    print('function {0};'.format(poly_list), file=fh)
    
    for p in polys_named:
        # p is a key-value pair, e.g., ('f1', 'x^2 - 1')
        print('{0} = {1};'.format(p[0], p[1]), file=fh)
    print('END', file=fh)
    
    # finish up
    fh.close()

def write_witness_data(witness_data, dirname, components=[], filename='witness_data'):
    """
    """
    fh = open(dirname + '/' + filename, 'w')
    
    if components:
        compids = [(c.codim, c.component_id) for c in components]
        nonempty_codims = len(set(c[0] for c in compids))
    else:
        compids = []
        for codim in witness_data:
            c = codim['codim']
            points = codim['points']
            cids = list(set([p['component number'] for p in points]))
            cids.sort()
            for id in cids:
                compids.append((c, id))
        nonempty_codims = len(witness_data)
    num_vars = len(components[0]['points']['coordinates'])
    fh.write('{0}\n'.format(num_vars))
    fh.write('{0}\n'.format(nonempty_codims))
    
    unique_codims = [c[0] for c in compids]
    for i in range(nonempty_codims):
        codim = unique_codims[i]
        fh.write('{0}\n'.format(codim))
        codim_components = [c[1] for c in compids if c[0] == codim]
        counter = 0
        while counter < len(witness_data):
            if witness_data[counter]['codim'] == codim:
                break
            counter += 1
        wd_codim = witness_data[counter]
        codim_points = [p for p in wd_codim['points'] if p['component number'] in codim_components]
        for p in codim_points:
            coordinates = p['coordinates']
            for c in coordinates:
                real,imag = c.as_real_imag()
                fh.write('{0} {1}\n'.format(real, imag))
                
    fh.close()
    
    # following block represents a single codim; repeated for each
    #codims = []
    #for i in range(nonempty_codims):
        #codim = int(lines[0])
        #codims.append({'codim':codim})
        #num_points = int(lines[1])
        #lines = lines[2:]
        ## following block represents a single point; repeated for each
        pts = []
        for j in range(num_points):
            prec = int(lines[0])

            lines = lines[1:]
            pt = []
            for k in range(num_vars):
                real,imag = lines[k].split(' ')
                pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))
            pt = spmatrix(pt)
            lines = lines[num_vars:]
            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines[0])
            lines = lines[1:]
            approx_pt = []
            for k in range(num_vars):
                real,imag = lines[k].split(' ')
                approx_pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))

            lines = lines[num_vars:]
            condition_number = float(lines[0])
            corank = float(lines[1]) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines[2])
            largest_zero_singular = float(lines[3])
            pt_type = int(lines[4])
            multiplicity = int(lines[5])
            component_number = int(lines[6])
            if component_number == -1:
                msg = "components in {0} have unclassified points".format(filename)
                raise UnclassifiedException(msg)
            deflations = int(lines[7])
            lines = lines[8:]
            pts.append({'coordinates':pt,
                        'corank':corank,
                        'condition number':condition_number,
                        'sigma_1':smallest_nonzero_singular,
                        'sigma_0':largest_zero_singular,
                        'type':pt_type,
                        'multiplicity':multiplicity,
                        'component number':component_number,
                        'deflations':deflations,
                        'last approximation':approx_pt})
        codims[-1]['points'] = pts

    # -1 designates the end of witness points
    lines = lines[1:]

    INT = 0
    DOUBLE = 1
    RATIONAL = 2

    # remaining data is related to slices, randomization,
    # homogenization, and patches
    num_format = int(lines[0])
    # previous line describes format for remainder of data
    lines = lines[1:]

    # the following block is repeated for each nonempty codim.
    # first, matrix A used for randomization
    # second, matrix W
    for i in range(nonempty_codims):
        num_rows, num_cols = lines[0].split(' ')
        num_rows = int(num_rows)
        num_cols = int(num_cols)
        AW_size = num_rows*num_cols

        lines = lines[1:]

        if AW_size == 0:
            A = None
            W = None
        else:
            A = lines[:AW_size]
            lines = lines[AW_size:]
            W = lines[:AW_size]
            lines = lines[AW_size:]

            A = [a.split(' ') for a in A] # A is complex-valued
            if num_format == INT:
                A = [Integer(a[0]) + I*Integer(a[1]) for a in A]
            elif num_format == DOUBLE:
                for j in range(len(A)):
                    a = A[j]
                    real,imag = a.split(' ')
                    A[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                A = [Rational(a[0]) + I*Rational(a[1]) for a in A]
            A = [A[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            A = spmatrix(A)

            W = [int(w) for w in W] # W is integer-valued
            W = [W[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            W = spmatrix(W)

        # third, a vector H used for homogenization
        # random if projective input
        H_size = int(lines[0])
        lines = lines[1:]
        H = lines[:H_size]
        H = [h.split(' ') for h in H] # H is complex-valued
        if num_format == INT:
            H = [Integer(h[0]) + I*Integer(h[1]) for h in H]
        elif num_format == DOUBLE:
            for j in range(len(H)):
                h = H[j]
                real,imag = h.split(' ')
                H[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            H = [Rational(h[0]) + I*Rational(h[1]) for h in H]

        H = spmatrix(H)
        lines = lines[H_size:]

        # fourth, a number homVarConst
        # 0 for affine, random for projective
        hvc = lines[0].split(' ')
        if num_format == INT:
            hvc = Integer(hvc[0]) + I*Integer(hvc[1])
        elif num_format == DOUBLE:
            real,imag = hvc
            hvc = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            hvc = Rational(hvc[0]) + I*Rational(hvc[1])

        lines = lines[1:]

        # fifth, matrix B for linear slice coefficients
        num_rows, num_cols = lines[0].split(' ')
        num_rows, num_cols = int(num_rows), int(num_cols)
        B_size = num_rows*num_cols
        lines = lines[1:]

        if B_size == 0:
            B = None
        else:
            B = lines[:B_size]
            lines = lines[B_size:]

            B = [b.split(' ') for b in B] # B is complex-valued
            if num_format == INT:
                B = [Integer(b[0]) + I*Integer(b[1]) for b in B]
            elif num_format == DOUBLE:
                for j in range(len(B)):
                    real,imag = B[j]
                    B[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                B = [Rational(b[0]) + I*Rational(b[1]) for b in B]
            B = [B[j:j+num_cols] for j in range(0,B_size,num_cols)]
            B = spmatrix(B)

        # sixth and finally, vector p for patch coefficients
        p_size = int(lines[0])
        lines = lines[1:]

        p = lines[:p_size]
        p = [q.split(' ') for q in p]
        if num_format == INT:
            p = [Integer(q[0]) + I*Integer(q[1]) for q in p]
        elif num_format == DOUBLE:
            for j in range(len(p)):
                real,imag = p[j]
                p[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            p = [Rational(q[0]) + I*Rational(q[1]) for q in p]

        p = spmatrix(p)
        codims[i]['A'] = A
        codims[i]['W'] = W
        codims[i]['H'] = H
        codims[i]['homVarConst'] = hvc
        codims[i]['slice'] = B
        codims[i]['p'] = p

    return codims