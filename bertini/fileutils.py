from __future__ import print_function

from fractions import Fraction as fraction
from os.path import isfile
from sys import stderr, stdout
import re

from mpmath import mpc, matrix

from naglib.fileutils import striplines

def parse_witness_data(filename):
    if not isfile(filename):
        return None

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
        num_points = int(lines[1])
        lines = lines[2:]
        # following block represents a single point; repeated for each
        pts = []
        for j in range(num_points):
            prec = int(lines[0])

            lines = lines[1:]
            pt = []
            for k in range(num_vars):
                coord = lines[k].split(' ')
                pt.append(mpc(*coord))
            pt = tuple(pt)
            lines = lines[num_vars:]
            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines[0])
            lines = lines[1:]
            approx_pt = []
            for k in range(num_vars):
                coord = lines[k].split(' ')
                approx_pt.append(mpc(*coord))

            lines = lines[num_vars:]
            condition_number = float(lines[0])
            corank = float(lines[1]) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines[2])
            largest_zero_singular = float(lines[3])
            pt_type = int(lines[4])
            multiplicity = int(lines[5])
            component_number = int(lines[6])
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
        codims.append((codim,pts))

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

        lines = lines[:2]

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
                A = [mpc(int(a[0]), int(a[1])) for a in A]
            elif num_format == DOUBLE:
                A = [mpc(float(a[0]), float(a[1])) for a in A]
            elif num_format == RATIONAL:
                A = [mpc(float(fraction(a[0])), float(fraction(a[1]))) for a in A]
            A = [A[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            A = matrix(A)

            W = [int(w) for w in W] # W is integer-valued
            W = [W[j:j+num_cols] for j in range(0,AW_size,num_cols)]
            W = matrix(W)

        # third, a vector H used for homogenization
        # random if projective input
        H_size = int(lines[0])
        lines = lines[1:]
        H = lines[:H_size]
        H = [h.split(' ') for h in H] # H is complex-valued
        if num_format == INT:
            H = [mpc(int(h[0]), int(h[1])) for h in H]
        elif num_format == DOUBLE:
            H = [mpc(float(h[0]), float(h[1])) for h in H]
        elif num_format == RATIONAL:
            H = [mpc(float(fraction(h[0])), float(fraction(h[1]))) for h in H]

        H = matrix(H)
        lines = lines[H_size:]

        # fourth, a number homVarConst
        # 0 for affine, random for projective
        hvc = lines[0].split(' ')
        if num_format == INT:
            hvc = mpc(int(hvc[0]), int(hvc[1]))
        elif num_format == DOUBLE:
            hvc = mpc(float(hvc[0]), float(hvc[1]))
        elif num_format == RATIONAL:
            hvc = mpc(float(fraction(hvc[0])), float(fraction(hvc[1])))

        lines = lines[1:]

        # fifth, matrix B for linear slice coefficients
        num_rows, num_cols = lines[0].split(' ')
        num_rows, num_cols = int(num_rows), int(num_cols)
        B_size = num_rows*num_cols

        if B_size == 0:
            B = None
        else:
            B = lines[:B_size]
            lines = lines[B_size:]

            B = [b.split(' ') for b in B] # B is complex-valued
            if num_format == INT:
                B = [mpc(int(b[0]), int(b[1])) for b in B]
            elif num_format == DOUBLE:
                B = [mpc(float(b[0]), float(b[1])) for b in B]
            elif num_format == RATIONAL:
                B = [mpc(float(fraction(b[0])), float(fraction(b[1]))) for b in B]
            B = [B[j:j+num_cols] for j in range(0,B_size,num_cols)]
            B = matrix(B)

        # sixth and finally, vector p for patch coefficients
        p_size = int(lines[0])
        lines = lines[1:]

        p = lines[:p_size]
        p = [q.split(' ') for q in p]
        if num_format == INT:
            p = [mpc(int(q[0]), int(q[1])) for q in p]
        elif num_format == DOUBLE:
            p = [mpc(float(q[0]), float(q[1])) for q in p]
        elif num_format == RATIONAL:
            p = [mpc(float(fraction(q[0])), float(fraction(q[1]))) for q in p]

        p = matrix(p)
        codims[i][1]['A'] = A
        codims[i][1]['W'] = W
        codims[i][1]['H'] = H
        codims[i][1]['homVarConst'] = hvc
        codims[i][1]['B'] = B
        codims[i][1]['p'] = p

    return codims


def read_input(filename):
    """Read in an input file and parse out the variables, functions, etc"""
    fh = open(filename, 'r')
    lines = striplines(fh.readlines())
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

def parselines(lines, **kwargs):
    """Returns a list of n-tuples of complex numbers, realized as ordered pairs"""
    lines = striplines(lines)
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
                    p0 = fraction(p[0]).limit_denominator(d_limit)
                    p1 = fraction(p[1]).limit_denominator(d_limit)
                elif rational:
                    p0 = fraction(p[0])
                    p1 = fraction(p[1])
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
                point = [(fraction(p[0]).limit_denominator(d_limit), fraction(p[1]).limit_denominator(d_limit)) for p in point]
            elif rational:
                point = [(fraction(p[0]), fraction(p[1])) for p in point]
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
                print('Bailing out', file=stderr)
                res = None
                return res
        elif isfile(filename):
            if kwargs['overwrite']:
                fh = open(filename, 'w')
            else:
                print('Bailing out', file=stderr)
                res = None
                return res
        else:
            fh = open(filename, 'w')
    else:
        usestdout = True
        fh = stdout

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

    points = parselines(lines,**kwargs)
    return points

def write_system(variables, params, functions, constants=None,
                 tracktype=0, phtpy=0, witness=0, filename='input'):

    # Python exponentiation: x**y; bertini exponentiation: x^y
    str_funcs = [re.sub(string=str(f), pattern=r'\*\*', repl='^') \
        for f in functions]

    fh = open(filename, 'w')
    fh.write('CONFIG\n')
    fh.write('TrackType:{0};\n'.format(tracktype))
    fh.write('ParameterHomotopy:{0};\n'.format(phtpy))
    fh.write('ConstructWitnessSet:{0};\n'.format(witness))
    fh.write('END\n')
    fh.write('INPUT\n')

    if phtpy:
        fh.write('parameter {0};\n'.format(','.join([str(p) for p in params])))
    # write out variable group
    variable_group = ','.join([str(v) for v in variables])
    fh.write('variable_group {0};\n'.format(variable_group))
    if constants:
        fh.write('constant {0};\n'.format(','.join(['rp{0}'.format(i) \
            for i in range(len(constants))])))

    function_group = ['f{0}'.format(i+1) for i in range(len(str_funcs))]
    fh.write('function {0};\n'.format(','.join([f for f in function_group])))
    if constants:
        for i in range(len(constants)):
            fh.write('rp{0} = {1};\n'.format(i, constants[i]))
    for i in range(len(str_funcs)):
        fh.write('f{0} = {1};\n'.format(i+1, str_funcs[i]))

    fh.write('\n')

    fh.write('END\n')
    fh.close()