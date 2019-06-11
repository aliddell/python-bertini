from __future__ import print_function

import os.path as op
import sys
from fractions import Fraction

from typing import List

import numpy as np

from naglib.exceptions import UnclassifiedException


def _line_to_complex(line: str) -> np.complex:
    real, imag = line.split(" ")
    return np.complex(float(real), float(imag))


def read_points(filename: str, tol: float=None) -> np.ndarray:
    """Read points from an output file.

    Parameters
    ----------
    filename : str
        Path to file to read.
    tol : float, optional
        If given, numbers smaller than this in absolute value will be set to 0.
        Otherwise, they will be left as is.

    Returns
    -------
    points : ndarray
        An m-by-n array of points, where m is the number of dims and n is the
        number of points.
    """

    if not op.isfile(filename):
        raise IOError(f"Can't find {filename}")

    with open(filename, "r") as io:
        lines = [l.strip() for l in io.readlines() if l != "\n"]

    n_points = int(lines[0])
    if n_points == 0:
        return np.zeros(0, dtype=np.complex) # return an empty array

    lines = lines[1:]
    n_lines = len(lines)
    n_dims = n_lines // n_points

    # TODO: peek ahead at the first point to get an idea of the data type
    # for now just assume float64
    points = np.zeros((n_dims, n_points), dtype=np.complex)

    for i in range(0, n_lines, n_dims):
        point_lines = lines[i:i+n_dims]
        points[:, i//n_dims] = np.array([_line_to_complex(p) for p in point_lines])

    if tol is not None:
        points.imag[np.abs(points.imag) < tol] = 0

    return points


def read_witness_data(filename: str):
    """Read a witness_data file.

    Parameters
    ----------
    filename : str
        Path to file to read.
    """

    if not op.isfile(filename):
        raise IOError(f"Can't find {filename}")

    with open(filename, "r") as io:
        lines = [l.strip() for l in io.readlines() if l != "\n"]

    n_dims, nonempty_codims = int(lines[0]), int(lines[1])
    # previous line includes additional homogenizing variable(s),
    # as appropriate
    lines = lines[2:]

    # following block represents a single codim; repeated for each
    codims = []

    for i in range(nonempty_codims):
        codim = int(lines[0])
        codims.append({"codim": codim})
        n_points = int(lines[1])
        lines = lines[2:]

        # following block represents a single point; repeated for each
        pts = []

        for j in range(n_points):
            prec = int(lines[0])
            lines = lines[1:]

            coords = np.array([_line_to_complex(p) for p in lines[:n_dims]])
            lines = lines[n_dims:]

            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines[0])
            lines = lines[1:]

            approx_pt = np.array([_line_to_complex(p) for p in lines[:n_dims]])
            lines = lines[n_dims:]

            condition_number = float(lines[0])
            corank = int(lines[1]) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines[2])
            largest_zero_singular = float(lines[3])
            pt_type = int(lines[4])
            multiplicity = int(lines[5])
            component_number = int(lines[6])

            if component_number == -1:
                raise UnclassifiedException(f"components in {filename} have unclassified points")

            deflations = int(lines[7])
            lines = lines[8:]

            pts.append({"coordinates": coords,
                        "corank": corank,
                        "condition number": condition_number,
                        "smallest nonzero": smallest_nonzero_singular,
                        "largest zero": largest_zero_singular,
                        "type": pt_type,
                        "multiplicity": multiplicity,
                        "component number": component_number,
                        "deflations": deflations,
                        "precision": prec,
                        "last approximation": approx_pt})
        codims[-1]["points"] = pts

    # -1 designates the end of witness points
    assert int(lines[0]) == -1
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
        n_rows, n_cols = [int(l) for l in lines[0].split()]
        AW_size = n_rows*n_cols

        lines = lines[1:]

        if AW_size == 0:
            A = None
            W = None
        else:
            A = lines[:AW_size]
            lines = lines[AW_size:]

            W = lines[:AW_size]
            lines = lines[AW_size:]

            # A is complex-valued
            if num_format == RATIONAL:
                A = [a.split() for a in A]
                A = np.array([np.complex(Fraction(a[0]), Fraction(a[1])) for a in A])
            else:
                A = np.array([_line_to_complex(a) for a in A])

            A = A.reshape((n_rows, n_cols))

            W = np.array([int(w) for w in W]).reshape((n_rows, n_cols)) # W is integer-valued

        ## WIP: PICK UP HERE

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

        H = Matrix(H)
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
        n_rows, n_cols = lines[0].split(' ')
        n_rows, n_cols = int(n_rows), int(n_cols)
        B_size = n_rows*n_cols
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
            B = [B[j:j+n_cols] for j in range(0,B_size,n_cols)]
            B = Matrix(B)

        # sixth and finally, vector p for patch coefficients
        p_size = int(lines[0])
        lines = lines[1:]

        p = lines[:p_size]
        lines = lines[p_size:]

        p = [q.split(' ') for q in p]
        if num_format == INT:
            p = [Integer(q[0]) + I*Integer(q[1]) for q in p]
        elif num_format == DOUBLE:
            for j in range(len(p)):
                real,imag = p[j]
                p[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
        elif num_format == RATIONAL:
            p = [Rational(q[0]) + I*Rational(q[1]) for q in p]

        p = Matrix(p)
        codims[i]['A'] = A
        codims[i]['W'] = W
        codims[i]['H'] = H
        codims[i]['homVarConst'] = hvc
        codims[i]['slice'] = B
        codims[i]['p'] = p

    return codims


def write_points(points: List[np.ndarray], filename: str="") -> None:
    """Print a set of points in Bertini output fashion, optionally to a file.

    Parameters
    ----------
    Keyword arguments:
    points : iterable
        Collection of points to print.
    filename : str, optional

    """

    if filename:
        io = open(filename, 'w')
    else:
        io = sys.stdout

    n_points = len(points)
    print(f"{n_points}\n", file=io)

    for p in points:
        real, imag = p.real, p.imag

        for i in range(real.size):
            print(f"{real[i]} {imag[i]}", file=io)

        print("", file=io)

    if filename:
        io.close()
