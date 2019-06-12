from collections import deque
from fractions import Fraction
import os.path as op
import sys
import re

from typing import List, Tuple

import numpy as np

#from naglib.exceptions import UnclassifiedException


def _line_to_complex(line: str) -> np.complex:
    real, imag = line.split(" ")
    return np.complex(float(real), float(imag))


def _line_to_complex_rat(line: str) -> np.complex:
    real, imag = line.split(" ")
    return np.complex(Fraction(real), Fraction(imag))


def read_input_file(input_file: str) -> Tuple:
    """Parse input file.

    Parameters
    ----------
    input_file : str
        Path to input file.
        
    Returns
    -------
    configs : dict
        Key-value pairs of parameters set in the CONFIG section.
    inputs : dict
        Values set in the INPUT section.
    misclines : list
        This will go away.
    """
    if not op.isfile(op.abspath(input_file)):
        raise IOError(f"Input file '{input_file}' does not exist")

    with open(input_file, "r") as fh:
        lines = deque([l.strip() for l in fh.readlines() if l != "\n"])

    configs = {}
    inputs = {"variable_group": deque(),
              "variable": deque(),
              "hom_variable_group": deque(),
              "pathvariable": deque(),
              "random": deque(),
              "constant": {},
              "function": {},
              "parameter": {}}
    misclines = []

    in_config = in_input = False
    vargroup_re = re.compile(r"^variable_group\s+", re.I)
    var_re = re.compile(r"^variable\s+", re.I)
    homvargroup_re = re.compile(r"^hom_variable_group\s+", re.I)
    pathvar_re = re.compile(r"^pathvariable\s+", re.I)
    random_re = re.compile(r"^random\s+", re.I)
    constant_re = re.compile(r"^constant\s+", re.I)
    function_re = re.compile(r"^function\s+", re.I)
    parameter_re = re.compile(r"^parameter\s+", re.I)

    while lines:
        line = lines.popleft().strip(" ;")

        if line.lower() == "config":
            in_config = True
            continue
        elif line.lower() == "input":
            in_input = True
            continue
        elif line.lower() == "end":
            in_input = in_config = False
            continue
        elif line.startswith("%"):
            continue
        
        if in_config:
            key, val = map(lambda x: x.strip(), line.split(":"))
            val = val.split("%")[0].strip(" ;") # remove comment, semicolon, trailing whitespace

            if key.lower() in ["condnumthreshold", "endpointfinitethreshold", "finaltol",
                               "pathtruncationthreshold", "samplefactor", "securitymaxnorm",
                               "slicetolbeforeeg", "slicetolduringeg"]:
                val = float(val)
            elif key.lower() in ["coeffbound", "degreebound", "maxcodimension", "mptype",
                                 "parameterhomotopy", "precision", "randomseed", "securitylevel",
                                 "sharpendigits", "sharpenonly", "specificcodimension",
                                 "tracktype", "useregeneration", "userhomotopy", "witnessgentype",
                                 "witnesssupersetonly"]:
                val = int(val)

            configs[key] = val

        elif in_input:
            if vargroup_re.match(line):
                line = vargroup_re.sub("", line)
                inputs["variable_group"].append(re.split(r",\s*", line))
            elif var_re.match(line):
                line = var_re.sub("", line)
                inputs["variable"].append(re.split(r",\s*", line))
            elif homvargroup_re.match(line):
                line = homvargroup_re.sub("", line)
                inputs["hom_variable_group"].append(re.split(r",\s*", line))
            elif pathvar_re.match(line):
                line = pathvar_re.sub("", line)
                inputs["pathvariable"].append(re.split(r",\s*", line))
            elif random_re.match(line):
                line = random_re.sub("", line)
                inputs["random"].append(re.split(r",\s*", line))
            elif constant_re.match(line):
                line = constant_re.sub("", line)
                constants = re.split(r",\s*", line)
                for c in constants:
                    inputs["constant"][c] = None
            elif function_re.match(line):
                line = function_re.sub("", line)
                functions = re.split(r",\s*", line)
                for f in functions:
                    inputs["function"][f] = None
            elif parameter_re.match(line):
                line = parameter_re.sub("", line)
                params = re.split(r",\s*", line)
                for p in params:
                    inputs["parameter"][p] = None
            else:
                terms = re.split(r";\s*", line)
                for term in terms:
                    # remove comments
                    term = term.split("%")[0].strip()
                    # split by =
                    term = re.split(r"\s*=\s*", term)
                    if len(term) != 2:
                        misclines.append("=".join(term))
                    else:
                        term, val = term
                        if term in inputs["constant"]:
                            inputs["constant"][term] = val
                        elif term in inputs["function"]:
                            inputs["function"][term] = val
                        elif term in inputs["parameter"]:
                            inputs["parameter"][term] = val
                        else:
                            misclines.append(term + "=" + val)
    
    return configs, inputs, misclines


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

    with open(filename, "r") as fh:
        lines = [l.strip() for l in fh.readlines() if l != "\n"]

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

    with open(filename, "r") as fh:
        lines = deque([l.strip() for l in fh.readlines() if l != "\n"])

    n_dims, nonempty_codims = int(lines.popleft()), int(lines.popleft())
    # previous line includes additional homogenizing variable(s),
    # as appropriate

    # following block represents a single codim; repeated for each
    codims = []

    for i in range(nonempty_codims):
        codim, n_points = int(lines.popleft()), int(lines.popleft())
        codims.append({"codim": codim})

        # following block represents a single point; repeated for each
        pts = []

        for j in range(n_points):
            prec = int(lines.popleft())

            coords = np.array([_line_to_complex(p) for p in lines[:n_dims]])
            lines = lines[n_dims:]

            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines.popleft())

            approx_pt = np.array([_line_to_complex(p) for p in lines[:n_dims]])
            lines = lines[n_dims:]

            condition_number = float(lines.popleft())
            corank = int(lines.popleft()) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines.popleft())
            largest_zero_singular = float(lines.popleft())
            pt_type = int(lines.popleft())
            multiplicity = int(lines.popleft())
            component_number = int(lines.popleft())

            if component_number == -1:
                raise UnclassifiedException(f"components in {filename} have unclassified points")

            deflations = int(lines.popleft())

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
    assert int(lines.popleft()) == -1

    RATIONAL = 2

    # remaining data is related to slices, randomization,
    # homogenization, and patches
    num_format = int(lines.popleft())

    # previous line describes format for remainder of data

    # the following block is repeated for each nonempty codim.
    # first, matrix A used for randomization
    # second, matrix W
    for i in range(nonempty_codims):
        n_rows, n_cols = [int(l) for l in lines.popleft().split(" ")]
        AW_size = n_rows*n_cols

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
                A = np.array([_line_to_complex_rat(a) for a in A])
            else:
                A = np.array([_line_to_complex(a) for a in A])

            A = A.reshape((n_rows, n_cols))

            W = np.array([int(w) for w in W]).reshape((n_rows, n_cols)) # W is integer-valued

        # third, a vector H used for homogenization
        # random if projective input
        H_size = int(lines.popleft())
        H = lines[:H_size]
        lines = lines[H_size:]

        # H is complex-valued
        if num_format == RATIONAL:
            H = np.array([_line_to_complex_rat(h) for h in H])
        else:
            H = np.array(_line_to_complex(h) for h in H)

        # fourth, a number homVarConst
        # 0 for affine, random for projective
        hvc = lines.popleft()
        if num_format == RATIONAL:
            hvc = _line_to_complex_rat(hvc)
        else:
            hvc = _line_to_complex(hvc)

        # fifth, matrix B for linear slice coefficients
        n_rows, n_cols = [int(l) for l in lines.popleft().split(" ")]
        B_size = n_rows*n_cols

        if B_size == 0:
            B = None
        else:
            B = lines[:B_size]
            lines = lines[B_size:]

            # B is complex-valued
            if num_format == RATIONAL:
                B = np.array([_line_to_complex_rat(b) for b in B])
            else:
                B = np.array([_line_to_complex(b) for b in B])

            B = B.reshape((n_rows, n_cols))

        # sixth and finally, vector p for patch coefficients
        p_size = int(lines.popleft())

        patch = lines[:p_size]
        lines = lines[p_size:]

        # patch is complex-valued
        if num_format == RATIONAL:
            patch = np.array([_line_to_complex_rat(p) for p in patch])
        else:
            patch = np.array([_line_to_complex(p) for p in patch])

        codims[i]["A"] = A
        codims[i]["W"] = W
        codims[i]["H"] = H
        codims[i]["homVarConst"] = hvc
        codims[i]["slice"] = B
        codims[i]["patch"] = patch

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
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    n_points = len(points)
    print(f"{n_points}\n", file=fh)

    for p in points:
        real, imag = p.real, p.imag

        for i in range(real.size):
            print(f"{real[i]} {imag[i]}", file=fh)

        print("", file=fh)

    if filename:
        fh.close()

read_witness_data("C:/Users/Alan/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/alan/BertiniLinux64_v1.6/examples/pos_dim/basic_pos_dim/witness_data")