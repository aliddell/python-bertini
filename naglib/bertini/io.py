from _io import TextIOWrapper
from collections import deque, OrderedDict
from fractions import Fraction
import os.path as op
import sys
import re

from typing import Tuple, Callable

import mpmath as mp
import numpy as np

from naglib.bertini.input_file import BertiniConfig, BertiniInput
from naglib.exceptions import UnclassifiedException


def _line_to_complex(line: str, multi: bool = False) -> complex:
    real, imag = line.split(" ")
    return mp.mpc(real, imag) if multi else np.complex(float(real), float(imag))


def _line_to_complex_rat(line: str) -> complex:
    real, imag = line.split(" ")
    # TODO: a Gaussian rational type that plays well with arr.real, arr.imag
    return np.complex(Fraction(real), Fraction(imag))


def parse_input_file(fh: TextIOWrapper, stop_if: Callable = None) -> Tuple[BertiniConfig, BertiniInput, list]:
    """Given an open file handle, read until stopping criterion is satisfied
    and try to parse an input file from it.

    Parameters
    ----------
    fh : _io.TextIOWrapper
        An open file handle.
    stop_if : function, optional
        Stop when the next line satisfies this function.

    Returns
    -------
    config : BertiniConfig
        Key-value pairs of parameters set in the CONFIG section.
    inputs : BertiniInput
        Values set in the INPUT section.
    """
    if stop_if is None:
        stop_if = lambda line: line != ""

    config = BertiniConfig()
    inputs = BertiniInput()
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

    line = fh.readline()
    while not stop_if(line):
        line = line.strip(" ;")

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
            key, val = map(lambda l: l.strip(), line.split(":"))
            val = val.split("%")[0].strip(" ;")  # remove comment, semicolon, trailing whitespace

            setattr(config, key.lower(), val)

        elif in_input:
            if vargroup_re.match(line):
                line = vargroup_re.sub("", line)
                inputs.variable_group.append(re.split(r",\s*", line))
            elif var_re.match(line):
                line = var_re.sub("", line)
                inputs.variable.append(re.split(r",\s*", line))
            elif homvargroup_re.match(line):
                line = homvargroup_re.sub("", line)
                inputs.hom_variable_group.append(re.split(r",\s*", line))
            elif pathvar_re.match(line):
                line = pathvar_re.sub("", line)
                inputs.pathvariable.append(re.split(r",\s*", line))
            elif random_re.match(line):
                line = random_re.sub("", line)
                inputs.random.append(re.split(r",\s*", line))
            elif constant_re.match(line):
                line = constant_re.sub("", line)
                constants = re.split(r",\s*", line)
                for c in constants:
                    inputs.constant[c] = None
            elif function_re.match(line):
                line = function_re.sub("", line)
                functions = re.split(r",\s*", line)
                for f in functions:
                    inputs.function[f] = None
            elif parameter_re.match(line):
                line = parameter_re.sub("", line)
                params = re.split(r",\s*", line)
                for p in params:
                    inputs.parameter[p] = None
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
                        if term in inputs.constant:
                            inputs.constant[term] = val
                        elif term in inputs.function:
                            inputs.function[term] = val
                        elif term in inputs.parameter:
                            inputs.parameter[term] = val
                        else:
                            inputs.subfunction[term] = val

        line = fh.readline()

    return config, inputs, misclines


def read_input_file(input_file: str) -> Tuple[BertiniConfig, BertiniInput, list]:
    """Given a path to an input file, parse its contents.

    Parameters
    ----------
    input_file : str
        Path to input file.

    Returns
    -------
    config : BertiniConfig
        Key-value pairs of parameters set in the CONFIG section.
    inputs : BertiniInput
        Values set in the INPUT section.
    """

    if not op.isfile(op.abspath(input_file)):
        raise IOError(f"Input file '{input_file}' not found")

    with open(input_file, "r") as fh:
        config, inputs, misclines = parse_input_file(fh)

    return config, inputs, misclines


def read_points_file(points_file: str, tol: float = None, multi: bool = False) -> np.ndarray:
    """Read points from an output file.

    Parameters
    ----------
    points_file : str
        Path to file to read.
    tol : float, optional
        If given, numbers smaller than this in absolute value will be set to 0.
        Otherwise, they will be left as is.
    multi : bool, optional
        Use a multiple precision type if true (infer precision from string length)

    Returns
    -------
    points : ndarray
        An m-by-n array of points, where m is the number of dims and n is the
        number of points.
    """

    if not op.isfile(points_file):
        raise IOError(f"Points file '{points_file}' not found")

    with open(points_file, "r") as fh:
        lines = [l.strip(";\n") for l in fh.readlines() if l != "\n"]

    n_points = int(lines.pop(0))
    if n_points == 0:
        return np.zeros(0, dtype=np.complex)  # return an empty array

    n_lines = len(lines)
    n_dims = n_lines // n_points

    if multi:  # peak ahead to get an idea of precision required
        r = lines[0].split(" ")[0]
        if "e" in r.lower():  # remove exponent
            r = r[:r.lower().index("e")]
        r = r.strip("-0").replace(".", "")

        prec_est = int(np.floor(len(r) * np.log2(10)))

        if mp.mp.prec < prec_est:  # only ever raise working precision
            mp.mp.prec = prec_est

        if mp.mp.prec == 53:  # default precision, just use doubles
            multi = False

    points = np.zeros((n_dims, n_points), dtype=mp.mpc if multi else np.complex)

    for i in range(0, n_lines, n_dims):
        point_lines = lines[i:i+n_dims]
        points[:, i//n_dims] = np.array([_line_to_complex(p, multi) for p in point_lines])

    if tol is not None:
        points.imag[np.abs(points.imag) < tol] = 0

    return points


def read_witness_data_file(witness_data_file: str) -> list:
    """Read a witness_data file.

    Parameters
    ----------
    witness_data_file : str
        Path to file to read.
    """

    if not op.isfile(witness_data_file):
        raise IOError(f"Witness data file '{witness_data_file}' not found")

    with open(witness_data_file, "r") as fh:
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
            coords = np.array([_line_to_complex(lines.popleft(), prec>52) for k in range(n_dims)])

            # the next point is the last approximation of the point
            # on the path before convergence
            prec = int(lines.popleft())
            approx_pt = np.array([_line_to_complex(lines.popleft(), prec>52) for k in range(n_dims)])

            condition_number = float(lines.popleft())
            corank = int(lines.popleft()) # corank of Jacobian at this point
            smallest_nonzero_singular = float(lines.popleft())
            largest_zero_singular = float(lines.popleft())
            pt_type = int(lines.popleft())
            multiplicity = int(lines.popleft())
            component_number = int(lines.popleft())

            if component_number == -1:
                raise UnclassifiedException(f"components in {witness_data_file} have unclassified points")

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

    MULTIPRECISION = 1
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
            # A is complex-valued
            if num_format == RATIONAL:
                A = np.array([_line_to_complex_rat(lines.popleft()) for k in range(AW_size)])
            else:
                A = np.array([_line_to_complex(lines.popleft(),
                                               num_format==MULTIPRECISION) for k in range(AW_size)])
            A = A.reshape((n_rows, n_cols))

            W = np.array([int(lines.popleft()) for k in range(AW_size)]).reshape((n_rows, n_cols)) # W is integer-valued

        # third, a vector H used for homogenization
        # random if projective input
        H_size = int(lines.popleft())

        # H is complex-valued
        if num_format == RATIONAL:
            H = np.array([_line_to_complex_rat(lines.popleft()) for k in range(H_size)])
        else:
            H = np.array(_line_to_complex(lines.popleft(),
                                          num_format==MULTIPRECISION) for k in range(H_size))

        # fourth, a number homVarConst
        # 0 for affine, random for projective
        if num_format == RATIONAL:
            hvc = _line_to_complex_rat(lines.popleft())
        else:
            hvc = _line_to_complex(lines.popleft())

        # fifth, matrix B for linear slice coefficients
        n_rows, n_cols = [int(l) for l in lines.popleft().split(" ")]
        B_size = n_rows*n_cols

        if B_size == 0:
            B = None
        else:
            # B is complex-valued
            if num_format == RATIONAL:
                B = np.array([_line_to_complex_rat(lines.popleft()) for k in range(B_size)])
            else:
                B = np.array([_line_to_complex(lines.popleft(),
                                               num_format==MULTIPRECISION) for k in range(B_size)])

            B = B.reshape((n_rows, n_cols))

        # sixth and finally, vector p for patch coefficients
        patch_size = int(lines.popleft())

        # patch is complex-valued
        if num_format == RATIONAL:
            patch = np.array([_line_to_complex_rat(lines.popleft()) for k in range(patch_size)])
        else:
            patch = np.array([_line_to_complex(lines.popleft(),
                                               num_format==MULTIPRECISION) for k in range(patch_size)])

        codims[i]["A"] = A
        codims[i]["W"] = W
        codims[i]["H"] = H
        codims[i]["homVarConst"] = hvc
        codims[i]["slice"] = B
        codims[i]["patch"] = patch

    return codims


def write_input_file(config: BertiniConfig, inputs: BertiniInput, input_file: str):
    """Write a Bertini input file.

    Parameters
    ----------
    inputs : BertiniInput
        Input values.
    config : BertiniConfig
        Config values.
    input_file : str
        Path to input file.
    """

    with open(input_file, "w") as fh:
        # config section
        print(config, file=fh)

        # input section
        print(inputs, file=fh)


def write_points_file(points: np.ndarray, points_file: str = "") -> None:
    """Print a set of points in Bertini output fashion, optionally to a file.

    Parameters
    ----------
    Keyword arguments:
    points : iterable
        Collection of points to print.
    points_file : str, optional
        Path to file to write to. Writes to stdout if not specified.
    """

    if points_file:
        fh = open(points_file, "w")
    else:
        fh = sys.stdout

    if points.ndim == 1:  # single point
        points = points.reshape(points.size, 1)

    n_dims, n_points = points.shape

    print(f"{n_points}\n", file=fh)

    for j in range(n_points):
        p = points[:, j]
        real, imag = p.real, p.imag

        for i in range(n_dims):
            print(f"{real[i]} {imag[i]}", file=fh)

        print("", file=fh)

    if points_file:
        fh.close()


def write_witness_data_file(witness_data: dict, witness_data_file: str):
        """
        """
        from sympy import Integer, Float
        fh = open(witness_data_file, "w")

        nonempty_codims = len(witness_data)
        num_vars = len(witness_data[0]['points'][0]['coordinates'])
        fh.write('{0}\n'.format(num_vars))
        fh.write('{0}\n'.format(nonempty_codims))

        for i in range(nonempty_codims):
            wd_codim = witness_data[i]
            fh.write('{0}\n'.format(wd_codim['codim']))
            codim_points = [p for p in wd_codim['points']]
            fh.write('{0}\n'.format(len(codim_points)))
            for p in codim_points:
                prec = p['precision']
                fh.write('{0}\n'.format(prec))

                coordinates = p['coordinates']
                for c in coordinates:
                    real,imag = c.as_real_imag()
                    fh.write('{0} {1}\n'.format(real, imag))

                fh.write('{0}\n'.format(prec))
                approx = p['last approximation']
                for a in approx:
                    real,imag = a.as_real_imag()
                    fh.write('{0} {1}\n'.format(real, imag))
                fh.write('{0}\n'.format(p['condition number']))
                fh.write('{0}\n'.format(p['corank']))
                fh.write('{0}\n'.format(p['smallest nonzero']))
                fh.write('{0}\n'.format(p['largest zero']))
                fh.write('{0}\n'.format(p['type']))
                fh.write('{0}\n'.format(p['multiplicity']))
                fh.write('{0}\n'.format(p['component number']))
                fh.write('{0}\n'.format(p['deflations']))
        fh.write('-1\n\n') # -1 designates the end of witness points

        h1 = witness_data[0]['H'][0]
        if type(h1) == Integer:
            numtype = 0
        elif type(h1) == Float:
            numtype = 1
        else:
            numtype = 2
        fh.write('{0}\n'.format(numtype))

        for i in range(nonempty_codims):
            wd_codim = witness_data[i]
            A   = wd_codim['A']
            W   = wd_codim['W']
            H   = wd_codim['H']
            hvc = wd_codim['homVarConst']
            B   = wd_codim['slice']
            P   = wd_codim['p']

            if A: # also W
                num_rows, num_cols = A.shape
                fh.write('{0} {1}\n'.format(num_rows, num_cols))
                for j in range(num_rows):
                    for k in range(num_cols):
                        real, imag = A[j,k].as_real_imag()
                        fh.write('{0} {1}\n'.format(real, imag))
                for j in range(num_rows):
                    for k in range(num_cols):
                        fh.write('{0}\n'.format(W[j,k])) # W is an *integer* matrix
            else:
                fh.write('1 0\n')

            fh.write('\n')
            h = len(H)
            fh.write('{0}\n'.format(h))
            for j in range(h):
                real, imag = H[j].as_real_imag()
                fh.write('{0} {1}\n'.format(real, imag))

            fh.write('\n')
            real,imag = hvc.as_real_imag()
            fh.write('{0} {1}\n'.format(real, imag))
            if B:
                num_rows, num_cols = B.shape
                fh.write('{0} {1}\n'.format(num_rows, num_cols))
                for j in range(num_rows):
                    for k in range(num_cols):
                        real, imag = B[j,k].as_real_imag()
                        fh.write('{0} {1}\n'.format(real, imag))
            else:
                fh.write('1 0\n')

            p = len(P)
            fh.write('{0}\n'.format(p))
            for j in range(p):
                real, imag = P[j].as_real_imag()
                fh.write('{0} {1}\n'.format(real, imag))

        fh.close()


def extract_error_message(output: str) -> str:
    """Extract Bertini error message.

    Parameters
    ----------
    output : str
        Bertini standard output.

    Returns
    -------
    err_message: str
        Just the relevant error text.
    """
    lines = output.splitlines()

    idx = None
    for l in lines:
        if l.startswith("ERROR"):
            idx = lines.index(l)

    if idx is None:
        return output
    else:
        lines = lines[idx:-1]  # remove "Bertini will now exit due to this error"

        # remove "ERROR: "
        lines[0] = lines[0].replace("ERROR: ", "")

        return "\n".join(lines)
