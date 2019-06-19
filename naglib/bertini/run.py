# -*- coding: utf-8 -*-
import os
import os.path as op
import multiprocessing
import subprocess
import sys
import tempfile

import numpy as np

from typing import Union

from naglib.bertini.input_file import BertiniConfig, BertiniInput
from naglib.bertini.io import (read_input_file, read_witness_data_file,
                               write_input_file, write_points_file)
from naglib.exceptions import BertiniError, NoBertiniException


def _which(exe: str) -> Union[str, None]:
    if sys.platform == "win32":
        which_exe = "where.exe"
    else:
        which_exe = "which"

    try:
        path = subprocess.check_output([which_exe, exe], stderr=subprocess.DEVNULL).splitlines()[0].decode()
    except subprocess.CalledProcessError:
        path = None

    return path


class BertiniRun(object):
    def __init__(self, config, inputs, **kwargs):
        """Construct a BertiniRun.

        Parameters
        ----------
        config : BertiniConfig
            Options to pass to Bertini.
        inputs : BertiniInput
            Bertini input section.
        **kwargs : Keyword arguments
            Arguments required for specific run types.
        """

        self.config = config

        if not isinstance(inputs, BertiniInput):
            raise TypeError("inputs must be an instance of BertiniInput")
        self._inputs = inputs

        # can use MPI for these types of run
        self._parallel = self.tracktype in (config.TZERODIM, config.TPOSDIM, config.TREGENEXT)

        if config.needs_start_points():
            if "start" not in kwargs:
                raise ValueError("specify a point or points to evaluate with the keyword argument 'start'")
            self.start = kwargs["start"]

        if config.needs_component():
            if "component" not in kwargs:
                raise ValueError("specify a component with the keyword argument 'component'")
            self._component = kwargs["component"]

        if config.needs_sample_count():
            if "sample" not in kwargs:
                raise ValueError("specify how many points to sample with the keyword argument 'sample'")
            self._sample = kwargs["sample"]

        if config.needs_projection_variables():
            if "projection" not in kwargs:
                raise ValueError("specify the variables onto which you wish to project with the keyword argument 'projection'")
            self._projection = kwargs["projection"]

        # a setting of ParameterHomotopy > 0 requires parameters
        if config.parameterhomotopy > 0 and not inputs.parameter:
            raise ValueError("you are attempting a parameter homotopy with no parameters")

        # conversely, the presence of parameters requires ParameterHomotopy > 0
        if inputs.parameter and config.parameterhomotopy == 0:
            raise ValueError("your system has parameters but you have not specified a parameter homotopy")

        # parameterhomotopy:2 requires start and final params
        if config.parameterhomotopy == 2:
            if "start_parameters" in kwargs:
                self.start_parameters = kwargs["start_parameters"]
            else:
                raise ValueError("you have selected parameterhomotopy:2 but you have not given start parameters")
            if "final_parameters" in kwargs:
                self.final_parameters = kwargs["final_parameters"]
            else:
                raise ValueError("you have selected parameterhomotopy:2 but you have not given final parameters")

        if "bertini_path" in kwargs:
            if not op.isfile(kwargs["bertini_path"]):
                raise OSError(f"didn't find Bertini at '{kwargs['bertini_path']}'")
            self._bertini = kwargs["bertini_path"]
        else:
            bertini = _which("bertini")
            if bertini is None:
                raise OSError("couldn't find a bertini executable and you didn't specify one")
            self._bertini = bertini

        if "mpi_path" in kwargs:
            if not op.isfile(kwargs["mpi_path"]):
                raise OSError(f"didn't find MPI executable at '{kwargs['mpi_path']}'")
            self._mpi = kwargs["mpi_path"]
        else:
            self._mpi = _which("mpirun")

        self._complete = False

    def _recover_components(self, witness_data):
        """
        """
        from sympy import sympify
        from naglib.core.algebra import LinearSlice
        from naglib.core.base import AffinePoint, ProjectivePoint
        from naglib.core.geometry import IrreducibleComponent
        from naglib.core.witnessdata import WitnessPoint, WitnessSet
        system = self._system
        variables = system.variables
#        if system.homvar:
#            proj_dim = len(variables)
            #homsys = system
#        else:
#            proj_dim = len(system.variables) + 1
            #homsys = system.homogenize()
        homvar = sympify('_homvar')
        while homvar in system.variables:
            homvar = sympify('_' + str(homvar))
        homsys  = system.homogenize(homvar)
        homvars = homsys.variables

        components = []

        for c in witness_data:
            codim       = c['codim']
            homVarConst = c['homVarConst']
            points      = c['points']
            coeffs      = c['slice']
            rand_mat    = c['A']
            homog_mat   = c['W']
            homog_vec   = c['H']
            hvc         = c['homVarConst']
            patch_coeff = c['p']

            comp_isprojective = homVarConst == 0

            hslice = None
            if coeffs:
                if comp_isprojective:
                    hslice = LinearSlice(coeffs, homvars, homvar)
                    if not system.homvar:
                        lslice = hslice.dehomogenize()
                else:
                    lslice = LinearSlice(coeffs, variables)
            else:
                lslice = None

            dim_list = {}

            hcoord = None
            for point in points:
                comp_id = point['component number']
                if comp_isprojective:
                    hcoord = point['coordinates']
                    if not system.homvar:
                        coord = ProjectivePoint(hcoord).dehomogenize()
                else:
                    coord = AffinePoint(point['coordinates'])

                wpoint = WitnessPoint(coord, comp_id,
                                      corank=point['corank'],
                                      condition_number=point['condition number'],
                                      smallest_nonzero=point['smallest nonzero'],
                                      largest_zero=point['largest zero'],
                                      point_type=point['type'],
                                      multiplicity=point['multiplicity'],
                                      deflations=point['deflations'],
                                      precision=point['precision'],
                                      last_approximation=point['last approximation'],
                                      homogeneous_coordinates=hcoord)

                if comp_id not in dim_list:
                    dim_list[comp_id] = []

                dim_list[comp_id].append(wpoint)

            for comp_id in dim_list.keys():
                ws = WitnessSet(system.copy(),
                                lslice,
                                dim_list[comp_id],
                                witness_data,
                                homogeneous_slice=hslice)
                component = IrreducibleComponent(ws, codim, comp_id,
                                                 randomization_matrix=rand_mat,
                                                 homogenization_matrix=homog_mat,
                                                 homogenization_vector=homog_vec,
                                                 homogenization_variable=hvc,
                                                 patch_coefficients=patch_coeff)

                components.append(component)

        return components

    def _recover_data(self):
        """
        recover the information pertinent to a run
        """

        if not self._complete:
            return

        dirname = self._dirname
        system = self._system
        tol = self._tol
        tracktype = self._tracktype
        main_data = dirname + '/main_data'
        fh = open(main_data, 'r')
        self._main_data = striplines(fh.readlines())
        fh.close()

        projective = not not system.homvar

        if tracktype == self.TEVALP:
            pass
        elif tracktype == self.TEVALPJ:
            pass
        elif tracktype == self.TNEWTP:
            pass
        elif tracktype == self.TNEWTPJ:
            pass
        elif tracktype == self.TZERODIM:
            finites = dirname + '/finite_solutions'
            startp = dirname + '/start_parameters'
            finite_solutions = read_points_file(finites, tol=tol, projective=projective)

            ptype  = self._parameter_homotopy['arg']
            if ptype == 1:
                start_parameters = read_points_file(startp, tol=tol, projective=projective)
                return finite_solutions, start_parameters

            return finite_solutions
        elif tracktype == self.TPOSDIM:
            wdfile = dirname + '/witness_data'
            self._witness_data = self.read_witness_data_file(wdfile)
            components = self._recover_components(self._witness_data)

            return components
        elif tracktype == self.TSAMPLE:
            wdfile = dirname + '/witness_data'
            self._witness_data = self.read_witness_data_file(wdfile)
            samplef = dirname + '/sampled'
            sampled = read_points_file(samplef, tol=tol, projective=projective)

            return sampled
        elif tracktype == self.TMEMTEST:
            from sympy import zeros

            wdfile = dirname + '/witness_data'
            self._witness_data = self.read_witness_data_file(wdfile)
            inmat = dirname + '/incidence_matrix'
            fh = open(inmat, 'r')
            lines = striplines(fh.readlines())
            fh.close()

            testcodim = self._component.codim
            testcid   = 0 # component_id should be 0 after write

            testp = self._start
            if type(testp) not in (list, tuple):
                testp = [testp]

            nonempty_codims = int(lines[0])
            lines = lines[1:]
            # gather nonempty codims with component count for each
            ccounts = lines[:nonempty_codims]
            lines = lines[nonempty_codims:]

            ccounts = [tuple([int(d) for d in c.split(' ')]) for c in ccounts]
            # ordered list of codims with component ids, for matrix
            cids = [(c[0], j) for c in ccounts for j in range(c[1])]
            colcount = len(cids)
            dex = cids.index((testcodim, testcid))

            numpoints = int(lines[0])
            lines = lines[1:]

            inmat = zeros(numpoints, colcount)
            # populate incidence matrix
            for i in range(numpoints):
                line = lines[i].split(' ')
                row = [int(l) for l in line]
                for j in range(colcount):
                    inmat[i,j] = row[j]

            if numpoints == 1:
                return inmat[0, dex] == 1
            else:
                ret = []
                for i in range(numpoints):
                    ret.append(inmat[i, dex] == 1)
                return ret

        elif tracktype == self.TPRINTWS:
            pointsfile = dirname + '/points.out'
            #sysfile    = dirname + '/sys.out'

            points = read_points_file(pointsfile, tol=tol, projective=projective)
            #TODO: parse linear system file and return a LinearSystem

            return points
        elif tracktype == self.TPROJECT:
            #TODO: implement
            pass
        elif tracktype == self.TISOSTAB:
            config = self._config
            ckeys = config.keys()
            lkeys = [k.lower() for k in ckeys]
            cws = None
            if 'constructwitnessset' in lkeys:
                for k in ckeys:
                    if k.lower() == 'constructwitnessset':
                        cws = k
                        break
            if cws and config[cws] == 1:
                wdfile = dirname + '/witness_data'
                self._witness_data = self.read_witness_data_file(wdfile)
                components = self._recover_components(self._witness_data)
                return components

            #TODO: read isosingular_summary and maybe output_isosingular
        elif tracktype == self.TREGENEXT:
            wdfile = dirname + '/witness_data'
            self._witness_data = self.read_witness_data_file(wdfile)
            components = self._recover_components(self._witness_data)

            return components

    def _recover_input(self):
        """
        Reads main_data and recovers the input file
        needed to reproduce a run
        """
        filename = self._dirname + "/main_data"
        key = "*************** input file needed to reproduce this run ***************\n"
        try:
            fh = open(filename, 'r')
            lines = fh.readlines()
            fh.close()
        except IOError:
            lines = []

        try:
            dex = lines.index(key)
        except ValueError:
            msg = 'no main_data file!'
            raise BertiniError(msg)
        inlines = lines[dex+1:]
        while inlines[0] == '\n':
            inlines = inlines[1:]

        dex = inlines.index('END;\n')
        inlines = inlines[:dex+1]

        return inlines

    def _write_files(self):
        from naglib.bertini.fileutils import fprint

        tracktype = self._tracktype
        dirname = self._dirname
        system = self._system
        if not op.exists(dirname):
            os.mkdir(dirname)

        ### write the system
        sysconfig = self._config.copy()
        sysconfig.update({'TrackType':tracktype})
        inputf = self._write_system(system, config=sysconfig)

        ### write out `start', `start_parameters', `final_parameters'
        if '_start' in dir(self):
            start = self._start
            if self._tracktype == self.TMEMTEST:
                startfile = dirname + '/member_points'
            else:
                startfile = dirname + '/start'
            fprint(start, startfile)
        if self._parameter_homotopy:
            phtpy = self._parameter_homotopy
            pkeys = phtpy.keys()
            if 'start parameters' in pkeys:
                startp = phtpy['start parameters']
                startpfile = dirname + '/start_parameters'
                fprint(startp, startpfile)
            if 'final parameters' in pkeys:
                finalp = phtpy['final parameters']
                finalpfile = dirname + '/final_parameters'
                fprint(finalp, finalpfile)

        ### write out component information
        if '_component' in dir(self):
            component = self._component
            #cid = component.component_id
            dim = component.dim
            if tracktype == self.TREGENEXT:
                self._write_system(component.system, 'iold', {'TrackType':1})
                witness_data = component._construct_witness_data()
                self._write_witness_data(witness_data, dirname, filename='wdold')
                instructions = ['1', 'iold', 'wdold', str(dim), '0']
                self._write_instructions(instructions)
            elif tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTWS, self.TPROJECT):
                witness_data = component._construct_witness_data()
                self._write_witness_data(witness_data, dirname)
                if tracktype == self.TSAMPLE:
                    sample = self._sample
                    instructions = [str(dim), '0', str(sample), '0', 'sampled']
                    self._write_instructions(instructions)
                elif tracktype == self.TPRINTWS:
                    instructions = [str(dim), '0', 'points.out', 'sys.out']
                    self._write_instructions(instructions)
                elif tracktype == self.TPROJECT:
                    instructions = [str(dim), '0']
                    self._write_instructions(instructions)

                    ### write out projection information
                    projection = self._projection
                    projnum = ['1' if x in projection else '0' for x in system.variables]
                    projnum = [' '.join(projnum)]
                    self._write_instructions(projnum, 'projection')

        return inputf

    def _write_instructions(self, lines, filename='instructions'):
        filename = self._dirname + '/' + filename
        fh = open(filename, 'w')
        for line in lines:
            fh.write(line + '\n')
        fh.close()

    def run(self, rerun_on_fail=False):
        # in case the user has changed any of these
        if self._parallel and self._mpi is not None:
            cmd = self._mpi
            arg = [cmd, '-np', str(multiprocessing.cpu_count()), self._bertini]
        else:
            arg = [self._bertini]

        self.setup()  # write files

        if op.isfile(op.join(self.dirname, "/instructions")):
            stdin = op.join(self.dirname, "/instructions")
        else:
            stdin = None

        arg += ["input"]

        os.chdir(self.dirname)
        if stdin is not None:
            stdin = open(stdin, "r")

        try:
            output = subprocess.check_output(arg, stdin=stdin, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            msg = naglib.system.proc_err_output(e.output)
            raise BertiniError(msg)

        if stdin is not None:
            stdin.close()

        self._complete = True
        self._output = output

    def setup(self, dirname: str = None):
        """

        Parameters
        ----------
        dirname : str
            Path to directory to run Bertini.

        Returns
        -------

        """
        # write input file
        if dirname is None:
            self._dirname = tempfile.mkdtemp()
        elif not op.isdir(dirname):
            os.makedirs(dirname)
            self._dirname = dirname
        else:
            self._dirname = dirname

        input_file = op.join(self._dirname, "input")
        write_input_file(self.config, self.inputs, input_file)

        if self.config.parameterhomotopy == 2:
            write_points_file(self.start_parameters, op.join(self._dirname, "start_parameters"))
            write_points_file(self.final_parameters, op.join(self._dirname, "final_parameters"))

    @property
    def bertini(self):
        return self._bertini

    @bertini.setter
    def bertini(self, val):
        self._bertini = val

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, val):
        if not isinstance(val, BertiniConfig):
            raise TypeError("config must be an instance of BertiniConfig")
        self._config = val

    @property
    def complete(self):
        return self._complete

    @property
    def dirname(self):
        return self._dirname

    @dirname.setter
    def dirname(self, name):
        self._dirname = name

    @property
    def inputs(self):
        return self._inputs

    @inputs.setter
    def inputs(self, val):
        if not isinstance(val, BertiniInput):
            raise TypeError("inputs must be an instance of BertiniInput")
        self._inputs = val

    @property
    def output(self):
        return self._output

    @property
    def start(self):
        if hasattr(self, "_start"):
            return self._start
        else:
            return np.zeros(0, dtype=np.complex)

    @start.setter
    def start(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError("expected a numpy array")
        if val.ndim == 1:
            val = val.reshape(val.size, 1)
        if val.shape[0] != self.config.ndims:
            raise ValueError(f"expected points of dimension {self.config.ndims} but you specified {val.shape[0]}")

        self._start = val.as_type(np.complex)

    @property
    def start_parameters(self):
        if hasattr(self, "_start_parameters"):
            return self._start_parameters
        else:
            return np.zeros(0, dtype=np.complex)

    @start_parameters.setter
    def start_parameters(self, val):
        if self.config.parameterhomotopy != 2:
            raise ValueError("only specify start_parameters for parameterhomotopy:2")
        if not isinstance(val, np.ndarray):
            raise TypeError("expected a numpy array")
        if val.size != len(self.inputs.parameter):
            raise ValueError(f"expected {len(self.inputs.parameter)} parameters but you specified {val.size}")

        self._start_parameters = val.as_type(np.complex)

    @property
    def final_parameters(self):
        if hasattr(self, "_final_parameters"):
            return self._final_parameters
        else:
            return np.zeros(0, dtype=np.complex)

    @final_parameters.setter
    def final_parameters(self, val):
        if self.config.parameterhomotopy != 2:
            raise ValueError("only specify final_parameters for parameterhomotopy:2")
        if not isinstance(val, np.ndarray):
            raise TypeError("expected a numpy array")
        if val.size != len(self.inputs.parameter):
            raise ValueError(f"expected {len(self.inputs.parameter)} parameters but you specified {val.size}")

        self._final_parameters = val.as_type(np.complex)

    @property
    def tracktype(self):
        return self.config.tracktype


def parameter_homotopy(system: list):
    pass
