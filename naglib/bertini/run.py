# -*- coding: utf-8 -*-
from collections import OrderedDict
import os
import os.path as op
import subprocess
import tempfile

import numpy as np

from naglib.bertini.input_file import BertiniConfig
from naglib.bertini.io import read_input_file, read_witness_data_file
from naglib.constants import TOL
from naglib.system import BERTINI, MPIRUN, PCOUNT
from naglib.exceptions import BertiniError, NoBertiniException


class BertiniRun(object):
    def __init__(self, config, **kwargs):
        """Construct a BertiniRun.

        Parameters
        ----------
        config : BertiniConfig
            Options to pass to Bertini.
        **kwargs
            Keyword arguments.

            variable_group : list
                Variable group(s)
            variable : list
                Variables
            hom_variable_group : list
                Homogeneous variables
            pathvariable : list
                Path variables
            random : list
                Collection of complex-valued random variables
            random_real : list
                Collection of complex-valued random variables
            constant : dict
                Key-value pairs of constant symbols and their values
            subfunction : dict
                Key-value pairs of subfunction symbols and their values
            parameter : OrderedDict
                Key-value pairs of parameter symbols and their values
            function : OrderedDict
                Key-value pairs of function symbols and their values
        """

        if not isinstance(config, BertiniConfig):
            raise TypeError("config must be an instance of BertiniConfig")
        self._config = config

        variable_group = kwargs["variable_group"] if "variable_group" in kwargs else []
        if isinstance(variable_group, tuple):
            variable_group = list(variable_group)
        if not isinstance(variable_group, list):
            raise TypeError(f"variable_group type '{type(variable_group)}' not understood")

        variable = kwargs["variable"] if "variable" in kwargs else []
        if isinstance(variable, tuple):
            variable = list(variable)
        if not isinstance(variable, list):
            raise TypeError(f"variable type '{type(variable)}' not understood")

        hom_variable_group = kwargs["hom_variable_group"] if "hom_variable_group" in kwargs else []
        if isinstance(hom_variable_group, tuple):
            hom_variable_group = list(hom_variable_group)
        if not isinstance(hom_variable_group, list):
            raise TypeError(f"hom_variable_group type '{type(hom_variable_group)}' not understood")

        pathvariable = kwargs["pathvariable"] if "pathvariable" in kwargs else []
        if isinstance(pathvariable, tuple):
            pathvariable = list(pathvariable)
        if not isinstance(pathvariable, list):
            raise TypeError(f"pathvariable type '{type(pathvariable)}' not understood")

        random = kwargs["random"] if "random" in kwargs else []
        if isinstance(random, tuple):
            random = list(random)
        if not isinstance(random, list):
            raise TypeError(f"random type '{type(random)}' not understood")

        random_real = kwargs["random_real"] if "random_real" in kwargs else []
        if isinstance(random_real, tuple):
            random_real = list(random_real)
        if not isinstance(random_real, list):
            raise TypeError(f"random_real type '{type(random_real)}' not understood")

        constant = kwargs["constant"] if "constant" in kwargs else {}
        if not isinstance(constant, dict):
            raise TypeError(f"constant type '{type(constant)}' not understood")

        subfunction = kwargs["subfunction"] if "subfunction" in kwargs else {}
        if not isinstance(subfunction, dict):
            raise TypeError(f"subfunction type '{type(subfunction)}' not understood")

        parameter = kwargs["parameter"] if "parameter" in kwargs else OrderedDict()
        if not isinstance(parameter, dict):
            raise TypeError(f"parameter type '{type(parameter)}' not understood")

        function = kwargs["function"] if "function" in kwargs else OrderedDict()
        if isinstance(function, dict):
            function = OrderedDict(function)
        if not isinstance(function, OrderedDict):
            raise TypeError(f"function type '{type(function)}' not understood")

        self._config = dict([(k.lower(), config[k]) for k in config])

        # check tracktype, set values accordingly
        if "tracktype" not in config:
            config["tracktype"] = 0

        if config["tracktype"] not in range(-4, 8):
            raise ValueError("TrackType must be an integer between -4 and 7 (inclusive)")

        self._parallel = self.tracktype in (self.TZERODIM, self.TPOSDIM, self.TREGENEXT)

        # check to see if tracktype jives with kwargs
        msg = ""
        ## start point(s) required
        if self.tracktype in (self.TEVALP, self.TEVALPJ, self.TNEWTP, self.TNEWTPJ, self.TMEMTEST, self.TISOSTAB) and "start" not in kkeys:
            msg = "specify a point or points to evaluate with the keyword argument 'start'"
        elif "start" in kkeys:
            start = kwargs["start"]
            if type(start) not in (list, tuple):
                start = [start]
            self._start = start
        ## component required
        if tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTWS, self.TPROJECT, self.TREGENEXT) and 'component' not in kkeys:
            msg = "specify a component with the keyword argument 'component'"
        elif 'component' in kkeys:
            self._component = kwargs['component']
        ## sample count required
        if tracktype == self.TSAMPLE and 'sample' not in kkeys:
            msg = "specify how many points to sample with the keyword argument `sample'"
        elif tracktype == self.TSAMPLE:
            self._sample = kwargs['sample']
        ## projection variables required
        if tracktype == self.TPROJECT and 'projection' not in kkeys:
            msg = "specify the variables onto which you wish to project with the keyword argument `projection'"
        elif 'projection' in kkeys:
            self._projection = kwargs['projection']
        if msg:
            raise KeyError(msg)

        # tolerance
        # TODO:make this mean something
        if 'tol' in kkeys:
            self._tol = kwargs['tol']
        else:
            self._tol = TOL

        # parameter homotopy
        self._parameter_homotopy = {'key':'', 'arg':0}
        if 'parameterhomotopy' in ckeys:
            ckeys2 = config.keys()
            for k in ckeys2:
                if k.lower() == 'parameterhomotopy':
                    self._parameter_homotopy['key'] = 'ParameterHomotopy'
                    self._parameter_homotopy['arg'] = config[k] # in (0,1,2)
                    break

            del ckeys2

        # ensure the system jives with the call for parameter homotopy
        msg = ''
        if not system.parameters and self._parameter_homotopy['arg'] > 0:
            msg = "you have attempted to define a parameter homotopy on a system with no parameters!"
        elif system.parameters and self._parameter_homotopy['arg'] <= 0:
            msg = "a parameterized system requires ParameterHomotopy either 1 or 2"
        elif tracktype != self.TZERODIM and self._parameter_homotopy['arg'] > 0:
            msg = "parameter homotopy only supported for zero-dimensional runs"

        if msg:
            raise KeyError(msg)

        if 'start' in kkeys:
            start = kwargs['start']
            if type(start) not in (tuple, list):
                start = [start]
            # this doesn't go in self._parameter_homotopy because other kinds of run use start files
            self._start = start
        if 'start_parameters' in kkeys:
            startp = kwargs['start_parameters']
            if type(startp) not in (tuple, list):
                startp = [startp]
            self._parameter_homotopy['start parameters'] = startp
        else:
            startp = None
        if 'final_parameters' in kkeys:
            finalp = kwargs['final_parameters']
            if type(finalp) not in (tuple, list):
                finalp = [finalp]
            self._parameter_homotopy['final parameters'] = finalp
        else:
            finalp = None

        # if a system specifies one of start parameters or final parameters it must specify the other
        if (startp and not finalp) or (finalp and not startp):
            msg = "specify both start parameters and final parameters or neither"
            raise BertiniError(msg)

        # user did not specify start or final parameters
        if 'parameterhomotopy' in ckeys and self._parameter_homotopy['arg'] > 1:
            if not (startp or finalp):
                msg = "specify start and/or final parameters with the keyword arguments `start_parameters' and/or `final_parameters'"
                raise KeyError(msg)

        self._dirname = tempfile.mkdtemp()
        self._bertini = BERTINI
        self._system = system
        self._config = config
        self._complete = False
        self._inputf = []

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

        from naglib.bertini.fileutils import read_points_file
        from naglib.utils import striplines
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
        from os.path import exists
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

    def rerun(self, config={}):
        if not self._complete:
            return self.run()
        else:
            self._config.update(config)
            return self.run()

    def run(self, rerun_on_fail=False):
        # in case the user has changed any of these
        if not BERTINI:
            raise NoBertiniException()

        self._bertini = BERTINI

        if self._parallel and MPIRUN:
            cmd = MPIRUN
            arg = [cmd, '-np', str(PCOUNT), self._bertini]
        else:
            arg = [self._bertini]

        dirname = self._dirname

        input_file = self._write_files()

        if op.exists(dirname + '/instructions'):
            stdin = dirname + '/instructions'
        else:
            stdin = None

        arg += [input_file]

        os.chdir(dirname)
        if stdin:
            stdin = open(stdin, 'r')
        try:
            output = subprocess.check_output(arg, stdin=stdin, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            msg = naglib.bertini.system.proc_err_output(e.output)
            raise BertiniError(msg)

        if stdin:
            stdin.close()

        self._complete = True
        self._output = output

        if rerun_on_fail:
            try:
                self._inputf = self._recover_input()
                data = self._recover_data()
            except:
                data = self.rerun()
        else:
            self._inputf = self._recover_input()
            data = self._recover_data()

        return data

    @property
    def bertini(self):
        return self._bertini
    @bertini.setter
    def bertini(self, bert):
        self._bertini = bert

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
    def inputf(self):
        return self._inputf

    @property
    def output(self):
        return self._output

    @property
    def tracktype(self):
        return self._config["tracktype"]
    @tracktype.setter
    def tracktype(self, val):
        if  val not in range(-4, 8):
            raise ValueError("tracktype must be an integer between -4 and 7 (inclusive)")

        self._config["tracktype"] = int(val)


def parameter_homotopy(system: list):
    pass
