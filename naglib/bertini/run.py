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
from naglib.bertini.io import (read_input_file, read_main_data_file, read_points_file, read_witness_data_file,
                               write_input_file, write_points_file, extract_error_message)
from naglib.exceptions import BertiniError


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


class BertiniRun:
    def __init__(self, config, inputs, **kwargs):
        """A single Bertini run.

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
                raise ValueError("specify the variables onto which you wish to project with the keyword argument "
                                 "'projection'")
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
            if kwargs["mpi_path"] is not None and not op.isfile(kwargs["mpi_path"]):
                raise OSError(f"didn't find MPI executable at '{kwargs['mpi_path']}'")
            self._mpi = kwargs["mpi_path"]
        else:
            self._mpi = _which("mpirun")

        self._complete = False

    def _recover_data(self):
        """Recover data from a run.

        Returns
        -------
        data : BertiniResult
        """
        if not self._complete:
            return

        result = read_main_data_file(op.join(self.dirname, "main_data"))

        if self.tracktype == self.config.TZERODIM:
            if op.isfile(op.join(self.dirname), "finite_solutions"):
                result.finite_solutions = read_points_file(op.join(self.dirname, "finite_solutions"),
                                                           multi=self.config.mptype != 0)
            else:
                result.finite_solutions = np.array(0, dtype=np.complex)

            if op.isfile(op.join(self.dirname), "real_finite_solutions"):
                result.real_finite_solutions = read_points_file(op.join(self.dirname, "real_finite_solutions"),
                                                                multi=self.config.mptype != 0)
            else:
                result.real_finite_solutions = np.array(0, dtype=np.complex)

            result.nonsingular_solutions = read_points_file(op.join(self.dirname, "nonsingular_solutions"),
                                                            multi=self.config.mptype != 0)
            result.singular_solutions = read_points_file(op.join(self.dirname, "singular_solutions"),
                                                         multi=self.config.mptype != 0)

            if self.config.parameterhomotopy == 1:
                result.start_parameters = read_points_file(op.join(self.dirname, "start_parameters"), multi=True)
                result.start = read_points_file(op.join(self.dirname, "start"), multi=self.config.mptype != 0)

        return result

    def _recover_input(self):
        """Read main_data and recover input needed to reproduce this run.
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

    # def _write_files(self):
    #     from naglib.bertini.fileutils import fprint
    #
    #     tracktype = self._tracktype
    #     dirname = self._dirname
    #     system = self._system
    #     if not op.exists(dirname):
    #         os.mkdir(dirname)
    #
    #     ### write the system
    #     sysconfig = self._config.copy()
    #     sysconfig.update({'TrackType':tracktype})
    #     inputf = self._write_system(system, config=sysconfig)
    #
    #     ### write out `start', `start_parameters', `final_parameters'
    #     if '_start' in dir(self):
    #         start = self._start
    #         if self._tracktype == self.TMEMTEST:
    #             startfile = dirname + '/member_points'
    #         else:
    #             startfile = dirname + '/start'
    #         fprint(start, startfile)
    #     if self._parameter_homotopy:
    #         phtpy = self._parameter_homotopy
    #         pkeys = phtpy.keys()
    #         if 'start parameters' in pkeys:
    #             startp = phtpy['start parameters']
    #             startpfile = dirname + '/start_parameters'
    #             fprint(startp, startpfile)
    #         if 'final parameters' in pkeys:
    #             finalp = phtpy['final parameters']
    #             finalpfile = dirname + '/final_parameters'
    #             fprint(finalp, finalpfile)
    #
    #     ### write out component information
    #     if '_component' in dir(self):
    #         component = self._component
    #         #cid = component.component_id
    #         dim = component.dim
    #         if tracktype == self.TREGENEXT:
    #             self._write_system(component.system, 'iold', {'TrackType':1})
    #             witness_data = component._construct_witness_data()
    #             self._write_witness_data(witness_data, dirname, filename='wdold')
    #             instructions = ['1', 'iold', 'wdold', str(dim), '0']
    #             self._write_instructions(instructions)
    #         elif tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTWS, self.TPROJECT):
    #             witness_data = component._construct_witness_data()
    #             self._write_witness_data(witness_data, dirname)
    #             if tracktype == self.TSAMPLE:
    #                 sample = self._sample
    #                 instructions = [str(dim), '0', str(sample), '0', 'sampled']
    #                 self._write_instructions(instructions)
    #             elif tracktype == self.TPRINTWS:
    #                 instructions = [str(dim), '0', 'points.out', 'sys.out']
    #                 self._write_instructions(instructions)
    #             elif tracktype == self.TPROJECT:
    #                 instructions = [str(dim), '0']
    #                 self._write_instructions(instructions)
    #
    #                 ### write out projection information
    #                 projection = self._projection
    #                 projnum = ['1' if x in projection else '0' for x in system.variables]
    #                 projnum = [' '.join(projnum)]
    #                 self._write_instructions(projnum, 'projection')
    #
    #     return inputf

    def _write_instructions(self, lines, filename='instructions'):
        filename = self._dirname + '/' + filename
        fh = open(filename, 'w')
        for line in lines:
            fh.write(line + '\n')
        fh.close()

    def run(self, dirname: str = None, tee: bool = True):
        """Run Bertini and collect results.

        Parameters
        ----------
        dirname : str, optional
            Path to working directory.
        tee : bool, optional
            Print to stdout as you go if true.
        Returns
        -------

        """
        # in case the user has changed any of these
        if self._parallel and self._mpi is not None:
            cmd = self._mpi
            arg = [cmd, '-np', str(multiprocessing.cpu_count()), self._bertini, "input"]
        else:
            arg = [self._bertini, "input"]

        self.setup(dirname)  # write files

        if op.isfile(op.join(self.dirname, "/instructions")):
            stdin = op.join(self.dirname, "/instructions")
        else:
            stdin = None

        os.chdir(self.dirname)
        if stdin is not None:
            stdin = open(stdin, "r")

        try:
            proc = subprocess.Popen(arg, stdin=stdin, stdout=subprocess.PIPE,
                                    universal_newlines=True)
        except subprocess.CalledProcessError as e:
            msg = extract_error_message(e.output)
            raise BertiniError(msg)

        output = []
        while True:
            line = proc.stdout.readline()
            if line == "" and proc.poll() is not None:
                break

            line = line.strip()
            output.append(line)
            if tee:
                print(line)

        if stdin is not None:
            stdin.close()

        self._complete = True
        self._output = output

        return self._recover_data()

    def setup(self, dirname: str = None):
        """Write input and instructions files to working directory.

        Parameters
        ----------
        dirname : str
            Path to directory to run Bertini.
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
            write_points_file(self.start_parameters.reshape(1, self.start_parameters.size), op.join(self._dirname, "start_parameters"))
            write_points_file(self.final_parameters.reshape(1, self.final_parameters.size), op.join(self._dirname, "final_parameters"))

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

        self._start = val.astype(np.complex)

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

        self._start_parameters = val.astype(np.complex)

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

        self._final_parameters = val.astype(np.complex)

    @property
    def tracktype(self):
        return self.config.tracktype


def parameter_homotopy(system: list):
    pass
