from context import op, ZERO_DIM_BASE

from collections import OrderedDict

import numpy as np

from bertini.io.io import write_input_file
from bertini.run import BertiniRun
from bertini.io.input_file import BertiniInput, BertiniConfig


class TestParameterHomotopy:
    def setup_method(self):
        self.working_dir = op.join(ZERO_DIM_BASE, "parameter_homotopy")
        self.inputs = BertiniInput(variable_group=[["x"]],
                                   parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                                   function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))

        self.config = BertiniConfig(parameterhomotopy=1, randomseed=10191)

        print("initial run")
        brun = BertiniRun(self.config, self.inputs, mpi_path=None)
        self.ab_initio_result = brun.run(dirname=self.working_dir)

    def test_ab_initio_result(self):
        finite_solutions = np.array([[0.78101599+0.18185757j,  0.33253178+1.4856318j,
                                      -0.71413828+1.01669647j, -0.53917013-0.04838107j,
                                      -0.65402653-0.7627159j ,  0.43428293-0.98042348j]])
        assert np.linalg.norm(finite_solutions - self.ab_initio_result.finite_solutions) < 1e-8

    def test_parameter_homotopy2(self):
        start_parameters = self.ab_initio_result.start_parameters
        final_parameters = np.array([1.1, 2.4, 0.8, 3.6, -0.52, -1.8, 4.4], dtype=np.complex)
        self.config.parameterhomotopy = 2

        brun = BertiniRun(self.config, self.inputs, start=self.ab_initio_result.nonsingular_solutions,
                          start_parameters=start_parameters, final_parameters=final_parameters, mpi_path=None)

        self.final_result = brun.run(dirname=self.working_dir)
        assert isinstance(self.final_result.config, BertiniConfig)
        assert isinstance(self.final_result.inputs, BertiniInput)

    def teardown_method(self):
        self.config.parameterhomotopy = 1
        write_input_file(self.config, self.inputs, op.join(self.working_dir, "input"))