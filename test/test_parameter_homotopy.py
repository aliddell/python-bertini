from context import op, ZERO_DIM_BASE

from collections import OrderedDict

import numpy as np

from naglib.bertini.io import read_points_file
from naglib.bertini.run import BertiniRun
from naglib.bertini.input_file import BertiniInput, BertiniConfig


class TestParameterHomotopy():
    def setup_method(self):
        self.working_dir = op.join(ZERO_DIM_BASE, "parameter_homotopy")
        self.inputs = BertiniInput(variable_group=[["x"]],
                                   parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                                   function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))

        self.config = BertiniConfig(parameterhomotopy=1, randomseed=10191)

        brun = BertiniRun(self.config, self.inputs)
        self.ab_initio_result = brun.run(dirname=self.working_dir)

    def test_ab_initio(self):
        finite_solutions = np.array([[0.78101599+0.18185757j,  0.33253178+1.4856318j,
                                      -0.71413828+1.01669647j, -0.53917013-0.04838107j,
                                      -0.65402653-0.7627159j ,  0.43428293-0.98042348j]])
        assert np.linalg.norm(finite_solutions - self.ab_initio_result.finite_solutions) < 1e-8

    # def test_parameter_homotopy2(self):



def test_parameter_homotopy2():
    inputs = BertiniInput(variable_group=[["x"]],
                          parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                          function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))

    start_parameters = read_points_file(op.join(ZERO_DIM_BASE, "parameter_homotopy", "start_parameters"))
    start_points = read_points_file(op.join(ZERO_DIM_BASE, "parameter_homotopy", "start"))
    config = BertiniConfig(parameterhomotopy=2, randomseed=10191)
    brun = BertiniRun(config, inputs)
