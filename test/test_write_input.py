# -*- coding: utf-8 -*-

from context import op, TESTBASE

from collections import deque, OrderedDict

class TestAbInitio:
    def setup(self):
        self.inputs = dict(variable_group=[["x"]],
                           variable=[],
                           hom_variable_group=[],
                           pathvariable=[],
                           random=[],
                           random_real=[],
                           constant={},
                           subfunction={},
                           parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                           function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))
        self.config = {"ParameterHomotopy": 1}
        self.input_file = op.join(TESTBASE, "input")
