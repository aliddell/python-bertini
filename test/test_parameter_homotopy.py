# -*- coding: utf-8 -*-

from context import op, TESTBASE

from naglib.bertini.io import read_witness_data

#class TestReadWitnessData:
#    def setup(self):
#

class TestAbInitio:
    def setup(self):
        self.inputs = dict(variable_group=deque(["x"]),
                           variable=deque(),
                           hom_variable_group=deque(),
                           pathvariable=deque(),
                           random=deque(),
                           random_real=deque(),
                           constant=OrderedDict(),
                           function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"),
                           parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                           subfunction=OrderedDict())
        self.config = {"ParameterHomotopy": 1}
        self.input_file = op.join(TESTBASE, "input")
