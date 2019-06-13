# -*- coding: utf-8 -*-

from context import op, TESTBASE

from naglib.bertini.io import read_witness_data

#class TestReadWitnessData:
#    def setup(self):
#

wd = read_witness_data(op.join(TESTBASE, "pos_dim", "basic_pos_dim", "witness_data"))