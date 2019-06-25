import os.path as op
import sys

print(__file__)
BASEDIR = op.abspath(op.join(op.dirname(__file__), ".."))
sys.path.insert(0, BASEDIR)

ZERO_DIM_BASE = op.join(BASEDIR, "test", "data", "zero_dim")
POS_DIM_BASE = op.join(BASEDIR, "test", "data", "pos_dim")
