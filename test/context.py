import os.path as op
import sys

print(__file__)
BASEDIR = op.abspath(op.join(op.dirname(__file__), ".."))
sys.path.insert(0, BASEDIR)
