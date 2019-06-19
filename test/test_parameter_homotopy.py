from context import op, BASEDIR

from collections import OrderedDict

from naglib.bertini.run import BertiniRun
from naglib.bertini.input_file import BertiniInput, BertiniConfig


def test_ab_initio():
    inputs = BertiniInput(variable_group=[["x"]],
                          parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                          function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))

    config = BertiniConfig(parameterhomotopy=1)
    brun = BertiniRun(config, inputs)
    brun.run(dirname=op.join(BASEDIR, "test", "param-htpy"))

