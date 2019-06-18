from context import op, TESTBASE

from collections import OrderedDict

from naglib.bertini.input_file import BertiniInput, BertiniConfig


def test_inputs_parameter_homotopy():
    inputs = BertiniInput(variable_group=[["x"]],
                          parameter=OrderedDict((f"a{i}", None) for i in range(7)),
                          function=OrderedDict(f="a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"))
    assert inputs.variable_group == [["x"]]
    assert "a0" in inputs.parameter and inputs.parameter["a0"] is None
    assert "a1" in inputs.parameter and inputs.parameter["a1"] is None
    assert "a2" in inputs.parameter and inputs.parameter["a2"] is None
    assert "a3" in inputs.parameter and inputs.parameter["a3"] is None
    assert "a4" in inputs.parameter and inputs.parameter["a4"] is None
    assert "a5" in inputs.parameter and inputs.parameter["a5"] is None
    assert "a6" in inputs.parameter and inputs.parameter["a6"] is None
    assert inputs.function["f"] == "a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"
