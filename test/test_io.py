from context import op, ZERO_DIM_BASE, POS_DIM_BASE

from bertini.io.io import read_input_file
from bertini.io.input_file import BertiniInput, BertiniConfig


def is_empty(x):
    return len(x) == 0


def test_read_input1():
    config, inputs, misclines = read_input_file(op.join(POS_DIM_BASE, "sampling", "input"))

    assert isinstance(config, BertiniConfig)
    assert config.tracktype == 2
    assert config.mptype == 2
    assert config.precision == 96
    assert config.coeffbound == 1000
    assert config.degreebound == 5
    assert config.ampmaxprec == 1024
    assert config.parameterhomotopy == 0
    assert config.randomseed == 0

    assert isinstance(inputs, BertiniInput)
    assert is_empty(inputs.constant)
    assert is_empty(inputs.variable)
    assert is_empty(inputs.hom_variable_group)
    assert is_empty(inputs.subfunction)
    assert is_empty(inputs.parameter)
    assert is_empty(inputs.random)
    assert is_empty(inputs.random_real)
    assert is_empty(inputs.pathvariable)
    assert inputs.ndims == 3
    assert len(inputs.variable_group) == 1
    assert inputs.variable_group[0] == ["x", "y", "z"]
    assert len(inputs.function) == 3
    assert inputs.function["f1"] == "(y-x^2)*(x^2+y^2+z^2-1)*(x-0.5)"
    assert inputs.function["f2"] == "(z-x^3)*(x^2+y^2+z^2-1)*(y-0.5)"
    assert inputs.function["f3"] == "(y-x^2)*(z-x^3)*(x^2+y^2+z^2-1)*(z-0.5)"

    assert isinstance(misclines, list) and len(misclines) == 0


def test_read_input2():
    config, inputs, misclines = read_input_file(op.join(ZERO_DIM_BASE, "parameter_homotopy", "input"))

    assert isinstance(config, BertiniConfig)
    assert config.tracktype == 0
    assert config.mptype == 2
    assert config.precision == 96
    assert config.coeffbound == 1000
    assert config.degreebound == 5
    assert config.ampmaxprec == 1024
    assert config.parameterhomotopy == 1
    assert config.randomseed == 10191

    assert isinstance(inputs, BertiniInput)
    assert is_empty(inputs.constant)
    assert is_empty(inputs.variable)
    assert is_empty(inputs.hom_variable_group)
    assert is_empty(inputs.subfunction)
    assert len(inputs.parameter) == 7
    assert inputs.parameter["a0"] is None
    assert inputs.parameter["a1"] is None
    assert inputs.parameter["a2"] is None
    assert inputs.parameter["a3"] is None
    assert inputs.parameter["a4"] is None
    assert inputs.parameter["a5"] is None
    assert inputs.parameter["a6"] is None
    assert is_empty(inputs.random)
    assert is_empty(inputs.random_real)
    assert is_empty(inputs.pathvariable)
    assert inputs.ndims == 1
    assert len(inputs.variable_group) == 1
    assert inputs.variable_group[0] == ["x"]
    assert len(inputs.function) == 1
    assert inputs.function["f"] == "a0+x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))"

    assert isinstance(misclines, list) and len(misclines) == 0
