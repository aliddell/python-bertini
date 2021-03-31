from copy import deepcopy
from collections import OrderedDict
from numbers import Number


def _list_of_str(obj):
    return isinstance(obj, list) and all(isinstance(x, str) for x in obj)


def _list_of_lists_of_str(obj):
    return isinstance(obj, list) and all(_list_of_str(l) for l in obj)


def _dict_of_type(obj, dtype):
    return isinstance(obj, dict) and _list_of_str(list(obj.keys())) \
           and all((isinstance(v, dtype) or v is None) for v in obj.values())


def _ordereddict_of_type(obj, dtype):
    return isinstance(obj, OrderedDict) and _list_of_str(list(obj.keys())) \
           and all((isinstance(v, dtype) or v is None) for v in obj.values())


PARAMETERS = OrderedDict(tracktype={"default": 0,
                                    "is valid": lambda x: x in range(-4, 8)},
                         tracktolbeforeeg={"default": 1e-11,
                                           "is valid": lambda x: isinstance(x, Number) and x > 0},
                         mptype={"default": 2,
                                 "is valid": lambda x: x in range(0, 3)},
                         precision={"default": 96,
                                    "is valid": lambda x: isinstance(x, int) and x >= 64},
                         coeffbound={"default": 1000.,  # default for user-defined homotopies only
                                     "is valid": lambda x: isinstance(x, float) and x > 0},
                         degreebound={"default": 5.,  # default for user-defined homotopies only
                                      "is valid": lambda x: isinstance(x, float) and x > 0},
                         ampmaxprec={"default": 1024,
                                     "is valid": lambda x: isinstance(x, int) and x >= 64},
                         finaltol={"default": 1.,
                                   "is valid": lambda x: isinstance(x, Number) and x > 0},
                         imagthreshold={"default": 1.,
                                        "is valid": lambda x: isinstance(x, Number) and x > 0},
                         parameterhomotopy={"default": 0,
                                            "is valid": lambda x: x in range(0, 3)},
                         randomseed={"default": 0,
                                     "is valid": lambda x: isinstance(x, int) and x >= 0})


def _validate_param(name, val):
    if name not in PARAMETERS:
        return val, False

    default_val = PARAMETERS[name]["default"]
    is_valid = PARAMETERS[name]["is valid"]
    if not isinstance(val, type(default_val)):
        try:
            val = type(default_val)(val)  # cast val to the same type as default_val
        except ValueError:  # failed to cast
            return val, False

    return val, is_valid(val)


INPUT_TYPES = OrderedDict(variable_group={"default": [],
                                          "is valid": _list_of_lists_of_str},
                          variable={"default": [],
                                    "is valid": _list_of_lists_of_str},
                          hom_variable_group={"default": [],
                                              "is valid": _list_of_lists_of_str},
                          pathvariable={"default": [],
                                        "is valid": _list_of_str},
                          random={"default": [],
                                  "is valid": _list_of_str},
                          random_real={"default": [],
                                       "is valid": _list_of_str},
                          constant={"default": {},
                                    "is valid": lambda x: _dict_of_type(x, Number)},
                          subfunction={"default": {},
                                       "is valid": lambda x: _dict_of_type(x, str)},
                          parameter={"default": OrderedDict(),
                                     "is valid": lambda x: _ordereddict_of_type(x, str)},
                          function={"default": OrderedDict(),
                                    "is valid": lambda x: _ordereddict_of_type(x, str)})


class BertiniInput(object):
    def __init__(self, **kwargs):
        for arg_name in INPUT_TYPES:
            if arg_name in kwargs:
                arg_val = kwargs.pop(arg_name)
                self.__setattr__(arg_name, arg_val)
            else:
                self.__setattr__(arg_name, deepcopy(INPUT_TYPES[arg_name]["default"]))

        if kwargs:
            key, _ = kwargs.popitem()
            raise AttributeError(f"BertiniInput has no attribute '{key}'")

        self._validate()

    def __str__(self):
        s = "INPUT\n"
        # name declarations
        if self.variable_group:
            for variable_group in self.variable_group:
                s += f"variable_group {','.join(variable_group)};\n"

        if self.variable:
            for variable in self.variable:
                s += f"variable {','.join(variable)};\n"

        if self.hom_variable_group:
            for hom_variable_group in self.hom_variable_group:
                s += f"hom_variable_group {','.join(hom_variable_group)};\n"

        if self.pathvariable:
            for pathvariable in self.pathvariable:
                s += f"pathvariable {','.join(pathvariable)};\n"

        if self.random:
            for random in self.random:
                s += f"random {','.join(random)};\n"

        if self.random_real:
            for random_real in self.random_real:
                s += f"random_real {','.join(random_real)};\n"

        if self.parameter:
            parameters = ",".join(self.parameter.keys())
            s += f"parameter {parameters};\n"

        if self.constant:
            constant = ",".join(self.constant.keys())
            s += f"constant {constant};\n"

        function = ",".join(self.function.keys())
        s += f"function {function};\n"

        # optional assignments (parameters, constants, &c.)
        if self.parameter:
            for key, val in self.parameter.items():
                if val is None:
                    continue
                s += f"{key} = {val};\n"

        if self.constant:
            for key, val in self.constant.items():
                s += f"{key} = {val};\n"

        if self.subfunction:
            for key, val in self.subfunction.items():
                s += f"{key} = {val};\n"

        # function definitions
        for key, val in self.function.items():
            s += f"{key} = {val};\n"

        s += "END;"

        return s

    def _validate(self):
        if self.variable and (self.variable_group or self.hom_variable_group):
            raise ValueError("specify EITHER variable OR some combination of variable_group and hom_variable_group")

    @property
    def variable_group(self):
        return self._variable_group

    @variable_group.setter
    def variable_group(self, val):
        is_valid = INPUT_TYPES["variable_group"]["is valid"]
        if not is_valid(val):
            raise ValueError("variable_group must be a list of lists of str")
        self._variable_group = val

    @property
    def variable(self):
        return self._variable

    @variable.setter
    def variable(self, val):
        is_valid = INPUT_TYPES["variable"]["is valid"]
        if not is_valid(val):
            raise ValueError("variable must be a list of lists of str")
        self._variable = val

    @property
    def hom_variable_group(self):
        return self._hom_variable_group

    @hom_variable_group.setter
    def hom_variable_group(self, val):
        is_valid = INPUT_TYPES["hom_variable_group"]["is valid"]
        if not is_valid(val):
            raise ValueError("hom_variable_group must be a list of lists of str")
        self._hom_variable_group = val

    @property
    def pathvariable(self):
        return self._pathvariable

    @pathvariable.setter
    def pathvariable(self, val):
        is_valid = INPUT_TYPES["pathvariable"]["is valid"]
        if not is_valid(val):
            raise ValueError("pathvariable must be a list of str")
        self._pathvariable = val

    @property
    def random(self):
        return self._random

    @random.setter
    def random(self, val):
        is_valid = INPUT_TYPES["random"]["is valid"]
        if not is_valid(val):
            raise ValueError("random must be a list of str")
        self._random = val

    @property
    def random_real(self):
        return self._random_real

    @random_real.setter
    def random_real(self, val):
        is_valid = INPUT_TYPES["random_real"]["is valid"]
        if not is_valid(val):
            raise ValueError("random_real must be a list of str")
        self._random_real = val

    @property
    def constant(self):
        return self._constant

    @constant.setter
    def constant(self, val):
        is_valid = INPUT_TYPES["constant"]["is valid"]
        if not is_valid(val):
            raise ValueError("constant must be a dict of numeric")
        self._constant = val

    @property
    def subfunction(self):
        return self._subfunction

    @subfunction.setter
    def subfunction(self, val):
        is_valid = INPUT_TYPES["subfunction"]["is valid"]
        if not is_valid(val):
            raise ValueError("subfunction must be a dict of str")
        self._subfunction = val

    @property
    def parameter(self):
        return self._parameter

    @parameter.setter
    def parameter(self, val):
        is_valid = INPUT_TYPES["parameter"]["is valid"]
        if not is_valid(val):
            raise ValueError("parameter must be an OrderedDict of str")
        self._parameter = val

    @property
    def function(self):
        return self._function

    @function.setter
    def function(self, val):
        is_valid = INPUT_TYPES["function"]["is valid"]
        if not is_valid(val):
            raise ValueError("function must be an OrderedDict of str")
        self._function = val

    @property
    def ndims(self):
        if self.variable:
            nd = len(sum(self.variable, []))
        else:
            nd = len(sum(self.variable_group, []) + sum(self.hom_variable_group, []))

        return nd
