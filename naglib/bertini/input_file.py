from collections import OrderedDict
from numbers import Number

DEFAULT_CONFIGS = OrderedDict(tracktype=0,
                              parameterhomotopy=0)
LEGAL_VALUES = {"tracktype": lambda x: x in range(-4, 8),
                "parameterhomotopy": range(0, 3)}

class BertiniConfig(object):
    TEVALP    = -4
    TEVALPJ   = -3
    TNEWTP    = -2
    TNEWTPJ   = -1
    TZERODIM  =  0 # parallel
    TPOSDIM   =  1 # parallel
    TSAMPLE   =  2
    TMEMTEST  =  3
    TPRINTWS  =  4
    TPROJECT  =  5
    TISOSTAB  =  6
    TREGENEXT =  7 # parallel

    def __init__(self, **kwargs):
        for arg_name in DEFAULT_CONFIGS:
            if arg_name in kwargs:
                arg_val = kwargs.pop(arg_name)
                self.__setattr__(arg_name, arg_val)
            else:
                self.__setattr__(arg_name, DEFAULT_CONFIGS[arg_name])

        if kwargs:
            key, _ = kwargs.popitem()
            raise AttributeError(f"'BertiniConfig' has no attribute '{key}'")

        self._validate()

    def __str__(self):
        s = "CONFIG\n"
        for key, val in DEFAULT_CONFIGS.items():
            instance_val = self.__getattribute__(key)
            if instance_val != val:
                s += f"{key}:{instance_val};\n"
        s += "END;"
        return s

    def _validate(self):
        """Ensure combinations of parameters play nicely."""
        pass

    @property
    def parameterhomotopy(self):
        return self._parameterhomotopy
    @parameterhomotopy.setter
    def parameterhomotopy(self, val):
        is_legal = LEGAL_VALUES["parameterhomotopy"]
        if isinstance(val, str) or (not isinstance(val, int) and isinstance(val, Number)):
            val = int(val)

        if not is_legal(val):
            raise ValueError(f"parameterhomotopy must take one of the following values: {','.join(map(str, legal_vals))}")
        self._parameterhomotopy = val

    @property
    def tracktype(self):
        return self._tracktype
    @tracktype.setter
    def tracktype(self, val):
        is_legal = LEGAL_VALUES["tracktype"]
        if isinstance(val, str) or (not isinstance(val, int) and isinstance(val, Number)):
            val = int(val)

        if not is_legal(val)
            raise ValueError(f"tracktype must take one of the following values: {','.join(map(str, legal_vals))}")
        self._tracktype = int(val)


class BertiniInput(object):
    def __init__(self):
        pass
