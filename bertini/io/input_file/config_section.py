from copy import deepcopy

from bertini.io.input_file.input_section import PARAMETERS, _validate_param


class BertiniConfig(object):
    TEVALP = -4
    TEVALPJ = -3
    TNEWTP = -2
    TNEWTPJ = -1
    TZERODIM = 0  # parallel
    TPOSDIM = 1  # parallel
    TSAMPLE = 2
    TMEMTEST = 3
    TPRINTWS = 4
    TPROJECT = 5
    TISOSTAB = 6
    TREGENEXT = 7  # parallel

    def __init__(self, **kwargs):
        for arg_name in PARAMETERS:
            if arg_name in kwargs:
                arg_val = kwargs.pop(arg_name)
                self.__setattr__(arg_name, arg_val)
            else:
                self.__setattr__(arg_name, deepcopy(PARAMETERS[arg_name]["default"]))

        if kwargs:
            key, _ = kwargs.popitem()
            raise AttributeError(f"BertiniConfig has no attribute '{key}'")

        self._validate()

    def __str__(self):
        s = "CONFIG\n"
        for key, val in PARAMETERS.items():
            instance_val = self.__getattribute__(key)
            if instance_val != val["default"]:
                s += f"{key}:{instance_val};\n"
        s += "END;"
        return s

    def _validate(self):
        """Ensure combinations of parameters play nicely."""
        if self.mptype != 1 and self.precision != PARAMETERS["precision"]["default"]:
            raise ValueError("you have set a non-default precision but have specified "
                             f"{'double' if self.mptype == 0 else 'adaptive'} precision")

    def needs_component(self):
        return self.tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTWS,
                                  self.TPROJECT, self.TREGENEXT)

    def needs_start_points(self):
        nsp = self.tracktype in (self.TEVALP, self.TEVALPJ, self.TNEWTP, self.TNEWTPJ, self.TMEMTEST,
                                 self.TISOSTAB) or self.parameterhomotopy == 2
        return nsp

    def needs_sample_count(self):
        return self.tracktype == self.TSAMPLE

    def needs_projection_variables(self):
        return self.tracktype == self.TPROJECT

    @property
    def ampmaxprec(self):
        return self._ampmaxprec

    @ampmaxprec.setter
    def ampmaxprec(self, val):
        val, is_valid = _validate_param("ampmaxprec", val)
        if not is_valid:
            raise ValueError("ampmaxprec must be an integer greater than or equal to 64")
        self._ampmaxprec = val

    @property
    def coeffbound(self):
        return self._coeffbound

    @coeffbound.setter
    def coeffbound(self, val):
        val, is_valid = _validate_param("coeffbound", val)
        if not is_valid:
            raise ValueError("coeffbound must be a positive double")
        self._coeffbound = val

    @property
    def degreebound(self):
        return self._degreebound

    @degreebound.setter
    def degreebound(self, val):
        val, is_valid = _validate_param("degreebound", val)
        if not is_valid:
            raise ValueError("degreebound must be a positive double")
        self._degreebound = val

    @property
    def finaltol(self):
        return self._finaltol

    @finaltol.setter
    def finaltol(self, val):
        val, is_valid = _validate_param("finaltol", val)
        if not is_valid:
            raise ValueError("finaltol must be a positive double")
        self._finaltol = val

    @property
    def imagthreshold(self):
        return self._finaltol

    @imagthreshold.setter
    def imagthreshold(self, val):
        val, is_valid = _validate_param("imagthreshold", val)
        if not is_valid:
            raise ValueError("imagthreshold must be a positive double")
        self._imagthreshold = val

    @property
    def mptype(self):
        return self._mptype

    @mptype.setter
    def mptype(self, val):
        val, is_valid = _validate_param("mptype", val)
        if not is_valid:
            raise ValueError(
                f"mptype must take one of the following values: {','.join(map(str, range(0, 3)))}")
        self._mptype = val

    @property
    def parameterhomotopy(self):
        return self._parameterhomotopy

    @parameterhomotopy.setter
    def parameterhomotopy(self, val):
        val, is_valid = _validate_param("parameterhomotopy", val)
        if not is_valid:
            raise ValueError(
                f"parameterhomotopy must take one of the following values: {','.join(map(str, range(0, 3)))}")
        self._parameterhomotopy = val

    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, val):
        val, is_valid = _validate_param("precision", val)
        if not is_valid:
            raise ValueError("precision must be an integer greater than or equal to 64")
        self._precision = val

    @property
    def randomseed(self):
        return self._randomseed

    @randomseed.setter
    def randomseed(self, val):
        val, is_valid = _validate_param("randomseed", val)
        if not is_valid:
            raise ValueError("randomseed must be a nonnegative integer")
        self._randomseed = val

    @property
    def tracktolbeforeeg(self):
        return self._tracktolbeforeeg

    @tracktolbeforeeg.setter
    def tracktolbeforeeg(self, val):
        val, is_valid = _validate_param("tracktolbeforeeg", val)
        if not is_valid:
            raise ValueError(f"tracktolbeforeeg must be a positive double")
        self._tracktolbeforeeg = val

    @property
    def tracktype(self):
        return self._tracktype

    @tracktype.setter
    def tracktype(self, val):
        val, is_valid = _validate_param("tracktype", val)
        if not is_valid:
            raise ValueError(f"tracktype must take one of the following values: {','.join(map(str, range(-4, 7)))}")
        self._tracktype = int(val)
