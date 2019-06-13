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

    def __init__(self):
        self._tracktype = self._tracktype_default()
        self._parameterhomotopy = self._parameterhomotopy_default()

    def _parameterhomotopy_default(self):
        return 0
    def _parameterhomotopy_legal(self):
        return int, range(2)
    def _tracktype_default(self):
        return 0
    def _tracktype_legal(self):
        return int, range(-4, 8)

    @property
    def parameterhomotopy(self):
        return self._parameterhomotopy
    @parameterhomotopy.setter
    def parameterhomotopy(self, val):
        legal_type, legal_vals = self._parameterhomotopy_legal()
        if val not in legal_vals:
            raise ValueError(f"parameterhomotopy must be of type {legal_type} and take one of the following values: {','.join(map(str, legal_vals))}")
        self._tracktype = legal_type(val)

    @property
    def tracktype(self):
        return self._tracktype
    @tracktype.setter
    def tracktype(self, val):
        legal_type, legal_vals = self._tracktype_legal()
        if val not in legal_vals:
            raise ValueError(f"tracktype must be of type {legal_type} and take one of the following values: {','.join(map(str, legal_vals))}")
        self._tracktype = legal_type(val)
