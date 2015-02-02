class Component(object):

    def __init__(self, dim, degree, witness):
        self._dim = dim
        self._degree = degree
        self._witness = witness

    def __str__(self):
        dim = self._dim
        deg = self._degree
        return '{0}-dim component of degree {1}'.format(dim, deg)

    @property
    def degree(self):
        return self._degree
    @degree.setter
    def degree(self, deg):
        self._degree = deg

class WitnessSet(object):

    def __init__(self, f, L, W):
        self._f = f
        self._L = L
        self._W = W