class Component(object):

    def __init__(self, dim, degree, witness):
        self._dim = dim
        self._degree = degree
        self._witness = witness

    def __str__(self):
        dim = self._dim
        deg = self._degree
        return '{0}-dim component of degree {1}'.format(dim, deg)

    def __repr__(self):
        return self.__str__()

    @property
    def degree(self):
        return self._degree
    @degree.setter
    def degree(self, deg):
        self._degree = deg

class WitnessPoint(object):

    def __init__(self, dim, component_id, pt, projective=True):
        """Initialize the WitnessPoint object.

        Keyword arguments:
        dim -- the dimension of the component to which the point belongs
        component_id -- the component number for this dimension, as assigned by Bertini
        pt -- the coordinates of the witness point (an mpc vector)
        projective -- True if the point is projective, otherwise False
        """
        self._dim = dim
        self._component_id = component_id
        if not hasattr(pt, '__iter__'):
            self._pt = (pt,)
        else:
            self._pt = tuple(pt)
        self._projective = projective

    def __str__(self):
        return '(' + ', '.join(str(p) for p in self._pt) + ')'

    def __repr__(self):
        return 'WitnessPoint({0},{1},{2},projective={3})'.format(self._dim,
                                                                 self._component_id,
                                                                 self._pt,
                                                                 self._projective)

    @property
    def dim(self):
        return self._dim

    @property
    def component_id(self):
        return self._component_id
    @component_id.setter
    def component_id(self, c):
        if type(c) != int:
            raise(TypeError('component_id needs to be an int'))
        self._component_id = c

class WitnessSet(object):

    def __init__(self, f, L, W):
        self._f = f
        self._L = L
        self._W = W