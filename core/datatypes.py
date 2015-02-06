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

    def __init__(self, dim, component_id, pt, isprojective=True):
        """Initialize the WitnessPoint object.

        Keyword arguments:
        dim -- the dimension of the component to which the point belongs
        component_id -- the component number for this dimension, as assigned by Bertini
        pt -- the coordinates of the witness point (an mpc vector)
        isprojective -- True if the point is projective, otherwise False
        """
        self._dim = dim
        self._component_id = component_id
        if not hasattr(pt, '__iter__'):
            self._pt = (pt,)
        else:
            self._pt = tuple(pt)
        self._isprojective = isprojective

    def __str__(self):
        return '(' + ', '.join(str(p) for p in self._pt) + ')'

    def __repr__(self):
        return 'WitnessPoint({0},{1},{2},projective={3})'.format(self._dim,
                                                                 self._component_id,
                                                                 self._pt,
                                                                 self._isprojective)

    @property
    def dim(self):
        return self._dim

    @property
    def component_id(self):
        return self._component_id
    @property
    def pt(self):
        return self._pt
    @property
    def isprojective(self):
        return self._isprojective
    
    def dehomogenize(self):
        """Dehomogenize the witness point
        """
        dim = self._dim
        component_id = self._component_id
        pt = self._pt
        if not self._isprojective:
            newpt = self
        else:
            first = pt[0]
            pt = pt[1:]
            # first should not be equal to 0
            pt = [pt[i]/first for i in range(len(pt))]
            newpt = WitnessPoint(dim, component_id, pt, False)
            
        return newpt

class WitnessSet(object):

    def __init__(self, f, L, W):
        self._f = f
        self._L = L
        self._W = W