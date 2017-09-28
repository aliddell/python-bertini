from sympy import I, Float, Rational, ShapeError, sympify, Matrix

from naglib.exceptions import ExitSpaceError, AffineInfinityException
from naglib.startup import TOL

def scalar_num(x):
    """
    Determine if x is a scalar type for purposes of multiplication
    """
    from sympy import Number
    x = sympify(str(x))
    re, im = x.as_real_imag()

    return isinstance(re, Number) and isinstance(im, Number)

class NAGobject(object):
    """
    A meta class. Nothing here (yet)
    """
    pass

class Point(NAGobject):
    """
    A point in affine or projective space
    """

    def __init__(self, coordinates):
        """
        Initialize the Point object
        """
        try:
            coordinates = list(coordinates)
        except TypeError:
            coordinates = [coordinates]

        coordinates = [sympify(c) for c in coordinates]
        self._coordinates = Matrix(coordinates)

    def __add__(self, other):
        """
        x.__add__(y) <==> x + y
        """
        cls = self.__class__
        if not isinstance(other, cls):
            t = type(other)
            msg = "unsupported operand type(s) for +: '{0}' and '{1}'".format(cls, t)
            raise TypeError(msg)

        sco = self._coordinates
        oco = other._coordinates

        if len(sco) != len(oco):
            msg = "dimension mismatch"
            raise ShapeError(msg)

        return cls(list(sco + oco))

    def __sub__(self, other):
        """
        x.__sub__(y) <==> x - y
        """
        cls = self.__class__
        if not isinstance(other, cls):
            t = type(other)
            msg = "unsupported operand type(s) for -: '{0}' and '{1}'".format(cls, t)
            raise TypeError(msg)

        return self + -other

    def __neg__(self):
        """
        x.__neg__() <==> -x
        """
        cls = self.__class__
        coordinates  = self._coordinates
        return cls(list(-coordinates))

    def __mul__(self, other):
        """
        x.__mul__(y) <==> x*y
        """
        cls = self.__class__
        if not scalar_num(other):
            t = type(other)
            msg = "unsupported operand type(s) for *: '{0}' and '{1}'".format(cls, t)
            raise TypeError(msg)
        else:
            coordinates = self._coordinates
            return cls(other * coordinates)

    def __rmul__(self, other):
        """
        x.__rmul__(y) <==> y*x
        """
        cls = self.__class__
        if not scalar_num(other):
            t = type(other)
            msg = "unsupported operand type(s) for *: '{0}' and '{1}'".format(t, cls)
            raise TypeError(msg)
        else:
            return self * other

    def __div__(self, other):
        """
        x.__div__(y) <==> x/y
        """
        return self.__truediv__(other)

    def __truediv__(self, other):
        """
        x.__truediv__(y) <==> x/y
        """
        cls = self.__class__
        if not scalar_num(other):
            t = type(other)
            msg = "unsupported operand type(s) for /: '{0}' and '{1}'".format(cls, t)
            raise TypeError(msg)
        else:
            if other == 0:
                msg = "division by zero"
                raise ZeroDivisionError(msg)
            coordinates = self._coordinates
            return cls((1.0/other) * coordinates)

    def __eq__(self, other):
        """
        x.__eq__(y) <==> x == y
        """
        cls = self.__class__
        if not isinstance(other, cls):
            return False

        sco = self._coordinates
        oco = other._coordinates

        return (sco == oco)

    def __len__(self):
        """
        x.__len__() <==> len(x)
        """
        coordinates = self._coordinates
        return len(coordinates)

    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        coordinates = self._coordinates
        return coordinates[key]

    def __setitem__(self, key, value):
        """
        x.__setitem__(y, z) <==> x[y] = z
        """
        if not scalar_num(value):
            msg = "must assign a number"
            raise TypeError(msg)

        coordinates = self._coordinates
        coordinates[key] = sympify(value)

    def __setslice__(self, i, j, sequence):
        """
        x.__setslice__(i,j,sequence) <==> x[i:j] = sequence
        """
        coordinates = self._coordinates

        sequence = sympify(list(sequence))
        for s in sequence:
            if not scalar_num(s):
                msg = "must assign a number"
                raise TypeError(msg)

        cocopy = coordinates[:]
        cocopy[i:j] = sequence
        self._coordinates = Matrix(cocopy)

    def cat(self, other):
        """
        concatenate other onto self
        """
        cls = self.__class__
        if not isinstance(other, cls):
            t = type(other)
            msg = "I don't know how to concatenate type {0} to type {1}".format(t, cls)
            raise TypeError(msg)
        scoords = list(self._coordinates)
        ocoords = list(other._coordinates)

        return cls(scoords + ocoords)

    def float(self, prec=None):
        cls = self.__class__
        coordinates = self._coordinates
        newcoords = []
        for c in coordinates:
            real,imag = c.as_real_imag()
            if prec:
                real = Float(real, prec)
                imag = Float(imag, prec)
            else:
                real = Float(real)
                imag = Float(imag)
            newcoords.append(real + I*imag)

        return cls(newcoords)

    def insert(self, index, item):
        if not scalar_num(item):
            msg = "only takes scalar number types"
            raise TypeError(msg)
        coordinates = list(self._coordinates)
        coordinates.insert(index, item)
        self._coordinates = Matrix(coordinates)

    def is_zero(self, tol=TOL*10):
        coordinates = self._coordinates
        summa = sum([abs(c)**2 for c in coordinates])**0.5
        return summa < tol

    def normalized(self):
        """
        Returns the normalized version of ``self''
        """
        coordinates = self._coordinates
        return self.__class__(coordinates.normalized())

    def pop(self, index=-1):
        coordinates = list(self._coordinates)
        popped = coordinates.pop(index)
        self._coordinates = Matrix(coordinates)

        return popped

    def rational(self):
        """
        Returns a rational approximation of self
        """
        cls = self.__class__
        coordinates = self._coordinates
        newcoords = []
        for c in coordinates:
            real,imag = c.as_real_imag()
            real = Rational(real)
            imag = Rational(imag)
            newcoords.append(real + I*imag)

        return cls(newcoords)

    def is_real(self, tol=1e-13):
        coordinates = self._coordinates
        return coordinates.as_real_imag()[1].norm() < tol

    @property
    def coordinates(self):
        return self._coordinates

class AffinePoint(Point):
    """
    Point object living in affine space
    """
    def __init__(self, coordinates):
        super(AffinePoint, self).__init__(coordinates)

    def __repr__(self):
        coordinates = [str(c.n()) for c in self._coordinates]
        repstr = 'AffinePoint(['
        if len(coordinates) > 6:
            repstr += ', '.join(coordinates[:3] + ['...'] + coordinates[-3:]) + '])'
        else:
            repstr = 'AffinePoint({0})'.format(', '.join(coordinates))

        return repstr

    def __str__(self):
        """
        x.__str__ <==> str(x)
        """
        coordinates = self._coordinates
        repstr = '[' + ', '.join([str(c) for c in coordinates]) + ']'

        return repstr

    def __abs__(self):
        """
        x.__abs__ <==> abs(x)
        """
        coordinates = self._coordinates
        return coordinates.norm()

    def homogenize(self):
        """
        """
        homcoordinates = [1] + list(self._coordinates)

        return ProjectivePoint(homcoordinates)

    def norm(self, ord=None):
        """
        Return the norm of AffinePoint
        Acceptable values:
        None <==> 2-norm
         inf <==> max(abs(x))
        -inf <==> min(abs(x))
          n  <==> sum(abs(x)**n)**(1./n)
        """
        coordinates = self._coordinates
        return coordinates.norm(ord)

    @property
    def dim(self):
        return len(self._coordinates)

class ProjectivePoint(Point):
    """
    Point object living in affine space
    """
    def __init__(self, coordinates, tol=TOL*10):
        super(ProjectivePoint, self).__init__(coordinates)

        if self.is_zero(tol):
            corstr = '[' + ' : '.join([str(c) for c in self._coordinates]) + ']'
            msg = "{0} is not in projective space".format(corstr)
            raise ExitSpaceError(msg)
        self._dim = len(self._coordinates) - 1
        self._tol = tol

    def __repr__(self):
        coordinates = [str(c.n()) for c in self._coordinates]
        repstr = 'ProjectivePoint(['
        if len(coordinates) > 6:
            repstr += ', '.join(coordinates[:3] + ['...'] + coordinates[-3:]) + '])'
        else:
            repstr = 'ProjectivePoint({0})'.format(', '.join(coordinates))

        return repstr

    def __str__(self):
        """
        x.__str__ <==> str(x)
        """
        coordinates = self._coordinates
        repstr = '[' + ' : '.join([str(c) for c in coordinates]) + ']'

        return repstr

    def __add__(self, other):
        """
        x.__add__(y) <==> x + y
        """
        res = super(ProjectivePoint, self).__add__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

        return res

    def __sub__(self, other):
        """
        x.__sub__(y) <==> x - y
        """
        res = super(ProjectivePoint, self).__sub__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

        return res

    def __mul__(self, other):
        """
        x.__mul__(y) <==> x*y
        """
        res = super(ProjectivePoint, self).__mul__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

        return res

    def __rmul__(self, other):
        """
        x.__rmul__(y) <==> y*x
        """
        res = super(ProjectivePoint, self).__rmul__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

    def __div__(self, other):
        """
        x.__div__(y) <==> x/y
        """
        res = super(ProjectivePoint, self).__div__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

        return res

    def __truediv__(self, other):
        """
        x.__truediv__(y) <==> x/y
        """
        res = super(ProjectivePoint, self).__truediv__(other)
        if res.is_zero(self._tol):
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

        return res

    def __eq__(self, other):
        """
        x.__eq__(u=y) <==> x == y
        """
        if not isinstance(other, ProjectivePoint):
            return False

        # (x0, ..., xn) ==  (y0, ..., yn) if there exists c s.t.
        # (x0, ..., xn) == c(y0, ..., yn)
        sho = self.canonical()
        oho = other.canonical()
        sco = sho._coordinates
        oco = oho._coordinates

        return sco == oco

    def __setitem__(self, key, value):
        """
        x.__setitem__(y, z) <==> x[y] = z
        """
        coordinates = self._coordinates.copy()
        super(ProjectivePoint, self).__setitem__(key, value)
        if self.is_zero(self._tol):
            self._coordinates = coordinates
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

    def __setslice__(self, i, j, sequence):
        """
        x.__setitem__(y, z) <==> x[y] = z
        """
        coordinates = self._coordinates.copy()
        super(ProjectivePoint, self).__setslice__(i, j, sequence)
        if self.is_zero(self._tol):
            self._coordinates = coordinates
            msg = 'you have left projective space'
            raise ExitSpaceError(msg)

    def at_infinity(self, tol=TOL):
        coordinates = self._coordinates
        c0 = coordinates[0]
        div_coord = [c0/c for c in coordinates[1:]]
        return sum([abs(c)**2 for c in div_coord])**(0.5) < tol

    def canonical(self):
        """
        Return ProjectivePoint with coordinates rescaled
        """
        coordinates = self._coordinates
        nz = 0
        while abs(coordinates[nz]) < TOL and nz < len(coordinates):
            nz += 1
        c1 = coordinates[nz]

        return ProjectivePoint(coordinates/c1)

    def dehomogenize(self):
        """
        Return Affine Point, if not at at infinity
        """
        if not self.at_infinity():
            coordinates = self._coordinates
            return AffinePoint([c/coordinates[0] for c in coordinates[1:]])
        else:
            msg = "cannot create affine point from {0}; point is at infinity".format(self)
            raise AffineInfinityException(msg)

    @property
    def dim(self):
        return len(self._coordinates)-1
    @property
    def tol(self):
        return self._tol
    @tol.setter
    def tol(self, t):
        self._tol = t
