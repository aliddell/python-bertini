from sympy import I, Float, Rational, Matrix

from naglib.exceptions import WitnessDataException
from naglib.core.base import NAGobject, Point, AffinePoint, ProjectivePoint

class WitnessPoint(Point):
    """
    A single witness point for an algebraic set
    """
    def __init__(self, coordinates, component_id, is_projective=False, **kwargs):
        """
        Initialize the WitnessPoint object.

        Keyword arguments:
        point        -- AffinePoint or ProjectivePoint
        component_id -- int, the id of the component to which self belongs
        """
        super(WitnessPoint, self).__init__(coordinates)
        self._component_id = component_id
        self._is_projective = is_projective
        if isinstance(coordinates, AffinePoint):
            self._is_projective = False
        elif isinstance(coordinates, ProjectivePoint):
            self._is_projective = True
            
        kkeys = kwargs.keys()
        
        if 'corank' in kkeys:
            self._corank = kwargs['corank']
        else:
            self._corank = -1
        
        if 'condition_number' in kkeys:
            self._condition_number = kwargs['condition_number']
        else:
            self._condition_number = 0
        
        if 'smallest_nonzero' in kkeys:
            self._smallest_nonzero = kwargs['smallest_nonzero']
        else:
            self._smallest_nonzero = 0
            
        if 'largest_zero' in kkeys:
            self._largest_zero = kwargs['largest_zero']
        else:
            self._largest_zero = 0
            
        if 'point_type' in kkeys:
            self._point_type = kwargs['point_type']
        else:
            self._point_type = 0
            
        if 'multiplicity' in kkeys:
            self._multiplicity = kwargs['multiplicity']
        else:
            self._multiplicity = 1
            
        if 'deflations' in kkeys:
            self._deflations = kwargs['deflations']
        else:
            self._deflations = 0
            
        if 'precision' in kkeys:
            self._precision = kwargs['precision']
        else:
            self._precision = 0
            
        if 'last_approximation' in kkeys:
            self._last_approximation = kwargs['last_approximation']
        else:
            self._last_approximation = Matrix([])
        
        # TODO: sanity check projective point

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        coordinates = [str(c.n()) for c in self._coordinates]
        component_id = self._component_id
        is_projective = self._is_projective
        
        if is_projective:
            repstr = 'ProjectivePoint(['
            if len(coordinates) > 6:
                repstr += ', '.join(coordinates[:3] + ['...'] + coordinates[-3:]) + '])'
            else:
                repstr = 'ProjectivePoint({0})'.format(', '.join(coordinates))
        else:
            repstr = 'AffinePoint(['
            if len(coordinates) > 6:
                repstr += ', '.join(coordinates[:3] + ['...'] + coordinates[-3:]) + '])'
            else:
                repstr = 'AffinePoint({0})'.format(', '.join(coordinates))
        
        repstr = 'WitnessPoint({0},{1})'.format(repstr, component_id)
        
        return repstr
    
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        coordinates = self._coordinates
        is_projective = self._is_projective
        
        if is_projective:
            repstr = '[' + ' : '.join([str(c) for c in coordinates]) + ']'
        else:
            repstr = '[' + ', '.join([str(c) for c in coordinates]) + ']'
            
        return repstr
        
    def as_point(self):
        coordinates = self._coordinates
        is_projective = self._is_projective
        if is_projective:
            return ProjectivePoint(coordinates)
        else:
            return AffinePoint(coordinates)
    
    def dehomogenize(self):
        """
        Dehomogenize the witness point
        
        Returns a new WitnessPoint object
        """
        cls = self.__class__
        component_id = self._component_id
        coordinates = self._coordinates
        is_projective = self._is_projective
        
        if is_projective:
            point = ProjectivePoint(coordinates)
        else:
            point = AffinePoint(coordinates)
            
        deh = cls(point, component_id)
        deh._is_projective = self._is_projective
        deh._condition_number = self._condition_number
        deh._corank = self._corank
        deh._deflations = self._deflations
        deh._last_approximation = self._last_approximation
        deh._largest_zero = self._largest_zero
        deh._multiplicity = self._multiplicity
        deh._point_type = self._point_type
        deh._precision = self._precision
        deh._smallest_nonzero = self._smallest_nonzero
        
        return deh
    
    def float(self, prec=None):
        cls = self.__class__
        component_id = self._component_id
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
        flo = cls(coordinates, component_id)
        
        flo._is_projective = self._is_projective
        flo._condition_number = self._condition_number
        flo._corank = self._corank
        flo._deflations = self._deflations
        flo._last_approximation = self._last_approximation
        flo._largest_zero = self._largest_zero
        flo._multiplicity = self._multiplicity
        flo._point_type = self._point_type
        if prec and prec <= self._precision:
            flo._precision = prec
        else:
            flo._precision = self._precision
        flo._smallest_nonzero = self._smallest_nonzero
        
        return flo
    
    def rational(self):
        """
        Returns a rational approximation of self
        """
        cls = self.__class__
        component_id = self._component_id
        coordinates = self._coordinates
        newcoords = []
        for c in coordinates:
            real,imag = c.as_real_imag()
            real = Rational(real)
            imag = Rational(imag)
            newcoords.append(real + I*imag)
        rat = cls(newcoords, component_id)
        
        rat._is_projective = self._is_projective
        rat._condition_number = self._condition_number
        rat._corank = self._corank
        rat._deflations = self._deflations
        rat._last_approximation = self._last_approximation
        rat._largest_zero = self._largest_zero
        rat._multiplicity = self._multiplicity
        rat._point_type = self._point_type
        rat._precision = self._precision
        rat._smallest_nonzero = self._smallest_nonzero
            
        return rat

    @property
    def condition_number(self):
        """
        The condition number of the Jacobian at this point
        """
        return self._condition_number
        
    @property
    def corank(self):
        """
        The corank of the Jacobian at this point
        """
        return self._corank
        
    @property
    def component_id(self):
        return self._component_id
        
    @property
    def deflations(self):
        """
        The number of deflations this point required
        """
        return self._deflations
        
    @property
    def last_approximation(self):
        coordinates = self._last_approximation
        is_projective = self._is_projective
        if is_projective:
            return ProjectivePoint(coordinates)
        else:
            return AffinePoint(coordinates)
        
    @property
    def largest_zero(self):
        """
        The largest `zero' singular value of the Jacobian at this point
        """
        return self._largest_zero
        
    @property
    def multiplicity(self):
        """
        The multiplicity of this point with respect to the system
        """
        return self._multiplicity
    
    @property
    def point_type(self):
        return self._point_type
        
    @property
    def precision(self):
        """
        The numerical precision required by this point
        """
        return self._precision
    
    @property
    def smallest_nonzero(self):
        """
        The smallest nonzero singular value of the Jacobian at this point
        """
        return self._smallest_nonzero

class WitnessSet(NAGobject):
    """
    A witness set for a component
    """
    def __init__(self, system, slice, witness_points, witness_data):
        """Initialize the WitnessSet
        
        Keyword arguments:
        system -- PolynomialSystem, system on which all the points vanish
        slice -- LinearSystem, system defining a generic linear space
        witness_points -- iterable of WitnessPoint objects,  witness point set, V(f) \cap V(L)
        """
        self._system = system
        self._slice = slice
        self._witness_data = witness_data
        try:
            self._witness_points = list(witness_points)
        except TypeError:
            self._witness_points = [witness_points]
            
        wp = self._witness_points[0]
        self._component_id = wp._component_id
        for w in self._witness_points:
            if w._component_id != self._component_id:
                msg = 'WitnessPoint {0} and WitnessPoint {1} do not lie on the same component'.format(wp, w)
                raise WitnessDataException(msg)
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        sy = self._system
        sl = self._slice
        wp = self._witness_points
        repstr = 'WitnessSet({0},{1},{2})'.format(repr(sy), repr(sl), repr(wp))
        
        return repstr
    
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        wp = self._witness_points
        repstr = '{' + ', '.join([str(p) for p in wp]) + '}'
        
        return repstr
        

    @property
    def system(self):
        return self._system
    @property
    def slice(self):
        return self._slice
    @property
    def witness_data(self):
        return self._witness_data
    @property
    def witness_points(self):
        return self._witness_points