from mpmath import zeros, matrix as mpmatrix
from sympy import sympify, ShapeError, Matrix as spmatrix

from naglib.exceptions import NonPolynomialException, NonLinearException

class NAGobject(object):
    """
    A meta class. Nothing here (yet)
    """
    pass

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, dim, component_id, degree, witness_set):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        dim -- int, the dimension of the component
        component_id -- int, the component number for this dimension, as assigned by Bertini
        degree -- int, the degree of the component
        witness_set -- WitnessSet, the witness set describing this component
        """
        self._dim = dim
        self._component_id = component_id
        self._degree = degree
        self._witness_set = witness_set

    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        dim = self._dim
        cid = self._component_id
        deg = self._degree
        return '{0}-dimensional {2} component of degree {1}'.format(dim,deg,cid)

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        dim = self._dim
        cid = self._component_id
        deg = self._degree
        wst = self._witness_set
        repstr = 'IrreducibleComponent({0},{1},{2},{3})'.format(dim,cid,deg,repr(wst))
        return self.__str__()

    @property
    def degree(self):
        return self._degree
    @degree.setter
    def degree(self, deg):
        self._degree = deg

class PolynomialSystem(NAGobject):
    """
    A polynomial system
    """
    def __init__(self, polynomials, variables=None, parameters=None):
        """
        Initialize the PolynomialSystem object
        
        Keyword arguments:
        polynomials -- iterable, symbolic polynomials
        variables -- iterable, the variables in the system;
                     if None, 'variables' will be taken to be all
                     free symbols in 'polynomials'
        parameters -- iterable, the parameters in the function,
                      if the function is parameterized; if None or empty,
                      assume the function is not parameterized and
                      'parameters' becomes an empty tuple
        """
        if not hasattr(polynomials, '__iter__'):
            self._polynomials = (sympify(polynomials),)
        else:
            self._polynomials = tuple(sympify(polynomials))
            
        # check if any polynomials are actually not polynomials
        for p in self._polynomials:
            if not p.is_polynomial():
                raise(NonPolynomialException(str(p)))
            
        if variables and not hasattr(variables, '__iter__'):
            self._variables = (sympify(variables),)
        elif variables:
            self._variables = tuple(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._polynomials:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = tuple(sympify(variable_strings))
            
        if parameters and not hasattr(parameters, '__iter__'):
            self._parameters = (sympify(parameters),)
        elif parameters:
            self._parameters = tuple(sympify(parameters))
        else:
            self._parameters = ()
        
        self._degree = tuple([p.as_poly().degree() for p in self._polynomials])
        self._num_variables = len(self._variables)
        self._num_polynomials = len(self._polynomials)
        self._shape = (self._num_polynomials, self._num_variables)
            
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        polynomials = self._polynomials
        # even up the lengths of the function strings
        fstrs = [str(f) for f in polynomials]
        strlens = [len(f) for f in fstrs]
        maxlen = max(strlens)
        fstrs = [f + ' '*(maxlen - len(f)) for f in fstrs]
        fstr = '\n'.join(['[{0}]'.format(f) for f in fstrs])
        return fstr
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        polynomials  = self._polynomials
        variables  = self._variables
        parameters = self._parameters
        repstr = 'PolynomialSystem({0},{1},{2})'.format(polynomials,variables,parameters)
        return repstr
    
    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        polynomials = self._polynomials
        return polynomials[key]
    
    def __eq__(self, other):
        """
        x.__eq__(y) <==> x == y
        """
        if not isinstance(other, PolynomialSystem):
            return False
        elif isinstance(other, LinearSystem) and self._parameters:
            return False
        elif isinstance(other, LinearSystem):
            spoly = self._polynomials
            opoly = other._polynomials
            svars = self._variables
            ovars = other._variables
            return (spoly == opoly and svars == ovars)
        else:
            spoly = self._polynomials
            opoly = other._polynomials
            svars = self._variables
            ovars = other._variables
            spara = self._parameters
            opara = other._parameters
        return (spoly == opoly and svars == ovars and spara == opara)
    
    def lmul(self, other):
        """
        x.lmul(y) <==> y*x
        
        Because SymPy and mpmath don't return NotImplemented, can't
        override __rmul__
        """
        # multiply by an mpmatrix on the right
        num_rows = self._shape[0]
        polys = self._polynomials
        if isinstance(other, mpmatrix):
            if other.cols != num_rows:
                raise(ValueError('dimensions not compatible for multiplication'))
            res_polys = []
            for i in range(other.rows):
                res = sympify(0)
                for j in range(num_rows):
                    res += polys[j]*other[i,j]
                res_polys.append(res)
            
            return PolynomialSystem(res_polys, self._variables, self._parameters)
        elif isinstance(other, spmatrix):
            if other.cols != num_rows:
                raise(ShapeError('Matrices size mismatch'))
            res_polys = []
            for i in range(other.rows):
                res = sympify(0)
                for j in range(num_rows):
                    res += polys[j]*other[i,j]
                res_polys.append(res)
            
            return PolynomialSystem(res_polys)
        
    @property
    def polynomials(self):
        return self._polynomials
    @property
    def variables(self):
        return self._variables
    @property
    def degree(self):
        return self._degree
    @property
    def shape(self):
        return self._shape
    
class LinearSystem(PolynomialSystem):
    """
    A linear system
    
    !!!Use this only for slicing!!!
    """
    def __init__(self, polynomials, variables=None):
        """
        Initialize the LinearSystem object
        
        Keyword arguments:
        polynomials -- iterable, symbolic polynomials
        variables -- iterable, the variables in the system;
                     if None, 'variables' will be taken to be all
                     free symbols in 'polynomials'
        """
        if not hasattr(polynomials, '__iter__'):
            self._polynomials = (sympify(polynomials),)
        else:
            self._polynomials = tuple(sympify(polynomials))
            
        # check if any polynomials are actually not polynomials
        for p in self._polynomials:
            if not p.is_polynomial():
                raise(NonPolynomialException(str(p)))
            elif p.as_poly().degree() > 1:
                raise(NonLinearException(str(p)))
            
        if variables and not hasattr(variables, '__iter__'):
            self._variables = (sympify(variables),)
        elif variables:
            self._variables = tuple(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._polynomials:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = tuple(sympify(variable_strings))
            
        self._num_variables = len(self._variables)
        self._num_polynomials = len(self._polynomials)
        self._shape = (self._num_polynomials, self._num_variables)
            
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        polynomials = self._polynomials
        # even up the lengths of the function strings
        fstrs = [str(f) for f in polynomials]
        strlens = [len(f) for f in fstrs]
        maxlen = max(strlens)
        fstrs = [f + ' '*(maxlen - len(f)) for f in fstrs]
        fstr = '\n'.join(['[{0}]'.format(f) for f in fstrs])
        return fstr
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        polynomials  = self._polynomials
        variables  = self._variables
        repstr = 'LinearSystem({0},{1})'.format(polynomials,variables)
        return repstr
    
    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        polynomials = self._polynomials
        return polynomials[key]
    
    def __eq__(self, other):
        """
        x.__eq__(y) <==> x == y
        """
        if not isinstance(other, PolynomialSystem):
            return False
        elif isinstance(other, PolynomialSystem) and other._parameters:
            return False
        else:
            spoly = self._polynomials
            opoly = other._polynomials
            svars = self._variables
            ovars = other._variables
        return (spoly == opoly and svars == ovars)
    
    def lmul(self, other):
        """
        x.lmul(y) <==> y*x
        
        Because SymPy and mpmath don't return NotImplemented, can't
        override __rmul__
        """
        # multiply by an mpmatrix on the right
        num_rows = self._shape[0]
        polys = self._polynomials
        if isinstance(other, mpmatrix):
            if other.cols != num_rows:
                raise(ValueError('dimensions not compatible for multiplication'))
            res_polys = []
            for i in range(other.rows):
                res = sympify(0)
                for j in range(num_rows):
                    res += polys[j]*other[i,j]
                res_polys.append(res)
            
            return LinearSystem(res_polys, self._variables)
        elif isinstance(other, spmatrix):
            if other.cols != num_rows:
                raise(ShapeError('Matrices size mismatch'))
            res_polys = []
            for i in range(other.rows):
                res = sympify(0)
                for j in range(num_rows):
                    res += polys[j]*other[i,j]
                res_polys.append(res)
            
            try:
                res = LinearSystem(res_polys)
            except NonLinearException:
                res = PolynomialSystem(res_polys)
            
            return res
        
    def matrix(self):
        """
        Create an mpmath matrix from the LinearSystem
        """
        
        polynomials = self._polynomials
        variables = self._variables
        shape = self._shape
        
        res = zeros(*shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                res[i,j] = polynomials[i].coeff(variables[j])
                
        return res
    
    @property
    def polynomials(self):
        return self._polynomials
    @property
    def variables(self):
        return self._variables
    @property
    def shape(self):
        return self._shape
    
class WitnessPoint(NAGobject):
    """
    A single witness point for an irreducible component
    """
    def __init__(self, dim, component_id, pt, isprojective=True):
        """
        Initialize the WitnessPoint object.

        Keyword arguments:
        dim -- int, the dimension of the component to which the point belongs
        component_id -- int, the component number for this dimension, as assigned by Bertini
        pt -- numeric iterable, the coordinates of the witness point
        isprojective -- boolean, True if the point is projective, otherwise False
        """
        self._dim = dim
        self._component_id = component_id
        if not hasattr(pt, '__iter__'):
            self._pt = mpmatrix((pt,))
        else:
            self._pt = mpmatrix(pt)
        self._isprojective = isprojective

    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        return '(' + ','.join(str(p) for p in self._pt) + ')'

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        dim = self._dim
        cid = self._component_id
        pt  = self._pt
        prj = self._isprojective
        repstr = 'WitnessPoint({0},{1},{2},{3})'.format(dim,cid,repr(pt),prj)
        return repstr
    
    def dehomogenize(self):
        """
        Dehomogenize the witness point
        
        Returns a new WitnessPoint object
        """
        dim = self._dim
        component_id = self._component_id
        pt = self._pt
        if not self._isprojective:
            newpt = WitnessPoint(dim, component_id, pt, False)
        else:
            first = pt[0]
            pt = pt[1:]/first # first should not be equal to 0 if homogenous
            newpt = WitnessPoint(dim, component_id, pt, False)
            
        return newpt

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

class WitnessSet(NAGobject):
    """
    A witness set for a component
    """
    def __init__(self, system, slice, witness_points, isprojective=True):
        """Initialize the WitnessSet
        
        Keyword arguments:
        system -- PolynomialSystem, system on which all the points vanish
        slice -- LinearSystem, system defining a generic linear space
        witness_points -- iterable of WitnessPoint objects,  witness point set, V(f) \cap V(L)
        isprojective -- boolean, True if the witness set describes a projective set
        """
        self._system = system
        self._slice = slice
        self._witness_points = set(witness_points)
        self._isprojective = isprojective
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        sy = self._system
        sl = self._slice
        wp = self._witness_points
        ip = self._isprojective
        repstr = 'WitnessSet({0},{1},{2},{3})'.format(repr(sy),repr(sl),repr(wp),ip)
        
        return repstr
    
    @property
    def system(self):
        return self._system
    @property
    def slice(self):
        return self._slices
    @property
    def witness_points(self):
        return self._witness_points
    @property
    def isprojective(self):
        return self._isprojective