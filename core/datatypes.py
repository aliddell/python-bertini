from mpmath import mpc, rand, matrix as mpmatrix, zeros as mpzeros
from sympy import sympify, ShapeError, Matrix as spmatrix, zeros as spzeros

from naglib.exceptions import NonPolynomialException, NonLinearException, AffineException

class NAGobject(object):
    """
    A meta class. Nothing here (yet)
    """
    pass

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, dim, component_id, degree, witness_set, isprojective=True):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        dim          -- int, the dimension of the component
        component_id -- int, the component number for this dimension, as
                        assigned by Bertini
        degree       -- int, the degree of the component
        witness_set  -- WitnessSet, the witness set describing this component
        isprojective -- optional boolean, True if the system is projective,
                        otherwise, False (default: True)
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
        repstr = 'IrreducibleComponent(\n{0},\n{1},\n{2},\n{3})'.format(dim,
                                                                        cid,
                                                                        deg,
                                                                        repr(wst))
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
    def __init__(self, polynomials, variables=None, parameters=None, isprojective=True):
        """
        Initialize the PolynomialSystem object
        
        Keyword arguments:
        polynomials  -- symbolic iterable, polynomials
        variables    -- optional symbolic iterable, the variables in the
                        system; if None or empty, 'variables' will be taken to
                        be all free symbols in 'polynomials' (default: None)
        parameters   -- optional symbolic iterable, the parameters in the
                        system, if it is parameterized; if None or empty,
                        assume the system is not parameterized and
                        'parameters' becomes None (default: None)
        isprojective -- optional boolean, True if the system is projective,
                        otherwise, False (default: True)
        """
        
        self._polynomials = spmatrix(sympify(polynomials))
            
        # check if any polynomials are actually not polynomials
        for p in self._polynomials:
            if not p.is_polynomial():
                raise(NonPolynomialException(str(p)))
            
        if variables:
            self._variables = spmatrix(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._polynomials:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = spmatrix(sympify(variable_strings))
            
        if parameters:
            self._parameters = spmatrix(sympify(parameters))
        else:
            self._parameters = None
        
        d = []
        for p in self._polynomials:
            if p.is_number:
                d.append(0)
            else:
                d.append(p.as_poly().degree())        
        self._degree = tuple(d)
        self._num_variables = len(self._variables)
        self._num_polynomials = len(self._polynomials)
        self._shape = (self._num_polynomials, self._num_variables)
            
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        polynomials = self._polynomials
        # even up the lengths of the polynomial strings
        fstrs = [str(f) for f in polynomials]
        strlens = [len(f) for f in fstrs]
        maxlen = max(strlens)
        fstrs = [' '*(maxlen - len(f)) + f for f in fstrs]
        fstr = '\n'.join(['[{0}]'.format(f) for f in fstrs])
        return fstr
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        polynomials  = self._polynomials
        variables  = self._variables
        parameters = self._parameters
        repstr = 'PolynomialSystem(\n{0},\n{1},\n{2})'.format(repr(polynomials),
                                                              repr(variables),
                                                              repr(parameters))
        
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
            res_polys = list(other * self._polynomials)
            res = PolynomialSystem(res_polys)
            
            return res
    
    def jacobian(self):
        """
        Returns the Jacobian, the polynomial system, and the variables,
        all as symbolic matrices
        """
        variables = self._variables
        polynomials = self._polynomials
        num_polynomials,num_variables = self._shape
        jac = spzeros(num_polynomials,num_variables)
        for i in range(num_polynomials):
            for j in range(num_variables):
                jac[i,j] = polynomials[i].diff(variables[j])
        
        return jac,polynomials,variables
    
    def rank(self):
        """
        Return a numeric value, the rank of the Jacobian at
        a 'generic' point.
        """
        polynomials = self._polynomials
        parameters = self._parameters
        variables = self._variables
        if parameters:
            allvars = variables.col_join(parameters)
        else:
            allvars = variables
        
        varsubs = mpzeros(len(allvars), 1)
        # compute sufficiently generic complex points
        for i in range(len(varsubs)):
            # rand()/rand() can vary magnitude satisfactorily
            try:
                re = rand()/rand()
                im = rand()/rand()
            except ZeroDivisionError:
                # try again
                return self.rank()
            varsubs[i] = mpc(re, im)
            
        jac = self.jacobian()[0]
        jac = jac.subs(zip(variables, varsubs))
        
        return jac.rank()
        
    def subs(self, *args, **kwargs):
        """
        Return a new PolynomialSystem with subs applied to
        each entry of 'polynomials'
        """
        polynomials = self._polynomials
        psubs = polynomials.subs(*args, **kwargs)
        ps = PolynomialSystem(psubs)
        
        return ps        
        
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
    def __init__(self, polynomials, variables=None, isprojective=True):
        """
        Initialize the LinearSystem object
        
        Keyword arguments:
        polynomials  -- symbolic iterable, the linear polynomials making
                        up the system
        variables    -- optional symbolic iterable, the variables in
                        the system; if None, 'variables' will be taken to
                        be all free symbols in 'polynomials' (default: None)
        isprojective -- optional boolean, True if the system is projective,
                        otherwise False (default: True)
        """
        self._polynomials = spmatrix(sympify(polynomials))
            
        # check if any polynomials are actually not polynomials
        for p in self._polynomials:
            if not p.is_polynomial():
                raise(NonPolynomialException(str(p)))
            elif p.as_poly().degree() > 1:
                raise(NonLinearException(str(p)))
            
        if variables:
            self._variables = spmatrix(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._polynomials:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = spmatrix(sympify(variable_strings))
        
        # check if any polynomial is affine
        for p in self._polynomials:
            psub = p.subs(zip(self._variables, spzeros(*self._variables.shape)))
            if psub != 0:
                raise(AffineException(str(p)))
            
        self._num_variables = len(self._variables)
        self._num_polynomials = len(self._polynomials)
        self._shape = (self._num_polynomials, self._num_variables)
            
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        polynomials = self._polynomials
        # even up the lengths of the polynomial strings
        fstrs = [str(f) for f in polynomials]
        strlens = [len(f) for f in fstrs]
        maxlen = max(strlens)
        fstrs = [' '*(maxlen - len(f)) + f for f in fstrs]
        fstr = '\n'.join(['[{0}]'.format(f) for f in fstrs])
        return fstr
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        polynomials  = self._polynomials
        variables  = self._variables
        repstr = 'LinearSystem(\n{0},\n{1})'.format(polynomials,variables)
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
            res_polys = list(other * self._polynomials)
            try:
                res = LinearSystem(res_polys)
            except NonLinearException:
                res = PolynomialSystem(res_polys)
            
            return res
        
    def matrix(self, symbolic=False):
        """
        Create an mpmath matrix from the LinearSystem
        
        Keyword arguments:
        symbolic -- optional boolean, return an mpmath matrix if False,
                    a SymPy matrix if True (default: False)
        """
        
        polynomials = self._polynomials
        variables = self._variables
        shape = self._shape
        
        if symbolic:
            res = spzeros(*shape)
        else:
            res = mpzeros(*shape)

        for i in range(shape[0]):
            for j in range(shape[1]):
                res[i,j] = polynomials[i].coeff(variables[j])
                
        return res
    
    def rank(self):
        """
        Return the rank of the linear system
        """
        m = self.matrix(symbolic=True)
        return m.rank()
    
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
        dim          -- int, the dimension of the component to which the point
                        belongs
        component_id -- int, the component number for this dimension, as
                        assigned by Bertini
        pt           -- numeric iterable, the coordinates of the witness point
        isprojective -- optional boolean, True if the system is projective,
                        otherwise False (default: True)
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
        repstr = 'WitnessPoint(\n{0},\n{1},\n{2},\n{3})'.format(dim,
                                                                cid,
                                                                repr(pt),
                                                                prj)
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