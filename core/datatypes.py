from sympy import sympify

class NAGobject(type):
    """
    A meta class. Nothing here (yet)
    """
    pass

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    
    Consists of 
    """
    def __init__(self, dim, component_id, degree, witness_set):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        dim -- the dimension of the component
        component_id -- the component number for this dimension, as assigned by Bertini
        degree -- the degree of the component
        witness_set -- the witness set (type WitnessSet) describing this component
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
        repstr = 'IrreducibleComponent({0},{1},{2},{3})'.format(dim,cid,deg,wst)
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
    def __init__(self, functions, variables=None, parameters=None):
        """Initialize the PolynomialSystem object
        
        Keyword arguments:
        functions -- an iterable of symbolic function arguments
        variables -- an iterable of the variables in the function;
                     if None, 'variables' will be taken to be all
                     free symbols in 'functions'
        parameters -- an iterable of the parameters in the function,
                      if the function is parameterized; if None or empty,
                      I assume the function is not parameterized and
                      parameters becomes an empty tuple
        """
        if not hasattr(functions, '__iter_'):
            self._functions = (sympify(functions),)
        else:
            self._functions = tuple(sympify(functions))
        if variables and not hasattr(variables, '__iter__'):
            self._variables = (sympify(variables),)
        elif variables:
            self._variables = tuple(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._functions:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = sympify(variable_strings)            
        if parameters and not hasattr(parameters, '__iter__'):
            self._parameters = (sympify(parameters),)
        elif parameters:
            self._parameters = tuple(sympify(parameters))
        else:
            self._parameters = ()
            
    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        functions = self._functions
        # even up the lengths of the function strings
        fstrs = [str(f) for f in functions]
        strlens = [len(f) for f in fstrs]
        maxlen = max(strlens)
        fstrs = [f + ' '*(maxlen - len(f)) for f in fstrs]
        fstr = '\n'.join(['[{0}]'.format(f) for f in fstrs])
    
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        functions  = self._functions
        variables  = self._variables
        parameters = self._parameters
        repstr = 'PolynomialSystem({0},{1},{2})'.format(functions,variables,parameters)
        return repstr
    
    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        functions = self._functions
        return functions[key]
    
class LinearSystem(PolynomialSystem):
    """
    A linear system
    """
    def __init__(self, functions, variables=None):
        """Initialize the PolynomialSystem object
        
        Keyword arguments:
        functions -- an iterable of symbolic function arguments
        variables -- an iterable of the variables in the function;
                     if None, 'variables' will be taken to be all
                     free symbols in 'functions'
        """
        if not hasattr(functions, '__iter_'):
            self._functions = (sympify(functions),)
        else:
            self._functions = tuple(sympify(functions))
        if variables and not hasattr(variables, '__iter__'):
            self._variables = (sympify(variables),)
        elif variables:
            self._variables = tuple(sympify(variables))
        else: # variables not specified
            variable_list = set()
            for f in self._functions:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list)
            variable_strings = sorted([str(v) for v in variable_list])
            self._variables = sympify(variable_strings)
        
    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        functions  = self._functions
        variables  = self._variables
        repstr = 'LinearSystem({0},{1},{2})'.format(functions,variables,parameters)
        return repstr
    
class WitnessPoint(NAGobject):
    """
    A single witness point for an irreducible component
    """
    def __init__(self, dim, component_id, pt, isprojective=True):
        """
        Initialize the WitnessPoint object.

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
        """
        x.__str__() <==> str(x)
        """
        return '(' + ', '.join(str(p) for p in self._pt) + ')'

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        dim = self._dim
        cid = self._component_id
        pt  = self._pt
        prj = self._isprojective
        repstr = 'WitnessPoint({0},{1},{2},{3})'.format(dim,cid,pt,prj)
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
            pt = pt[1:]
            # first should not be equal to 0
            pt = [pt[i]/first for i in range(len(pt))]
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
    def __init__(self, f, L, W, isprojective=True):
        """Initialize the WitnessSet
        
        Keyword arguments:
        f -- the system (type PolynomialSystem) on which all the points vanish
        L -- the system (type LinearSystem) defining a generic linear space
        W -- the witness point set, V(f) \cap V(L)
        """
        self._f = f
        self._L = L
        self._W = W