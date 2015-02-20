from __future__ import print_function

from mpmath import mpc, rand, matrix as mpmatrix, zeros as mpzeros
from sympy import sympify, ShapeError, Matrix as spmatrix, zeros as spzeros

from naglib.exceptions import BertiniError, NonPolynomialException, NonLinearException, AffineException, NonHomogeneousException

class NAGobject(object):
    """
    A meta class. Nothing here (yet)
    """
    pass

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, system, dim, component_id, witness_set, dirname, isprojective=False):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        system       -- PolynomialSystem, defines the algebraic set of which
                        self is an irreducible component
        dim          -- int, the dimension of the component
        component_id -- int, the component number for this dimension, as
                        assigned by Bertini
        witness_set  -- WitnessSet, the witness set describing this component
        isprojective -- optional boolean, True if the system is projective,
                        otherwise, False (default: True)
        """
        self._system = system
        self._dim = dim
        self._component_id = component_id
        self._witness_set = witness_set
        self._degree = len(witness_set.witness_points)
        self._dirname = dirname
        self._isprojective = isprojective

    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        dim = self._dim
        cid = self._component_id
        deg = self._degree
        return '{0}-dimensional component ({2}) of degree {1}'.format(dim,deg,cid)

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
        return repstr
    
    def __eq__(self, other):
        """
        x.__eq__(y) <==> x == y
        """
        if not isinstance(other, IrreducibleComponent):
            return False
        
        ssys = self._system
        sdim = self._dim
        scid = self._component_id
        osys = other.system
        odim = other.dim
        ocid = other.component_id
        eq = (ssys == osys and sdim == odim and scid == ocid)
        
        return eq
    
    def sample(self, numpoints=1, usebertini=True):
        """
        Sample a point from self
        """
        dim     = self._dim
        comp_id = self._component_id
        system  = self._system
        
        points = None
        if usebertini:
            from naglib.bertini.fileutils import write_input, read_points
            from naglib.bertini.sysutils import call_bertini as call
            
            dirname = self._dirname
            inputfile = dirname + '/input'
            instructions = dirname + '/instructions'
            sampled = dirname + '/sampled'
            
            # instructions to Bertini (redirect stdin)
            fh = open(instructions, 'w')
            print('{0}'.format(dim), file=fh) # sample from dimension 'dim'
            print('{0}'.format(comp_id), file=fh) # sample from component 'comp_id'
            print('{0}'.format(numpoints), file=fh) # sample 'numpoints' points
            print('0', file=fh) # write point to a file
            print(sampled, file=fh) # write to file 'sampled'
            fh.close()
        
            config = {'filename':inputfile, 'TrackType':2}
            write_input(system=system, config=config)
            call(input_file=inputfile, stdin=instructions)
            points = read_points(sampled)
        else:
            raise(NotImplementedError('nothing to use but bertini yet'))
        
        return points
    
    @property
    def system(self):
        return self._system
    @property
    def degree(self):
        return self._degree
    @property
    def dim(self):
        return self._dim
    @property
    def component_id(self):
        return self._component_id
    @property
    def witness_set(self):
        return self._witness_set
    
class Point(NAGObject):
    """
    A point in affine or projective space
    """
    def __init__(self, coordinates, isprojective=False):
        """
        Initialize the Point object
        """
        # ensure not actually a Point object
        if isinstance(coordinates, Point):
            return coordinates
        
        if not hasattr(coordinates, '__iter__'):
            coordinates = [coordinates]
        else:
            coordinates = list(coordinates)
        
        coordinates = [sympify(c) for c in coordinates]
        self._coordinates = spmatrix(coordinates)
        self._isprojective = isprojective
        if isprojective:
            if coordinates[0] == 0:
                raise(AffineException("projective points can't be 0 in the first entry"))
            self._dim = len(self._coordinates) - 1
        else:
            self._dim = len(self._coordinates)
        

class PolynomialSystem(NAGobject):
    """
    A polynomial system
    """
    def __init__(self, polynomials, variables=None, parameters=None, isprojective=False):
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
                        'parameters' becomes () (default: None)
        isprojective -- optional boolean, True if the system is projective,
                        otherwise, False (default: True)
        """
        # check if 'polynomials' is a string and/or contains '^' for exp
        from re import sub as resub

        if type(polynomials) == str:
            polynomials = [resub(r'\^', r'**', polynomials)]

        self._polynomials = spmatrix(sympify(polynomials))
            
        # check if any polynomials are actually not polynomials
        for p in self._polynomials:
            if not p.is_polynomial():
                raise(NonPolynomialException(str(p)))
            
        if parameters:
            self._parameters = spmatrix(sympify(parameters))
        else:
            self._parameters = ()
        
        if variables:
            self._variables = spmatrix(sympify(variables))
        else: # variables not specified
            variable_list = set()
            param_set = set(self._parameters)
            for f in self._polynomials:
                variable_list = variable_list.union(f.free_symbols)
            variable_list = list(variable_list - param_set)
            variable_strings = [str(v) for v in variable_list]
            variable_strings.sort()
            self._variables = spmatrix(sympify(variable_strings))
            
        self._isprojective = isprojective
        
        d = []
        params = list(self._parameters)
        ones = [1 for param in params]
        paramsubs = zip(params,ones)
        
        # keep parameters out of degree calculation
        for polyn in self._polynomials:
            p = polyn.subs(paramsubs)
            if p.is_number:
                deg = 0
            else:
                deg = p.as_poly().total_degree()
            # check if polynomial is homogeneous
            if isprojective:
                terms = p.as_ordered_terms()
                for t in terms:
                    if t.is_number and deg != 0:
                        raise(NonHomogeneousException(p))
                    elif t.is_number:
                        pass
                    elif t.as_poly().degree() != deg:
                        raise(NonHomogeneousException(p))
            d.append(deg)
            
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
        variables    = self._variables
        parameters   = self._parameters
        isprojective = self._isprojective
        repstr = 'PolynomialSystem(\n{0},\n{1},\n{2}\n{3})'.format(repr(polynomials),
                                                                   repr(variables),
                                                                   repr(parameters),
                                                                   repr(isprojective))
        
        return repstr
    
    def __getitem__(self, key):
        """
        x.__getitem__(y) <==> x[y]
        """
        polynomials = self._polynomials
        return polynomials[key]
    
    def __add__(self, other):
        """
        x.__add__(y) <==> x + y
        """
        spolct = self._shape[0]
        opolct = other.shape[0]
        isprojective = self._isprojective
        if not isinstance(other, PolynomialSystem):
            othertype = type(other)
            errmsg = "unsupported operand type(s) for +: " \
                "'PolynomialSystem' and '{0}'".format(othertype)
            raise(TypeError(errmsg))
        elif other.isprojective != isprojective:
            errmsg = 'attempting to add systems of different types'
            raise(ValueError(errmsg))
        elif spolct != opolct:
            errmsg = "don't know how to add systems of different sizes"
            raise(ValueError(errmsg))
        else:
            svariables = self._variables
            sparams    = self._parameters
            spoly      = self._polynomials
            ovariables = other.variables
            oparams    = other.parameters
            opoly      = other.polynomials
            
            variables  = set(svariables).union(set(ovariables))
            parameters = set(sparams).union(set(oparams))
            variables  = [str(v) for v in variables]
            parameters = [str(p) for p in parameters]
            variables.sort()
            parameters.sort()
            
            variables  = sympify(variables)
            parameters = sympify(parameters)
            
            polynomials = spoly + opoly
            
            return PolynomialSystem(polynomials, variables, parameters, isprojective)
        
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
        
    def assign_parameters(self, params):
        """
        Set params as parameters in self
        """
        if not hasattr(params, '__iter__'):
            params = [params]
        params = set(sympify(params))
        
        sparams = set(self._parameters).union(params)
        svars   = set(self._variables) - sparams
        
        str_pars = sorted([str(p) for p in sparams])
        str_vars = sorted([str(v) for v in svars])
        
        self._parameters = spmatrix(sympify(str_pars))
        self._variables  = spmatrix(sympify(str_vars))
        
        # update shape
        self._shape = (len(self._polynomials), len(self._variables))
        
    def homogenize(self):
        """
        Homogenize the system
        
        If already homogeneous, return self
        """
        polynomials = self._polynomials
        polynomials = [p.expand() for p in polynomials]
        variables = self._variables
        parameters = self._parameters
        degree = self._degree
        
        if self._isprojective:
            return self
        
        # create a new homogenizing variable
        p0 = 'p0'
        while p0 in variables:
            p0 += '0'
            # shouldn't last too long
        p0 = spmatrix([sympify(p0)])
        
        homvars = variables.row_insert(0, p0)
        polyterms = [p.as_ordered_terms() for p in polynomials]
        hompolys = []
        
        p0 = p0[0]
        for i in range(len(polyterms)):
            pdeg = degree[i]
            polyterm = polyterms[i]
            hompolyterm = []
            for term in polyterm:
                if term.is_number:
                    hompolyterm.append(term*p0**pdeg)
                else:
                    tdeg = term.as_poly().degree()
                    hompolyterm.append(term*p0**(pdeg-tdeg))
            hompolys.append(sum(hompolyterm))
        
        hompolys = spmatrix(hompolys)
        return PolynomialSystem(hompolys, homvars, parameters, isprojective=True)
    
    def dehomogenize(self, homvar=None):
        """
        Dehomogenize the system
        
        If already nonhomogeneous, return self
        """
        hompolys = self._polynomials
        hompolys = spmatrix([p.expand() for p in hompolys])
        homvars = self._variables
        parameters = self._parameters
        degree = self._degree
        
        if not self._isprojective:
            return self
        
        if not homvar:
            homvar = homvars[0]
            variables = spmatrix(homvars[1:])
        elif homvar in homvars:
            homvars = list(homvars)
            dex = homvars.index(homvar)
            homvars.pop(dex)
            variables = spmatrix(homvars)
        else:
            raise(ValueError('homogenizing variable {0} not found'.format(homvar)))
        
        polynomials = hompolys.subs({homvar:1})
        
        return PolynomialSystem(polynomials, variables, parameters, isprojective=False)
        
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
                real = rand()/rand()
                imag = rand()/rand()
            except ZeroDivisionError:
                # try again
                return self.rank()
            varsubs[i] = mpc(real, imag)
            
        jac = self.jacobian()[0]
        jac = jac.subs(zip(variables, varsubs))
        
        # TODO: allow user to specify tolerance (what is 'zero')
        return jac.rank()
    
    def solve(self, start_params=None, final_params=None, start=None, usebertini=True):
        """
        Solve the system. If non-square, return the NID
        
        If the system has parameters and you don't supply
        final parameters, solve the system with random
        start parameters, return the solutions and the start parameters.
        If you do supply final parameters, solve the system with 
        random start parameters, then solve again and return the
        start parameters
        """
        rank = self.rank()
        polynomials = self._polynomials
        variables   = self._variables
        parameters  = self._parameters
        
        if usebertini:
            from tempfile import mkdtemp
            from naglib import TEMPDIR as basedir
            from naglib.bertini.sysutils import call_bertini
            from naglib.bertini.fileutils import write_input, read_points, fprint
            from naglib.bertini.data import compute_NID
            
            if rank != len(polynomials) or rank != len(variables):
                return compute_NID(self)
            elif parameters and not start_params and not final_params:
                dirname  = mkdtemp(prefix=basedir)
                filename = dirname + '/input'
                config   = {'filename': filename,
                            'TrackType': 0,
                            'ParameterHomotopy':1}
                
                write_input(self, config)
                call_bertini(filename)
                points = read_points(dirname + '/finite_solutions', as_set=False)
                sp = read_points(dirname + '/start_parameters', as_set=False)
                
                return points, sp
            elif parameters and not start_params:
                dirname  = mkdtemp(prefix=basedir)
                filename = dirname + '/input'
                config   = {'filename': filename,
                            'TrackType': 0,
                            'ParameterHomotopy':1}
                write_input(self, config)
                
                call_bertini(filename)
                
                start_params = read_points(dirname + '/start_parameters', as_set=False)
                
                config   = {'filename': filename,
                            'TrackType': 0,
                            'ParameterHomotopy':2}
                write_input(self, config)
                
                fprint(final_params, dirname + '/final_parameters')
                
                call_bertini(filename)
                
                points = read_points(dirname + '/finite_solutions', as_set=False)
                return points, start_params
            elif parameters and not final_params:
                return self.subs(zip(parameters, start_params)).solve()
                #dirname  = mkdtemp(prefix=basedir)
                #filename = dirname + '/input'
                #config   = {'filename': filename,
                            #'TrackType': 0,
                            #'ParameterHomotopy':1}
                #write_input(self, config)
                
                #call_bertini(filename)
                
                #start_params = read_points(dirname + '/start_parameters', as_set=False)
                
                #config   = {'filename': filename,
                            #'TrackType': 0,
                            #'ParameterHomotopy':2}
                #write_input(self, config)
                
                #fprint(start_params, dirname + '/final_parameters')
                
                #call_bertini(filename)
                
                #points = read_points(dirname + '/finite_solutions', as_set=False)
                #return points, start_params
            elif parameters and not start:
                raise(BertiniError("'start' does not exist!!!"))
            elif parameters: # and start_params and final_params
                dirname  = mkdtemp(prefix=basedir)
                filename = dirname + '/input'
                config   = {'filename': filename,
                            'TrackType': 0,
                            'ParameterHomotopy':2}
                
                fprint(start, dirname + '/start')
                fprint(start_params, dirname + '/start_parameters')
                fprint(final_params, dirname + '/final_parameters')

                write_input(self, config)
                call_bertini(filename)
                points = read_points(dirname + '/finite_solutions', as_set=False)
                
                return points
            else:
                dirname  = mkdtemp(prefix=basedir)
                filename = dirname + '/input'
                config   = {'filename': filename,
                            'TrackType': 0}
                
                write_input(self, config)
                call_bertini(filename)
                points = read_points(dirname + '/finite_solutions', as_set=False)
                
                return points
        else:
            raise(NotImplementedError('nothing to use but bertini yet'))
        
    def subs(self, *args, **kwargs):
        """
        Return a new PolynomialSystem with subs applied to
        each entry of 'polynomials'
        """
        polynomials = self._polynomials
        psubs = polynomials.subs(*args, **kwargs)
        ps = PolynomialSystem(psubs, isprojective=self._isprojective)
        
        return ps
        
    @property
    def polynomials(self):
        return self._polynomials
    @property
    def variables(self):
        return self._variables
    @property
    def parameters(self):
        return self._parameters
    @property
    def isprojective(self):
        return self._isprojective
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
    def __init__(self, polynomials, variables=None, isprojective=False):
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
                errmsg = 'polynomial {0} is not linear'.format(p)
                raise(AffineException(errmsg))
            
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
    def __init__(self, dim, component_id, pt, isprojective=False):
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
    def __init__(self, system, slice, witness_points, isprojective=False):
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
        repstr = 'WitnessSet(\n{0},\n{1},\n{2},\n{3})'.format(repr(sy),
                                                              repr(sl),
                                                              repr(wp),
                                                              ip)
        
        return repstr
    
    @property
    def system(self):
        return self._system
    @property
    def slice(self):
        return self._slice
    @property
    def witness_points(self):
        return self._witness_points
    @property
    def isprojective(self):
        return self._isprojective
