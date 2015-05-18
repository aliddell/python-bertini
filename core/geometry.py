from __future__ import print_function

from sympy import Matrix

from naglib.bertini.sysutils import BertiniRun
from naglib.exceptions import BertiniError
from naglib.core.base import NAGobject

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, witness_set, codim, component_id, **kwargs):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        witness_set  -- WitnessSet, the witness set describing this component
        codim        -- int, the codimension of the component
        component_id -- int, the component number for this dimension
        
        Optional keyword arguments:
        """
        self._witness_set = witness_set
        self._codim = codim
        self._component_id = component_id
        self._degree = len(witness_set.witness_points)
        
        # optional keyword arguments
        kkeys = kwargs.keys()
        
        if 'randomization_matrix' in kkeys:
            self._randomization_matrix = kwargs['randomization_matrix']
        else:
            self._randomization_matrix = Matrix([])
            
        if 'homogenization_matrix' in kkeys:
            self._homogenization_matrix = kwargs['homogenization_matrix']
        else:
            self._homogenization_matrix = Matrix([])
            
        if 'homogenization_vector' in kkeys:
            self._homogenization_vector = kwargs['homogenization_vector']
        else:
            self._homogenization_vector = Matrix([])
            
        if 'homogenization_variable' in kkeys:
            self._homogenization_variable = kwargs['homogenization_variable']
        else:
            self._homogenization_variable = None
            
        if 'patch_coefficients' in kkeys:
            self._patch_coefficients = kwargs['patch_coefficients']
        else:
            self._patch_coefficients = Matrix([])

    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        codim = self._codim
        cid = self._component_id
        deg = self._degree
        return '{0}-dimensional irreducible component ({2}) of degree {1}'.format(codim,deg,cid)

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        codim = self._codim
        cid = self._component_id
        wst = self._witness_set
        repstr = 'IrreducibleComponent({0},{1},{2})'.format(repr(wst),codim,cid)
        return repstr
        
    def __eq__(self, other):
        """
        x.__eq__(y) <==> x == y
        """
        if type(other) != IrreducibleComponent:
            return False
        
        sp = self.sample(1)
        op = other.sample(1)
        
        sco = self.contains(op)
        ocs = other.contains(sp)
        
        return sco and ocs
        
    def _construct_witness_data(self):
        pass
        
    def contains(self, other):
        """
        Return True if self contains other
        """
        if type(other) not in (list, tuple):
            msg = "cannot understand data type"
            raise TypeError(msg)
            
        system = self._witness_set.system
        test_run = BertiniRun(system,
                              tracktype=BertiniRun.TMEMTEST,
                              component=self,
                              start=other)
        return test_run.run()
    
    def sample(self, numpoints=1, usebertini=True):
        """
        Sample points from self
        """
        if numpoints < 1:
            msg = "sample at least one point"
            raise BertiniError(msg)
        system  = self.witness_set.system
        
        points = None
        if usebertini:
            sample_run = BertiniRun(system, BertiniRun.TSAMPLE, sample=numpoints, component=self)
            points = sample_run.run()
        else:
            msg = "nothing to use yet but Bertini"
            raise NotImplementedError(msg)
        
        return points
    
    @property
    def codim(self):
        return self._codim
    @property
    def component_id(self):
        return self._component_id
    @property
    def degree(self):
        return self._degree
    @property
    def dim(self):
        variables = self.system.variables
        return len(variables) - self._codim
    @property
    def homogenization_matrix(self):
        return self._homogenization_matrix
    @property
    def homogenization_variable(self):
        return self._homogenization_variable
    @property
    def homogenization_vector(self):
        return self._homogenization_vector
    @property
    def is_projective(self):
        return not not self._homogenization_variable
    @property
    def patch_coefficients(self):
        return self._patch_coefficients
    @property
    def randomization_matrix(self):
        return self._randomization_matrix
    @property
    def system(self):
        return self._witness_set.system
    @property
    def witness_data(self):
        return self._witness_set._witness_data
    @property
    def witness_set(self):
        return self._witness_set