from __future__ import print_function

from naglib.bertini.sysutils import BertiniRun
from naglib.exceptions import BertiniError
from naglib.startup import TEMPDIR as basedir
from naglib.core.base import NAGobject, Point, AffinePoint, ProjectivePoint
from naglib.core.witnessdata import WitnessPoint, WitnessSet

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, witness_set, codim, component_id):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        witness_set  -- WitnessSet, the witness set describing this component
        codim        -- int, the codimension of the component
        component_id -- int, the component number for this dimension
        """
        self._witness_set = witness_set
        self._codim = codim
        self._component_id = component_id
        self._degree = len(witness_set.witness_points)

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
        
    def contains(self, other):
        """
        Return True if self contains other
        """
        if not isinstance(other, Point):
            msg = "cannot understand data type"
            raise TypeError(msg)
            
        system = self._witness_set.system
        test_run = BertiniRun(system,
                              tracktype=BertiniRun.TMEMTEST,
                              component=self,
                              start=other)
        return test_run.run()
        
        #return eq
    def equals(self, other):
        pass
    
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
    def system(self):
        return self._witness_set.system
    @property
    def witness_data(self):
        return self._witness_set._witness_data
    @property
    def witness_set(self):
        return self._witness_set