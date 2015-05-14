from __future__ import print_function

from naglib.bertini.sysutils import BertiniRun
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
        if type(other) not in (tuple, list, AffinePoint, ProjectivePoint):
            msg = "cannot understand data type"
            raise TypeError(msg)
        elif isinstance(other, Point):
            other = [other]
        
        numpoints = len(other)
        from tempfile import mkdtemp
        from naglib.bertini.fileutils import striplines, write_input, fprint
        from naglib.bertini.sysutils import call_bertini
        system = self._witness_set.system
        dirname = mkdtemp(prefix=basedir)
        inputf = dirname + '/input'
        fprint(other, dirname + '/member_points')
        config = {'filename':inputf, 'TrackType':3}
        write_input(system, config)
        self.write_witness(dirname)
        call_bertini(inputf)
        fh = open(dirname + '/output_membership', 'r')
        lines = striplines(fh.readlines())
        indices = [lines.index('Testing {0}'.format(i)) for i in range(numpoints)] + [len(lines)]
        truevals = [False for i in range(numpoints)]
        for i in range(numpoints):
            if indices[i+1] - indices[i] > 4:
               truevals[i] = True
        
        return truevals
        
        #return eq
    def equals(self, other):
        pass
    
    def sample(self, numpoints=1, usebertini=True):
        """
        Sample a point from self
        """
        dim     = self.dim
        comp_id = self._component_id
        system  = self.witness_set.system
        
        points = None
        if usebertini:
            sample_run = BertiniRun(system, BertiniRun.TSAMPLE, sample=numpoints, components=[self])
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