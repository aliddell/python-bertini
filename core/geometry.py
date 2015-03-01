from __future__ import print_function

from naglib.startup import TEMPDIR as basedir
from naglib.core.base import NAGobject
from naglib.core.witnessdata import WitnessPoint, WitnessSet

class IrreducibleComponent(NAGobject):
    """
    An irreducible component of an algebraic set
    """
    def __init__(self, witness_set, dim, component_id, wdinfo):
        """
        Initialize the IrreducibleComponent object.
        
        Keyword arguments:
        witness_set  -- WitnessSet, the witness set describing this component
        dim          -- int, the dimension of the component
        component_id -- int, the component number for this dimension
        wdinfo       -- list, the lines of the witness_data file
        """
        self._witness_set = witness_set
        self._dim = dim
        self._component_id = component_id
        self._degree = len(witness_set.witness_points)
        self._wdinfo = wdinfo

    def __str__(self):
        """
        x.__str__() <==> str(x)
        """
        dim = self._dim
        cid = self._component_id
        deg = self._degree
        return '{0}-dimensional irreducible component ({2}) of degree {1}'.format(dim,deg,cid)

    def __repr__(self):
        """
        x.__repr__() <==> repr(x)
        """
        dim = self._dim
        cid = self._component_id
        deg = self._degree
        wst = self._witness_set
        repstr = 'IrreducibleComponent({0},{1},{2},{3})'.format(dim,cid,deg,repr(wst))
        return repstr
    
    #def __eq__(self, other):
        #"""
        #x.__eq__(y) <==> x == y
        #"""
        #if not isinstance(other, IrreducibleComponent):
            #return False
        
        #ssys = self._system
        #sdim = self._dim
        #scid = self._component_id
        #osys = other.system
        #odim = other.dim
        #ocid = other.component_id
        #eq = (ssys == osys and sdim == odim and scid == ocid)
        
        #return eq
    def equals(self, other, strict=True):
        pass
    
    def sample(self, numpoints=1, usebertini=True):
        """
        Sample a point from self
        """
        dim     = self._dim
        comp_id = self._component_id
        system  = self.witness_set.system
        wdinfo  = self._wdinfo
        
        points = None
        if usebertini:
            from tempfile import mkdtemp
            from naglib.bertini.fileutils import write_input, read_points
            from naglib.bertini.sysutils import call_bertini as call
            
            dirname = mkdtemp(prefix=basedir)
            wdfile = dirname + '/witness_data'
            inputfile = dirname + '/input'
            instructions = dirname + '/instructions'
            sampled = dirname + '/sampled'
            
            # write out witness_data file
            fh = open(wdfile, 'w')
            for line in wdinfo:
                print(line, file=fh)
            fh.close()
            
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
            msg = "nothing to use yet but Bertini"
            raise NotImplementedError(msg)
        
        return points
    
    @property
    def system(self):
        return self._witness_set.system
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