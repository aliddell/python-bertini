from __future__ import print_function

from fractions import Fraction as fraction
from os import chdir
from os.path import isfile
import re
from tempfile import mkdtemp

from mpmath import matrix as mpmatrix
from sympy import sympify, Matrix as spmatrix

from naglib.core.datatypes import LinearSystem, IrreducibleComponent, WitnessPoint, WitnessSet
from naglib.bertini.sysutils import call_bertini
from naglib.misc import striplines
from naglib.bertini.fileutils import write_input, parse_witness_data

def compute_NID(system):
    """
    Compute the numerical irreducible decomposition of
    PolynomialSystem system
    
    Returns an iterable of IrreducibleComponent
    """
    dirname = mkdtemp()
    input_file = dirname + '/input'
    config = {'filename':input_file, 'TrackType':1}

    write_input(system, config)
    call_bertini(input_file)

    components = __get_components(dirname, system)

    return components

def __get_components(dirname, system):
    variables = system.variables
    sys_isprojective = system.isprojective
    if sys_isprojective:
        proj_dim = len(variables)
    else:
        proj_dim = len(variables) + 1
    
    witness_file = dirname + '/witness_data'
    witness_data = parse_witness_data(witness_file)
    
    components = []
    
    for c in witness_data:
        codim       = c['codim']
        homVarConst = c['homVarConst']
        points      = c['points']
        B           = c['B']
        
        comp_isprojective = homVarConst == 0
        if comp_isprojective and sys_isprojective:
            dim = proj_dim - codim
            if B is not None:
                B = sympify(B.tolist())
        elif comp_isprojective: # and not sys_isprojective
            dim = proj_dim - codim - 1
            if B is not None:
                homcol = B.column(0)
                B1 = B.tolist()
                for row in range(len(B1)):
                    B1[row] = mpmatrix(B1[row][1:])
                    B1[row] = list(B1[row]/homcol[row])
                B = spmatrix(sympify(B1))
            
        dim_list = {}
        
        for point in points:
            comp_id = point['component number']
            coord = point['coordinates']
            if comp_isprojective and not sys_isprojective:
                # dehomogenize
                coord = coord[1:]/coord[0]
            pt = WitnessPoint(dim, comp_id, coord, comp_isprojective)
            
            if not dim_list.has_key(comp_id):
                dim_list[comp_id] = []
                
            dim_list[comp_id].append(pt)
        
        if B is not None:
            slice = LinearSystem(B*variables, variables)
        else:
            slice = None
        for comp in dim_list.keys():
            ws = WitnessSet(system, slice, dim_list[comp], comp_isprojective)
            component = IrreducibleComponent(system, dim, comp, ws, dirname, comp_isprojective)
            components.append(component)
            
        return components