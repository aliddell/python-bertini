from __future__ import print_function

from fractions import Fraction as fraction
from os import chdir
from os.path import isfile
import re
from tempfile import mkdtemp

from sympy import sympify, Matrix as spmatrix

from naglib import TEMPDIR as basedir
from naglib.core import IrreducibleComponent, AffinePoint, ProjectivePoint, WitnessPoint, WitnessSet
from naglib.core.algebra import LinearSlice
from naglib.core.misc import striplines
from naglib.bertini.sysutils import call_bertini
from naglib.bertini.fileutils import write_input, parse_witness_data

def get_components(dirname, system):
    #variables = system.variables
    #sys_homvar = system.homvar
    if system.homvar:
        proj_dim = len(system.variables)
        #homsys = system
    else:
        proj_dim = len(system.variables) + 1
        #homsys = system.homogenize()
    homsys  = system.homogenize()
    homvars = homsys.variables
    homvar  = homsys.homvar
    
    witness_file = dirname + '/witness_data'
    witness_data,wdinfo = parse_witness_data(witness_file)
    
    components = []
    
    for c in witness_data:
        codim       = c['codim']
        homVarConst = c['homVarConst']
        points      = c['points']
        coeffs      = c['slice']
        
        comp_isprojective = homVarConst == 0
        if coeffs is not None:
            if comp_isprojective:
                slice = LinearSlice(coeffs, homvars, homvar)
            if not system.homvar:
                slice = slice.dehomogenize()
        else:
            slice = None
        if comp_isprojective and system.homvar:
            dim = proj_dim - codim
            #if slice is not None:
                #slice = sympify(slice.tolist())
        elif comp_isprojective: # and not sys_homvar
            dim = proj_dim - codim - 1
            #if slice is not None:
                #homcol = slice.column(0)
                #B1 = slice.tolist()
                #for row in range(len(B1)):
                    #B1[row] = spmatrix(B1[row][1:])
                    #B1[row] = list(B1[row]/homcol[row])
                #slice = spmatrix(sympify(B1))
            
        dim_list = {}
        
        for point in points:
            comp_id = point['component number']
            if comp_isprojective:
                coord = ProjectivePoint(point['coordinates'])
                if not system.homvar:
                    coord = coord.dehomogenize()
            else:
                coord = AffinePoint(point['coordinates'])
                
            pt = WitnessPoint(coord, comp_id)
            
            if not dim_list.has_key(comp_id):
                dim_list[comp_id] = []
                
            dim_list[comp_id].append(pt)
        
        for comp_id in dim_list.keys():
            ws = WitnessSet(system.copy(), slice, dim_list[comp_id])
            component = IrreducibleComponent(ws, dim, comp_id, wdinfo)
            components.append(component)
            
    return components

def compute_NID(system):
    """
    Compute the numerical irreducible decomposition of
    PolynomialSystem system
    
    Returns an iterable of IrreducibleComponent
    """
    dirname = mkdtemp(prefix=basedir)
    input_file = dirname + '/input'
    config = {'filename':input_file, 'TrackType':1}

    write_input(system, config)
    call_bertini(input_file)

    components = get_components(dirname, system)

    return components

