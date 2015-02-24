from __future__ import absolute_import, print_function

def scalar_num(x):
    """
    Determine if x is a scalar type for purposes of multiplication
    """
    from numbers import Number as pynum
    from sympy import Number as spnum
    
    if isinstance(x, pynum) or isinstance(x, spnum):
        return True
    else:
        return False
    
def striplines(lines, nonempty=True):
    if nonempty:
        return [l.strip() for l in lines if l != '\n']
    else:
        return [l.strip() for l in lines]