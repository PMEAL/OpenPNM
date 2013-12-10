
"""
module ThroatDiameter
===============================================================================

"""
import scipy as sp

def constant(value,**params):
    r"""
    Assigns specified constant value
    """
    print('constant')
    return value

def cylinder(**params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    print('cylinder: nothing yet')

def cuboid(**params):
    r"""
    Calculate throat diameter from seeds for a cuboidal throat
    """
    print('cuboid: nothing yet')
    
def voronoi(**params):
    r"""
    Calculate throat diameter from analysis of Voronoi facets
    """
    print('voronoi: nothing yet')