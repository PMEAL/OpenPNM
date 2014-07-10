r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp

def cylinder(network,throats,**kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network['throat.diameter'][throats]
    L = network['throat.length'][throats]
    value = _sp.constants.pi/(D)*L
    return value

def cuboid(network,throats,**kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network['throat.diameter'][throats]
    L = network['throat.length'][throats]
    value = 4*D*L
    return value
    
def voronoi(network,throats,**kwargs):
    r"""
    Calculate surface area from perimeter and lenght - 
    perimeter calculated when throat area is calculated so must be run in correct order
    """
    P = network['throat.perimeter'][throats]
    L = network['throat.length'][throats]
    value = P*L
    return value