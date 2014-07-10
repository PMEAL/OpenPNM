r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp

def cylinder(network,
             throats,
             throat_diameter='throat.diameter',
             throat_length='throat.length',
             **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = network[throat_diameter][throats]
    L = network[throat_length][throats]
    value = _sp.constants.pi/(D)*L
    return value

def cuboid(network,
           throats,
           throat_diameter='throat.diameter',
           throat_length='throat.length',
           **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = network[throat_diameter][throats]
    L = network[throat_length][throats]
    value = 4*D*L
    return value

def extrusion(network,
              throats,
              throat_perimeter='throat.perimeter',
              throat_length='throat.length',
              **kwargs):
    r"""
    Calculate surface area from perimeter and length - 
    perimeter calculated when throat area is calculated so must be run in correct order
    """
    P = network[throat_perimeter][throats]
    L = network[throat_length][throats]
    value = P*L
    return value
    
def voronoi(network,
            throats,
            throat_perimeter='throat.perimeter',
            throat_length='throat.length',
            **kwargs):
    r"""
    Calculate surface area from perimeter and length - 
    perimeter calculated when throat area is calculated so must be run in correct order
    """
    P = network[throat_perimeter][throats]
    L = network[throat_length][throats]
    value = P*L
    return value