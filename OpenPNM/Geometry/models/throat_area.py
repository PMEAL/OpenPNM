r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def cylinder(network,
             throats,
             throat_diameter = 'throat.diameter',
             **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    diams = network[throat_diameter][throats]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cuboid(network,
           throats,
           throat_diameter = 'throat.diameter',
           **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    diams = network[throat_diameter][throats]
    value = (diams)**2
    return value