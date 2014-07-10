r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def cylinder(network,throats,**kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    diams = network['throat.diameter'][throats]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cuboid(network,throats,**kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    diams = network['throat.diameter'][throats]
    value = (diams)**2
    return value