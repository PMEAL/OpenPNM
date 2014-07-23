r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def cylinder(geometry,
             throat_diameter='throat.diameter',
             **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    diams = geometry[throat_diameter]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cuboid(geometry,
           throat_diameter='throat.diameter',
           **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    diams = geometry[throat_diameter]
    value = (diams)**2
    return value