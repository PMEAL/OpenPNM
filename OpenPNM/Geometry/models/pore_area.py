r"""
===============================================================================
pore_area -- Models for cross-sectional area of a pore body
===============================================================================

"""
import scipy as _sp

def spherical(geometry,
              pore_diameter='pore.diameter',
              **kwargs):
    r"""
    Calculate cross-sectional area for a spherical pore
    """
    diams = geometry[pore_diameter]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cubic(geometry,
          pore_diameter='pore.diameter',
          **kwargs):
    r"""
    Calculate cross-sectional pore area for a cubic pore
    """
    diams = geometry[pore_diameter]
    value = (diams)**2
    return value