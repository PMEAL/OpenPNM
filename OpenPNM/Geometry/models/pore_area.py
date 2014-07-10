r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def spherical(network,
              pores,
              pore_diameter='pore.diameter',
              **kwargs):
    r"""
    Calculate cross-sectional area for a spherical pore
    """
    diams = network[pore_diameter][pores]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cubic(network,
          pore,
          pore_diameter='pore.diameter',
          **kwargs):
    r"""
    Calculate cross-sectional pore area for a cubic pore
    """
    diams = network[pore_diameter][pore]
    value = (diams)**2
    return value