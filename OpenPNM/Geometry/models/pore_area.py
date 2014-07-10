r"""
===============================================================================
Submodule -- throat_area
===============================================================================

"""
import scipy as _sp

def spherical(network,pores,**kwargs):
    r"""
    Calculate cross-sectional area for a spherical pore
    """
    diams = network['pore.diameter'][pores]
    value = _sp.constants.pi/4*(diams)**2
    return value

def cubic(network,pore,**kwargs):
    r"""
    Calculate cross-sectional pore area for a cubic pore
    """
    diams = network['pore.diameter'][pore]
    value = (diams)**2
    return value