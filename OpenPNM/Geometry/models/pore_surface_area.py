r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp

def sphere(geometry,
             pore_diameter='pore.diameter',
             **kwargs):
    r"""
    Calculate internal surface area for a spherical pore
    """
    D = geometry[pore_diameter]
    value = 2*_sp.constants.pi*D
    return value

def cube(geometry,
           pore_diameter='pore.diameter',
           **kwargs):
    r"""
    Calculate internal surface area for a cubic pore
    """
    D = geometry[pore_diameter]
    value = 6*D**2
    return value

