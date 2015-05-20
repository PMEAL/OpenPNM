r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp


def sphere(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate internal surface area for a spherical pore
    """
    R = geometry[pore_diameter]/2
    value = 4*_sp.constants.pi*R**2
    return value


def cube(geometry, pore_diameter='pore.diameter', **kwargs):
    r"""
    Calculate internal surface area for a cubic pore
    """
    D = geometry[pore_diameter]
    value = 6*D**2
    return value
