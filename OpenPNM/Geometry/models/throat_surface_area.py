r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp


def cylinder(geometry, throat_diameter='throat.diameter',
             throat_length='throat.length', **kwargs):
    r"""
    Calculate throat area for a cylindrical throat
    """
    D = geometry[throat_diameter]
    L = geometry[throat_length]
    value = _sp.constants.pi*D*L
    return value


def cuboid(geometry, throat_diameter='throat.diameter',
           throat_length='throat.length', **kwargs):
    r"""
    Calculate throat area for a cuboid throat
    """
    D = geometry[throat_diameter]
    L = geometry[throat_length]
    value = 4*D*L
    return value


def extrusion(geometry, throat_perimeter='throat.perimeter',
              throat_length='throat.length', **kwargs):
    r"""
    Calculate surface area from perimeter and length -
    perimeter calculated when throat area is calculated so must be run in
    correct order
    """
    P = geometry[throat_perimeter]
    L = geometry[throat_length]
    value = P*L
    return value
