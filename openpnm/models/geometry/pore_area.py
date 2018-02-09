r"""
===============================================================================
pore_area -- Models for cross-sectional area of a pore body
===============================================================================

"""
import scipy as _sp


def spherical(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.

    """
    diams = target[pore_diameter]
    value = _sp.pi/4*(diams)**2
    return value


def cubic(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    target : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the array on the Geometry object containing the
        pore diameter values necessary to find the area.


    """
    diams = target[pore_diameter]
    value = diams**2
    return value
