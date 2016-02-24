r"""
===============================================================================
Submodule -- throat_surface_area
===============================================================================

"""
import scipy as _sp


def cylinder(geometry, throat_diameter='throat.diameter',
             throat_length='throat.length', **kwargs):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The object containing the geometrical properties of the throats

    throat_diameter : string
        Dictionary key to the throat diameter array.  Default is
        'throat.diameter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.
    """
    D = geometry[throat_diameter]
    L = geometry[throat_length]
    value = _sp.constants.pi*D*L
    return value


def cuboid(geometry, throat_diameter='throat.diameter',
           throat_length='throat.length', **kwargs):
    r"""
    Calculate surface area for a cuboid throat

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The object containing the geometrical properties of the throats

    throat_diameter : string
        Dictionary key to the throat diameter array.  Default is
        'throat.diameter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.
    """
    D = geometry[throat_diameter]
    L = geometry[throat_length]
    value = 4*D*L
    return value


def extrusion(geometry, throat_perimeter='throat.perimeter',
              throat_length='throat.length', **kwargs):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    geometry : OpenPNM Geometry object
        The object containing the geometrical properties of the throats

    throat_perimeter : string
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.

    """
    P = geometry[throat_perimeter]
    L = geometry[throat_length]
    value = P*L
    return value
