r"""
"""
import numpy as _np


def cylinder(target, throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    target : GenericGeometry
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    throat_diameter : str
        Dictionary key to the throat diameter array. Default is
        'throat.diameter'.
    throat_length : str
        Dictionary key to the throat length array. Default is 'throat.length'.

    Returns
    -------
    value : ndarray
        Array containing throat surface area values.

    """
    return _np.pi * target[throat_diameter] * target[throat_length]


def cuboid(target, throat_diameter='throat.diameter',
           throat_length='throat.length'):
    r"""
    Calculate surface area for a cuboid throat

    Parameters
    ----------
    target : GenericGeometry
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    throat_diameter : str
        Dictionary key to the throat diameter array.  Default is
        'throat.diameter'.
    throat_length : str
        Dictionary key to the throat length array.  Default is 'throat.length'.

    Returns
    -------
    value : ndarray
        Array containing throat surface area values.

    """
    return 4 * target[throat_diameter] * target[throat_length]


def extrusion(target, throat_perimeter='throat.perimeter',
              throat_length='throat.length'):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    target : GenericGeometry
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    throat_perimeter : str
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.
    throat_length : str
        Dictionary key to the throat length array.  Default is 'throat.length'.

    Returns
    -------
    value : ndarray
        Array containing throat surface area values.

    """
    return target[throat_perimeter] * target[throat_length]


def rectangle(target, throat_length='throat.length'):
    r"""
    Calculate surface area for a rectangular throat

    Only suitable for true 2D simulations

    Parameters
    ----------
    target : GenericGeometry
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    throat_length : str
        Dictionary key to the throat length array.  Default is 'throat.length'.

    """
    return 2 * target[throat_length]
