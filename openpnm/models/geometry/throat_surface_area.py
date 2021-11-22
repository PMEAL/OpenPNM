import numpy as _np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_surface_area',
                     sections=['Parameters', 'Returns'])
def cylinder(target, throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    thorat_area : str
        Name of the dictionary key on ``target`` where the array containing
        throat area values is stored
    throat_length : str
        Name of the dictionary key on ``target`` where the array containing
        throat length values is stored

    Returns
    -------
    areas : ndarray
        A numpy ndarray containing throat surface area values

    """
    return _np.pi * target[throat_diameter] * target[throat_length]


@docstr.dedent
def cuboid(target, throat_diameter='throat.diameter',
           throat_length='throat.length'):
    r"""
    Calculate surface area for a cuboid throat

    Parameters
    ----------
    %(models.geometry.throat_surface_area.parameters)s

    Returns
    -------
    %(models.geometry.throat_surface_area.returns)s

    """
    return 4 * target[throat_diameter] * target[throat_length]


def extrusion(target, throat_perimeter='throat.perimeter',
              throat_length='throat.length'):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.

    throat_perimeter : string
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat surface area values.

    """
    return target[throat_perimeter] * target[throat_length]


def rectangle(target, throat_length='throat.length'):
    r"""
    Calculate surface area for a rectangular throat

    Only suitable for true 2D simulations

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.

    throat_length : string
        Dictionary key to the throat length array.  Default is 'throat.length'.

    """
    return 2 * target[throat_length]
