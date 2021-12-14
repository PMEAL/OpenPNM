import numpy as _np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_surface_area',
                     sections=['Parameters', 'Returns'])
@docstr.dedent
def cylinder(target, throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    thorat_area : str
        Name of the dictionary key on ``target`` where the array containing
        throat area values is stored

    Returns
    -------
    surface_areas : ndarray
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
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    throat_perimeter : string
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.

    Returns
    -------
    %(models.geometry.throat_surface_area.returns)s

    """
    return target[throat_perimeter] * target[throat_length]


def rectangle(target, throat_length='throat.length'):
    r"""
    Calculate surface area for a rectangular throat

    Only suitable for true 2D simulations

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s

    Returns
    -------
    %(models.geometry.throat_surface_area.returns)s

    """
    return 2 * target[throat_length]
