import numpy as _np
from openpnm.utils import Docorator


__all__ = ["cylinder",
           "cuboid",
           "extrusion",
           "rectangle"]
docstr = Docorator()


@docstr.get_sections(base='models.geometry.throat_surface_area',
                     sections=['Parameters', 'Returns'])
@docstr.dedent
def cylinder(network, throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    thorat_area : str
        Name of the dictionary key on ``network`` where the array containing
        throat area values is stored

    Returns
    -------
    surface_areas : ndarray
        A numpy ndarray containing throat surface area values

    """
    return _np.pi * network[throat_diameter] * network[throat_length]


@docstr.dedent
def cuboid(network, throat_diameter='throat.diameter',
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
    return 4 * network[throat_diameter] * network[throat_length]


def extrusion(network, throat_perimeter='throat.perimeter',
              throat_length='throat.length'):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tlen)s
    throat_perimeter : str
        Dictionary key to the throat perimeter array.  Default is
        'throat.perimeter'.

    Returns
    -------
    %(models.geometry.throat_surface_area.returns)s

    """
    return network[throat_perimeter] * network[throat_length]


def rectangle(network, throat_length='throat.length'):
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
    return 2 * network[throat_length]
