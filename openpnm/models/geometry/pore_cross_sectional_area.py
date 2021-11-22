from numpy import pi as _pi
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.get_sections(base='models.geometry.pore_cross_section',
                     sections=['Parameters', 'Returns'])
def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere.

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Name of the dictionary key on ``target`` where the array containing
        pore diameter values is stored

    Returns
    -------
    areas : ndarray
        A numpy ndarry containing pore cross-sectional area values

    """
    D = target[pore_diameter]
    return _pi/4 * D**2


def cone(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cone.

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    """
    D = target[pore_diameter]
    return _pi / 4 * D**2


def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    """
    D = target[pore_diameter]
    return D**2


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a circle

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return target[pore_diameter]


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a square

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return target[pore_diameter]
