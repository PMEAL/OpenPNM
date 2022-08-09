from numpy import pi as _pi
from openpnm.utils import Docorator

__all__ = [
    "cone",
    "sphere",
    "cube",
    "circle",
    "square"
]

docstr = Docorator()


@docstr.get_sections(base='models.geometry.pore_cross_section',
                     sections=['Parameters', 'Returns'])
@docstr.dedent
def sphere(network, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.pdia)s


    Returns
    -------
    areas : ndarray
        A numpy ndarry containing pore cross-sectional area values

    """
    D = network[pore_diameter]
    return _pi/4 * D**2


def cone(network, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cone

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    """
    D = network[pore_diameter]
    return _pi / 4 * D**2


def cube(network, pore_diameter='pore.diameter'):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    %(models.geometry.pore_cross_section.parameters)s

    Returns
    -------
    %(models.geometry.pore_cross_section.returns)s

    """
    D = network[pore_diameter]
    return D**2


def circle(network, pore_diameter='pore.diameter'):
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
    return network[pore_diameter]


def square(network, pore_diameter='pore.diameter'):
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
    return network[pore_diameter]
