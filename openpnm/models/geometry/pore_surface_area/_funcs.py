import numpy as _np
from openpnm.utils import Docorator

__all__ = ["sphere",
           "circle",
           "cube",
           "square"
           ]
docstr = Docorator()


@docstr.get_sections(base='models.geometry.pore_surface_area',
                     sections=['Parameters', 'Returns', 'Notes'])
@docstr.dedent
def sphere(
    target,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area'
):
    r"""
    Calculates internal surface area of pore bodies assuming they are spheres,
    then subtracts the areas of the neighboring throats

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.pdia)s
    throat_cross_sectional_area : str
        Name of the dictionary key on ``target`` where the array containing
        throat cross sectional area values is stored

    Returns
    -------
    surface_areas : ndarray
        Numpy ndarry containing pore surface area values

    Notes
    -----
    This function subtracts the area of the neighboring throats in a
    crude way, by simply considering the throat cross-sectional area,
    thus not accounting for the actual curvature of the intersection.

    """
    network = target.project.network
    R = target[pore_diameter] / 2
    Asurf = 4 * _np.pi * R**2
    Tn = network.find_neighbor_throats(pores=target.Ps, flatten=False)
    Tsurf = _np.array([network[throat_cross_sectional_area][Ts].sum() for Ts in Tn])
    value = Asurf - Tsurf
    return value


@docstr.dedent
def circle(
    target,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area'
):
    r"""
    Calculates internal surface area of pore bodies assuming they are
    circular then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(models.geometry.pore_surface_area.parameters)s

    Returns
    -------
    %(models.geometry.pore_surface_area.returns)s

    Notes
    -----
    %(models.geometry.pore_surface_area.notes)s

    """
    network = target.project.network
    R = target[pore_diameter] / 2
    Asurf = 2 * _np.pi * R
    Tn = network.find_neighbor_throats(pores=target.Ps, flatten=False)
    Tsurf = _np.array([network[throat_cross_sectional_area][Ts].sum() for Ts in Tn])
    value = Asurf - Tsurf
    return value


def cube(
    target,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area'
):
    r"""
    Calculates internal surface area of pore bodies assuming they are cubes
    then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(models.geometry.pore_surface_area.parameters)s

    Returns
    -------
    %(models.geometry.pore_surface_area.returns)s

    """
    network = target.project.network
    D = target[pore_diameter]
    Tn = network.find_neighbor_throats(pores=target.Ps, flatten=False)
    Tsurf = _np.array([network[throat_cross_sectional_area][Ts].sum() for Ts in Tn])
    value = 6 * D**2 - Tsurf
    return value


def square(
    target,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area'
):
    r"""
    Calculates internal surface area of pore bodies assuming they are
    squares then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(models.geometry.pore_surface_area.parameters)s

    Returns
    -------
    %(models.geometry.pore_surface_area.returns)s

    """
    network = target.project.network
    D = target[pore_diameter]
    Tn = network.find_neighbor_throats(pores=target.Ps, flatten=False)
    Tsurf = _np.array([network[throat_cross_sectional_area][Ts].sum() for Ts in Tn])
    value = 4 * D - Tsurf
    return value
