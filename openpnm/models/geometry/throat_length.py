r"""
"""
import numpy as _np
from numpy.linalg import norm as _norm
from openpnm.utils import logging as _logging
_logger = _logging.getLogger(__name__)


def ctc(target):
    r"""
    Calculates throat length assuming point-like pores, i.e. center-to-center
    distance between pores. Also, this model assumes that pores and throat
    centroids are colinear.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    Returns
    -------
    value : ndarray
        Array containing throat length values.

    """
    network = target.project.network
    throats = network.throats(target.name)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    value = _norm(C1 - C2, axis=1)
    return value


def classic(target, pore_diameter='pore.diameter'):
    r"""
    Finds throat length as the pore-to-pore center distance, less the radii of
    each pore.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    throats = network.throats(target.name)
    cn = network['throat.conns'][throats]
    ctc_dist = ctc(target)
    value = ctc_dist - network[pore_diameter][cn].sum(axis=1) / 2
    return value


def spheres_and_cylinders(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are spheres and throats are
    cylinders.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.spheres_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def circles_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are circles and throats are
    rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.circles_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def cones_and_cylinders(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cones and throats are
    cylinders.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cones_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def trapezoids_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are trapezoids and throats are
    rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.trapezoids_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def pyramids_and_cuboids(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are pyramids and throats are
    cuboids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.pyramids_and_cuboids(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def cubes_and_cuboids(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are cubes and throats are
    cuboids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.cubes_and_cuboids(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]


def squares_and_rectangles(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter'
):
    r"""
    Finds throat length assuming pores are squares and throats are
    rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : str
        Dictionary key of the pore diameter values.

    throat_diameter : str
        Dictionary key of the throat diameter values.

    Returns
    -------
    ndarray
        Array containing throat length values.

    """
    from openpnm.models.geometry import conduit_lengths
    out = conduit_lengths.squares_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )
    return out[:, 1]
