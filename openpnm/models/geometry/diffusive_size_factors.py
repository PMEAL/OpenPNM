r"""
The diffusive size factor is the geometrical part of the pre-factor in
Fick's law:

.. math::

    n_A = \frac{A}{L} \Delta C_A
        = S_{diffusive} D_{AB} \Delta C_A

Thus :math:`S_{diffusive}` represents the combined effect of the area and
length of the *conduit*, which consists of a throat and 1/2 of the pore
on each end.

"""

import numpy as _np
import openpnm.models.geometry.conduit_lengths as _conduit_lengths
import openpnm.geometry.GenericGeometry as _GenericGeometry
from .misc import _get_conduit_diameters

__all__ = [
    "spheres_and_cylinders",
    "circles_and_rectangles",
    "cones_and_cylinders",
    "trapezoids_and_rectangles",
    "pyramids_and_cuboids",
    "cubes_and_cuboids",
    "squares_and_rectangles",
    "intersecting_cones",
    "intersecting_trapezoids",
    "intersecting_pyramids",
    "ncylinders_in_series"
]


def spheres_and_cylinders(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive shape coefficient for conduits assuming pores are
    spheres and throats are cylinders.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 2 / (D1 * _np.pi) * _np.arctanh(2 * L1 / D1)
    F2 = 2 / (D2 * _np.pi) * _np.arctanh(2 * L2 / D2)
    Ft = Lt / (_np.pi / 4 * Dt ** 2)

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def circles_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive shape coefficient for conduits assuming pores are
    circles and throats are rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.circles_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 0.5 * _np.arcsin(2 * L1 / D1)
    F2 = 0.5 * _np.arcsin(2 * L2 / D2)
    Ft = Lt / Dt

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def cones_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive shape coefficient assuming pores are truncated cones
    and throats are cylinders.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : string
        Dictionary key of the pore diameter values.
    throat_diameter : string
        Dictionary key of the throat diameter values.

    Returns
    -------
    A dictionary containing the diffusive_shape_coefficient, which can be
    accessed via the dict keys 'pore1', 'pore2', and 'throat'.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.cones_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = 4 * L1 / (D1 * Dt * _np.pi)
    F2 = 4 * L2 / (D2 * Dt * _np.pi)
    Ft = Lt / (_np.pi * Dt**2 / 4)

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def trapezoids_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Compute diffusive shape coefficient for conduits assuming pores are
    trapezoids and throats are rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.trapezoids_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = L1 * _np.log(Dt / D1) / (Dt - D1)
    F2 = L2 * _np.log(Dt / D2) / (Dt - D2)
    Ft = Lt / Dt

    # Edge case where Di = Dt
    mask = _np.isclose(D1, Dt)
    F1[mask] = (L1 / D1)[mask]
    mask = _np.isclose(D2, Dt)
    F2[mask] = (L2 / D2)[mask]

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def pyramids_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive size factor for conduits of truncated pyramids and
    cuboids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.pyramids_and_cuboids(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = L1 / (D1 * Dt)
    F2 = L2 / (D2 * Dt)
    Ft = Lt / Dt**2

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def cubes_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1, 1],
    throat_aspect=[1, 1, 1],
):
    r"""
    Computes diffusive shape coefficient for conduits assuming pores are
    cubes and throats are cuboids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.cubes_and_cuboids(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = L1 / D1**2
    F2 = L2 / D2**2
    Ft = Lt / Dt**2

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def squares_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1],
    throat_aspect=[1, 1],
):
    r"""
    Computes diffusive size factor for conduits assuming pores are squares
    and throats are rectangles.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    L1, Lt, L2 = _conduit_lengths.squares_and_rectangles(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A) dx, x = [0, Li]
    F1 = L1 / D1
    F2 = L2 / D2
    Ft = Lt / Dt

    return {"pore1": 1 / F1, "throat": 1 / Ft, "pore2": 1 / F2}


def intersecting_cones(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None
):
    r"""
    Computes diffusive size factor for conduits of intersecting cones.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    raise NotImplementedError


def intersecting_trapezoids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes diffusive shape coefficient for conduits of intersecting
    trapezoids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    raise NotImplementedError


def intersecting_pyramids(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes diffusive size factor for conduits of intersecting pyramids.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    raise NotImplementedError


def ncylinders_in_series(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    n=5
):
    r"""
    Computes diffusive size factors for conduits of spheres and cylinders
    with the spheres approximated as N cylinders in series.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Dictionary key pointing to the pore diameter values.
    throat_diameter : str
        Dictionary key pointing to the throat diameter values.
    n : int
        Number of cylindrical divisions for each pore

    Notes
    -----
    The diffusive size factor is the geometrical part of the pre-factor in
    Fick's law:

    .. math::

        n_A = \frac{A}{L} \Delta C_A
            = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)
    # Ensure throats are never bigger than connected pores
    Dt = _np.minimum(Dt, 0.99 * _np.minimum(D1, D2))
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    ).T
    dL1 = _np.linspace(0, L1, num=n)
    dL2 = _np.linspace(0, L2, num=n)
    r1 = D1 / 2 * _np.sin(_np.arccos(dL1 / (D1 / 2)))
    r2 = D2 / 2 * _np.sin(_np.arccos(dL2 / (D2 / 2)))

    gtemp = (_np.pi * r1**2 / (L1 / n)).T
    F1 = 1 / _np.sum(1 / gtemp, axis=1)
    gtemp = (_np.pi * r2 ** 2 / (L2 / n)).T
    F2 = 1 / _np.sum(1 / gtemp, axis=1)
    Ft = (_np.pi * (Dt / 2)**2 / (Lt)).T

    return {"pore1": F1, "throat": Ft, "pore2": F2}
