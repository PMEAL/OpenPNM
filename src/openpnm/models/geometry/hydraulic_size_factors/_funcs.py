import numpy as _np
import openpnm.models.geometry.conduit_lengths as _conduit_lengths
from openpnm.models.geometry import _geodocs


__all__ = [
    "spheres_and_cylinders",
    "circles_and_rectangles",
    "cones_and_cylinders",
    "intersecting_cones",
    "hybrid_cones_and_cylinders",
    "trapezoids_and_rectangles",
    "intersecting_trapezoids",
    "hybrid_trapezoids_and_rectangles",
    "pyramids_and_cuboids",
    "intersecting_pyramids",
    "hybrid_pyramids_and_cuboids",
    "cubes_and_cuboids",
    "squares_and_rectangles",
    "ncylinders_in_series"
]


@_geodocs
def spheres_and_cylinders(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    spheres and throats are cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------
    size_factors : ndarray
        Array (Nt by 3) containing conduit values for each element
        of the pore-throat-pore conduits. The array is formatted as
        ``[pore1, throat, pore2]``.

    Notes
    -----
    The hydraulic size factor is the geometrical part of the pre-factor in
    Stoke's flow:

    .. math::

        Q = \frac{A^2}{8 \pi \mu L} \Delta P
          = \frac{S_{hydraulic}}{\mu} \Delta P

    Thus :math:`S_{hydraulic}` represents the combined effect of the area
    and length of the *conduit*, which consists of a throat and 1/2 of the
    pores on each end.

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    a = 4 / (D1**3 * _np.pi**2)
    b = 2 * D1 * L1 / (D1**2 - 4 * L1**2) + _np.arctanh(2 * L1 / D1)
    F1 = a * b
    a = 4 / (D2**3 * _np.pi**2)
    b = 2 * D2 * L2 / (D2**2 - 4 * L2**2) + _np.arctanh(2 * L2 / D2)
    F2 = a * b
    Ft = Lt / (_np.pi / 4 * Dt**2)**2

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / (2 * _np.pi)

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = 1 / (16 * _np.pi**2 * It * Ft)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def circles_and_rectangles(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    circles and throats are rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.circles_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / (D1**2 * _np.sqrt(D1**2 - 4 * L1**2))
    F2 = L2 / (D2**2 * _np.sqrt(D2**2 - 4 * L2**2))
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return _np.vstack([S1, St, S2]).T


@_geodocs
def cones_and_cylinders(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated cones and throats are cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 16 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3 * _np.pi**2))
    F2 = 16 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3 * _np.pi**2))
    Ft = Lt / (_np.pi * Dt**2 / 4) ** 2

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / (2 * _np.pi)

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = 1 / (16 * _np.pi**2 * It * Ft)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def intersecting_cones(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors of intersecting cones.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.intersecting_cones(
        network=network,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 16 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3 * _np.pi**2))
    F2 = 16 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3 * _np.pi**2))

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = 1 / (2 * _np.pi)

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)
    St = _np.full(len(Lt), _np.inf)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def hybrid_cones_and_cylinders(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated cones and throats are cylinders.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.hybrid_cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 16 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3 * _np.pi**2))
    F2 = 16 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3 * _np.pi**2))
    Ft = Lt / (_np.pi * Dt**2 / 4) ** 2

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / (2 * _np.pi)

    mask = Lt == 0.0
    if mask.any():
        inv_F_t = _np.zeros(len(Ft))
        inv_F_t[~mask] = 1/Ft[~mask]
        inv_F_t[mask] = _np.inf
    else:
        inv_F_t = 1/Ft
    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = inv_F_t / (16 * _np.pi**2 * It)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def trapezoids_and_rectangles(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    trapezoids and throats are rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / 2 * (D1 + Dt) / (D1 * Dt)**2
    F2 = L2 / 2 * (D2 + Dt) / (D2 * Dt)**2
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return _np.vstack([S1, St, S2]).T


@_geodocs
def intersecting_trapezoids(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors for conduits of intersecting
    trapezoids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.intersecting_trapezoids(
        network=network,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / 2 * (D1 + Dt) / (D1 * Dt)**2
    F2 = L2 / 2 * (D2 + Dt) / (D2 * Dt)**2

    # S is 1 / (12 * F)
    S1 = 1 / (F1 * 12)
    S2 = 1 / (F2 * 12)
    St = _np.full(len(Lt), _np.inf)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def hybrid_trapezoids_and_rectangles(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    trapezoids and throats are rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.hybrid_trapezoids_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / 2 * (D1 + Dt) / (D1 * Dt)**2
    F2 = L2 / 2 * (D2 + Dt) / (D2 * Dt)**2
    Ft = Lt / Dt**3

    mask = Lt == 0.0
    if mask.any():
        inv_F_t = _np.zeros(len(Ft))
        inv_F_t[~mask] = 1/Ft[~mask]
        inv_F_t[mask] = _np.inf
    else:
        inv_F_t = 1/Ft

    # S is 1 / (12 * F)
    S1 = 1 / (F1 * 12)
    St = inv_F_t / 12
    S2 = 1 / (F2 * 12)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def pyramids_and_cuboids(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated pyramids and throats are cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.pyramids_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 1 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3))
    F2 = 1 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3))
    Ft = Lt / Dt**4

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / 6

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = 1 / (16 * _np.pi**2 * It * Ft)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def intersecting_pyramids(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors of intersecting pyramids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.intersecting_pyramids(
        network=network,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 1 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3))
    F2 = 1 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3))

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = 1 / 6

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)
    St = _np.full(len(Lt), _np.inf)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def hybrid_pyramids_and_cuboids(
    network,
    pore_diameter="pore.diameter",
    throat_coords="throat.coords"
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated pyramids and throats are cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.hybrid_pyramids_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = 1 / 3 * (L1 * (D1**2 + D1 * Dt + Dt**2) / (D1**3 * Dt**3))
    F2 = 1 / 3 * (L2 * (D2**2 + D2 * Dt + Dt**2) / (D2**3 * Dt**3))
    Ft = Lt / Dt**4

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / 6

    mask = Lt == 0.0
    if mask.any():
        inv_F_t = _np.zeros(len(Ft))
        inv_F_t[~mask] = 1/Ft[~mask]
        inv_F_t[mask] = _np.inf
    else:
        inv_F_t = 1/Ft

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = inv_F_t / (16 * _np.pi**2 * It)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def cubes_and_cuboids(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1, 1],
    throat_aspect=[1, 1, 1],
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are cubes
    and throats are cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.cubes_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^2) dx, x = [0, Li]
    F1 = L1 / D1**4
    F2 = L2 / D2**4
    Ft = Lt / Dt**4

    # I is the integral of (y^2 + z^2) dA, divided by A^2
    I1 = I2 = It = 1 / 6

    # S is 1 / (16 * pi^2 * I * F)
    S1 = 1 / (16 * _np.pi**2 * I1 * F1)
    St = 1 / (16 * _np.pi**2 * It * Ft)
    S2 = 1 / (16 * _np.pi**2 * I2 * F2)

    return _np.vstack([S1, St, S2]).T


@_geodocs
def squares_and_rectangles(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1],
    throat_aspect=[1, 1],
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    squares and throats are rectangles.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    L1, Lt, L2 = _conduit_lengths.squares_and_rectangles(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / D1**3
    F2 = L2 / D2**3
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return _np.vstack([S1, St, S2]).T


@_geodocs
def ncylinders_in_series(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    n=5,
):
    r"""
    Computes hydraulic size factors for conduits of spheres and cylinders
    with the spheres approximated as N cylinders in series.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s
    n : int
        Number of cylindrical divisions for each pore

    Returns
    -------

    Notes
    -----

    """
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split('.', 1)[-1]).T
    # Ensure throats are never bigger than connected pores
    Dt = _np.minimum(Dt, 0.99 * _np.minimum(D1, D2))
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T
    dL1 = _np.linspace(0, L1, num=n)
    dL2 = _np.linspace(0, L2, num=n)
    r1 = D1 / 2 * _np.sin(_np.arccos(dL1 / (D1 / 2)))
    r2 = D2 / 2 * _np.sin(_np.arccos(dL2 / (D2 / 2)))

    gtemp = (_np.pi * r1 ** 4 / (8 * L1 / n)).T
    F1 = 1 / _np.sum(1 / gtemp, axis=1)
    gtemp = (_np.pi * r2 ** 4 / (8 * L2 / n)).T
    F2 = 1 / _np.sum(1 / gtemp, axis=1)
    Ft = (_np.pi * (Dt / 2) ** 4 / (8 * Lt)).T

    return _np.vstack([F1, Ft, F2]).T
