import numpy as np

from openpnm.models.geometry import _geodocs

__all__ = [
    "spheres_and_cylinders",
    "circles_and_rectangles",
    "cones_and_cylinders",
    "intersecting_cones",
    "hybrid_cones_and_cylinders",
    "trapezoids_and_rectangles",
    "hybrid_trapezoids_and_rectangles",
    "intersecting_trapezoids",
    "pyramids_and_cuboids",
    "intersecting_pyramids",
    "hybrid_pyramids_and_cuboids",
    "cubes_and_cuboids",
    "squares_and_rectangles",
]


@_geodocs
def spheres_and_cylinders(
    network,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Calculates conduit lengths in the network assuming pores are spheres
    and throats are cylinders.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------
    lengths : ndarray
        Array (Nt by 3) containing conduit values for each element
        of the pore-throat-pore conduits. The array is formatted as
        ``[pore1, throat, pore2]``.

    """
    L_ctc = _get_L_ctc(network)
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split(".", 1)[-1]).T

    # Handle the case where Dt > Dp
    if (Dt > D1).any() or (Dt > D2).any():
        _raise_incompatible_data()

    # If spheres do not overlap:
    L1 = np.sqrt(D1**2 - Dt**2) / 2
    L2 = np.sqrt(D2**2 - Dt**2) / 2
    Lt = L_ctc - (L1 + L2)

    if np.any(Lt < 0):  # Find pores that touch/overlap
        # Find diameter of overlap between spheres
        d = L_ctc
        R1 = D1/2
        R2 = D2/2
        # Check distance to the intersection
        L1_int = (d**2 - R2**2 + R1**2) / (2*d)
        if np.any(L1_int < 0) or np.any(L1_int > L_ctc):
            raise Exception('The pores overlap too much')
        D_int = 2/(2*d) * np.sqrt(4*d**2 * R1**2 - (d**2 - R2**2 + R1**2)**2)
        mask = D_int > Dt
        if np.any(mask):
            d = d[mask]
            R1 = R1[mask]
            R2 = R2[mask]
            L1[mask] = (d**2 - R2**2 + R1**2) / (2*d)
            L2[mask] = L_ctc[mask] - L1[mask]
            Lt[mask] = 1e-15

    return np.vstack((L1, Lt, L2)).T


@_geodocs
def circles_and_rectangles(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are circles
    and throats are rectangles.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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
    return spheres_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter,
    )


@_geodocs
def cones_and_cylinders(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are cones
    and throats are cylinders.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    L_ctc = _get_L_ctc(network)
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split(".", 1)[-1]).T

    L1 = D1 / 2
    L2 = D2 / 2

    # Handle throats w/ overlapping pores
    _L1 = (4 * L_ctc**2 + D1**2 - D2**2) / (8 * L_ctc)
    mask = L_ctc - 0.5 * (D1 + D2) < 0
    L1[mask] = _L1[mask]
    L2[mask] = (L_ctc - L1)[mask]
    Lt = np.maximum(L_ctc - (L1 + L2), 1e-15)
    return np.vstack((L1, Lt, L2)).T


@_geodocs
def intersecting_cones(network, pore_coords="pore.coords", throat_coords="throat.coords"):
    r"""
    Calculates conduit lengths in the network assuming pores are cones
    that intersect. Therefore, the throat is the cross sectional plane where
    two pores meet and has negligible/zero volume.

    A conduit is defined as ( 1/2 pore - 1/2 pore ).

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    """
    P12 = network["throat.conns"]
    p_coords = network[pore_coords]
    t_coords = network[throat_coords]
    L1 = np.sqrt(np.sum(((p_coords[P12[:, 0]] - t_coords)) ** 2, axis=1))
    L2 = np.sqrt(np.sum(((p_coords[P12[:, 1]] - t_coords)) ** 2, axis=1))
    Lt = np.zeros(len(network.Ts))
    return np.vstack((L1, Lt, L2)).T


@_geodocs
def hybrid_cones_and_cylinders(
    network, pore_diameter="pore.diameter", throat_coords="throat.coords"
):
    r"""
    Calculates conduit lengths in the network assuming pores are cones
    and throats are cylinders.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    """
    L_ctc = _get_L_ctc(network)
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split(".", 1)[-1]).T

    L1 = D1 / 2
    L2 = D2 / 2
    Lt = np.maximum(L_ctc - (L1 + L2), 1e-15)
    # Handle intersecting pores
    XYp = network.coords[network.conns]
    XYt = network[throat_coords]
    _L1, _L2 = np.linalg.norm(XYp - XYt[:, None], axis=2).T
    mask1 = _L1 < L1
    mask2 = _L2 < L2
    mask = np.logical_or(mask1, mask2)
    if mask.any():
        L1[mask] = _L1[mask]
        L2[mask] = _L2[mask]
        Lt[mask] = 0.0
    return np.vstack((L1, Lt, L2)).T


@_geodocs
def trapezoids_and_rectangles(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are
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
    return cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter,
    )


@_geodocs
def intersecting_trapezoids(
    network, pore_coords="pore.coords", throat_coords="throat.coords"
):
    r"""
    Calculates conduit lengths in the network assuming pores are
    intersecting trapezoids.

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return intersecting_cones(
        network=network,
        pore_coords=pore_coords,
        throat_coords=throat_coords,
    )


@_geodocs
def hybrid_trapezoids_and_rectangles(
    network, pore_diameter="pore.diameter", throat_coords="throat.coords"
):
    r"""
    Calculates conduit lengths in the network assuming pores are
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
    return hybrid_cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords,
    )


@_geodocs
def pyramids_and_cuboids(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are truncated
    pyramids and throats are cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    return cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter,
    )


@_geodocs
def intersecting_pyramids(
    network, pore_coords="pore.coords", throat_coords="throat.coords"
):
    r"""
    Calculates conduit lengths in the network assuming pores are
    intersecting pyramids.

    Parameters
    ----------
    %(network)s
    %(Pcoords)s
    %(Tcoords)s

    Returns
    -------

    """
    return intersecting_cones(
        network=network,
        pore_coords=pore_coords,
        throat_coords=throat_coords,
    )


@_geodocs
def hybrid_pyramids_and_cuboids(
    network, pore_diameter="pore.diameter", throat_coords="throat.coords"
):
    r"""
    Calculates conduit lengths in the network assuming pores are truncated
    pyramids that intersect. Therefore, the throat is the cross sectional plane where
    two pores meet and has negligible/zero volume.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Tcoords)s

    Returns
    -------

    """
    return hybrid_cones_and_cylinders(
        network=network,
        pore_diameter=pore_diameter,
        throat_coords=throat_coords,
    )


@_geodocs
def cubes_and_cuboids(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are cubes
    and throats are cuboids.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    L_ctc = _get_L_ctc(network)
    D1, Dt, D2 = network.get_conduit_data(pore_diameter.split(".", 1)[-1]).T

    L1 = D1 / 2
    L2 = D2 / 2
    Lt = L_ctc - (L1 + L2)

    # Handle overlapping pores
    mask = (Lt < 0) & (L1 > L2)
    L2[mask] = (L_ctc - L1)[mask]
    mask = (Lt < 0) & (L2 > L1)
    L1[mask] = (L_ctc - L2)[mask]
    Lt = np.maximum(Lt, 1e-15)
    return np.vstack((L1, Lt, L2)).T


@_geodocs
def squares_and_rectangles(
    network, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are squares
    and throats are rectangles.

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
    return cubes_and_cuboids(
        network=network,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter,
    )


# Dealing with errors and exceptions
def _raise_incompatible_data():
    msg = (
        "'spheres_and_cylinders' can only be applied when throat diameter is"
        " smaller than that of adjacent pores."
    )
    raise Exception(msg)


def _get_L_ctc(network):
    """Returns throat spacing if it exists, otherwise calculates it."""
    try:
        L_ctc = network["throat.spacing"]
    except KeyError:
        P12 = network["throat.conns"]
        C1 = network["pore.coords"][P12[:, 0]]
        C2 = network["pore.coords"][P12[:, 1]]
        L_ctc = np.linalg.norm(C1 - C2, axis=1)
    return L_ctc
