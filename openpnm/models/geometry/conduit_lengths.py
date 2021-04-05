import numpy as _np
import openpnm.geometry.GenericGeometry as _GenericGeometry
from .misc import _get_conduit_diameters
from .throat_length import ctc as _ctc

__all__ = [
    "spheres_and_cylinders",
    "circles_and_rectangles",
    "cones_and_cylinders",
    "trapezoids_and_rectangles",
    "pyramids_and_cuboids",
    "cubes_and_cuboids",
    "squares_and_rectangles"
]


def spheres_and_cylinders(
    target: _GenericGeometry,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Calculates conduit lengths in the network assuming pores are spheres
    and throats are cylinders.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    """
    L_ctc = _ctc(target)
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)

    # Handle the case where Dt > Dp
    if (Dt > D1).any() or (Dt > D2).any():
        _raise_incompatible_data()
    L1 = _np.sqrt(D1**2 - Dt**2) / 2
    L2 = _np.sqrt(D2**2 - Dt**2) / 2

    # Handle throats w/ overlapping pores
    _L1 = (4 * L_ctc**2 + D1**2 - D2**2) / (8 * L_ctc)
    mask = L_ctc - 0.5 * (D1 + D2) < 0
    L1[mask] = _L1[mask]
    L2[mask] = (L_ctc - L1)[mask]

    Lt = _np.maximum(L_ctc - (L1 + L2), 1e-15)

    return _np.vstack((L1, Lt, L2)).T


def circles_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are circles
    and throats are rectangles.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return spheres_and_cylinders(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    )


def cones_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are cones
    and throats are cylinders.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    """
    L_ctc = _ctc(target)
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)

    L1 = D1 / 2
    L2 = D2 / 2

    # Handle throats w/ overlapping pores
    _L1 = (4 * L_ctc**2 + D1**2 - D2**2) / (8 * L_ctc)
    mask = L_ctc - 0.5 * (D1 + D2) < 0
    L1[mask] = _L1[mask]
    L2[mask] = (L_ctc - L1)[mask]

    Lt = _np.maximum(L_ctc - (L1 + L2), 1e-15)

    return _np.vstack((L1, Lt, L2)).T


def trapezoids_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are
    trapezoids and throats are rectangles.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    L_ctc = _ctc(target)
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)

    L1 = D1 / 2
    L2 = D2 / 2

    # Handle throats w/ overlapping pores
    _L1 = (4 * L_ctc**2 + D1**2 - D2**2) / (8 * L_ctc)
    mask = L_ctc - 0.5 * (D1 + D2) < 0
    L1[mask] = _L1[mask]
    L2[mask] = (L_ctc - L1)[mask]

    Lt = _np.maximum(L_ctc - (L1 + L2), 1e-15)

    return _np.vstack((L1, Lt, L2)).T


def pyramids_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are truncated
    pyramids and throats are cuboids.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    """
    return cones_and_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )


def cubes_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are cubes
    and throats are cuboids.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    """
    L_ctc = _ctc(target)
    D1, Dt, D2 = _get_conduit_diameters(target, pore_diameter, throat_diameter)

    L1 = D1 / 2
    L2 = D2 / 2
    Lt = L_ctc - (L1 + L2)

    # Handle overlapping pores
    mask = (Lt < 0) & (L1 > L2)
    L2[mask] = (L_ctc - L1)[mask]
    mask = (Lt < 0) & (L2 > L1)
    L1[mask] = (L_ctc - L2)[mask]
    Lt = _np.maximum(Lt, 1e-15)

    return _np.vstack((L1, Lt, L2)).T


def squares_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter"
):
    r"""
    Calculates conduit lengths in the network assuming pores are squares
    and throats are rectangles.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

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

    Returns
    -------
    ndarray
        Array (Nt by 3) containing conduit lengths' values for conduits in
        the given geometry object. The array is formatted as:
            [pore1, throat, pore2].T

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return cubes_and_cuboids(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )


# Dealing with errors and exceptions
def _raise_incompatible_data():
    raise Exception(
        "'spheres_and_cylinders' can only be applied when throat diameter is"
        " smaller than that of adjacent pores."
    )
