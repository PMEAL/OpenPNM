import numpy as _np
import openpnm.models.geometry.conduit_lengths as _conduit_lengths
from openpnm.utils import Docorator


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
docstr = Docorator()


@docstr.get_sections(base='models.geometry.hydraulic_size_factor',
                     sections=['Parameters', 'Returns', 'Notes'])
@docstr.dedent
def spheres_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    spheres and throats are cylinders.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.pdia)s
    %(models.geometry.tdia)s

    Returns
    -------
    size_factors : dict
        A dictionary containing the hydraulic size factors, which can be
        accessed via the dict keys 'pore1', 'pore2', and 'throat'.

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
    D1, Dt, D2 = target.get_conduit_data(poreprop=pore_diameter,
                                         throatprop=throat_diameter).T
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        target=target,
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

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def circles_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    circles and throats are rectangles.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    D1, Dt, D2 = target.get_conduit_data(poreprop=pore_diameter,
                                         throatprop=throat_diameter).T
    L1, Lt, L2 = _conduit_lengths.circles_and_rectangles(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / (D1**2 * _np.sqrt(D1**2 - 4 * L1**2))
    F2 = L2 / (D2**2 * _np.sqrt(D2**2 - 4 * L2**2))
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def cones_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated cones and throats are cylinders.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    D1, Dt, D2 = target.get_conduit_data(poreprop=pore_diameter,
                                         throatprop=throat_diameter).T
    L1, Lt, L2 = _conduit_lengths.cones_and_cylinders(
        target=target,
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

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def trapezoids_and_rectangles(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    trapezoids and throats are rectangles.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = target.get_conduit_data(poreprop=pore_diameter,
                                         throatprop=throat_diameter).T
    L1, Lt, L2 = _conduit_lengths.cones_and_cylinders(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / 2 * (D1 + Dt) / (D1 * Dt)**2
    F2 = L2 / 2 * (D2 + Dt) / (D2 * Dt)**2
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def pyramids_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
):
    r"""
    Computes hydraulic size factors for conduits assuming pores are
    truncated pyramids and throats are cuboids.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    D1, Dt, D2 = target.get_conduit_data(
        poreprop=pore_diameter,
        throatprop=throat_diameter
    ).T
    L1, Lt, L2 = _conduit_lengths.pyramids_and_cuboids(
        target=target,
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

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def cubes_and_cuboids(
    target,
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
    %(models.geometry.hydraulic_size_factor.parameters)s
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    D1, Dt, D2 = target.get_conduit_data(
        poreprop=pore_diameter,
        throatprop=throat_diameter
    ).T
    L1, Lt, L2 = _conduit_lengths.cubes_and_cuboids(
        target=target,
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

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def squares_and_rectangles(
    target,
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
    %(models.geometry.hydraulic_size_factor.parameters)s
    pore_aspect : list
        Aspect ratio of the pores
    throat_aspect : list
        Aspect ratio of the throats

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    D1, Dt, D2 = target.get_conduit_data(
        poreprop=pore_diameter,
        throatprop=throat_diameter
    ).T
    L1, Lt, L2 = _conduit_lengths.squares_and_rectangles(
        target=target,
        pore_diameter=pore_diameter,
        throat_diameter=throat_diameter
    ).T

    # Fi is the integral of (1/A^3) dx, x = [0, Li]
    F1 = L1 / D1**3
    F2 = L2 / D2**3
    Ft = Lt / Dt**3

    # S is 1 / (12 * F)
    S1, St, S2 = [1 / (Fi * 12) for Fi in [F1, Ft, F2]]

    return {"pore1": S1, "throat": St, "pore2": S2}


@docstr.dedent
def intersecting_cones(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes hydraulic size factors for conduits of intersecting cones.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    raise NotImplementedError


@docstr.dedent
def intersecting_trapezoids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes hydraulic size factors for conduits of intersecting
    trapezoids.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    raise NotImplementedError

@docstr.dedent
def intersecting_pyramids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
):
    r"""
    Computes hydraulic size factors for conduits of intersecting pyramids.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s
    midpoint : str, optional
        Dictionary key of the midpoint values.

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    raise NotImplementedError


@docstr.dedent
def ncylinders_in_series(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    n=5,
):
    r"""
    Computes hydraulic size factors for conduits of spheres and cylinders
    with the spheres approximated as N cylinders in series.

    Parameters
    ----------
    %(models.geometry.hydraulic_size_factor.parameters)s
    n : int
        Number of cylindrical divisions for each pore

    Returns
    -------
    %(models.geometry.hydraulic_size_factor.returns)s

    Notes
    -----
    %(models.geometry.hydraulic_size_factor.notes)s

    """
    D1, Dt, D2 = target.get_conduit_data(
        poreprop=pore_diameter,
        throatprop=throat_diameter
    ).T
    # Ensure throats are never bigger than connected pores
    Dt = _np.minimum(Dt, 0.99 * _np.minimum(D1, D2))
    L1, Lt, L2 = _conduit_lengths.spheres_and_cylinders(
        target=target,
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

    return {"pore1": F1, "throat": Ft, "pore2": F2}
