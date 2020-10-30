import numpy as _np
from numpy import pi as _pi
from numpy import arctanh as _atanh
from numpy import sqrt as _sqrt
from numpy.linalg import norm as _norm


def _calc_conduit_length_spheres_cylinders_2D(
    target, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    return _calc_conduit_length_spheres_cylinders(
        target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
    )


def _calc_conduit_length_spheres_cylinders(
    target, pore_diameter="pore.diameter", throat_diameter="throat.diameter"
):
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]
    coords = network["pore.coords"]
    # throat endpoints
    C1 = coords[cn[:, 0]]
    C2 = coords[cn[:, 1]]
    ctc = _norm(C1 - C2, axis=1)
    L = ctc + 1e-15
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _np.zeros_like(L)
    L2 = _np.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _np.sqrt(D1[~mask] ** 2 - Dt[~mask] ** 2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _np.sqrt(D2[~mask] ** 2 - Dt[~mask] ** 2) / 2
    unit_vec_P1T = (coords[cn[:, 1]] - coords[cn[:, 0]]) / L[:, None]
    unit_vec_P2T = -1 * unit_vec_P1T
    # Find throat endpoints
    EP1 = coords[cn[:, 0]] + L1[:, None] * unit_vec_P1T
    EP2 = coords[cn[:, 1]] + L2[:, None] * unit_vec_P2T
    # Handle throats w/ overlapping pores
    L1 = (4 * L ** 2 + D1 ** 2 - D2 ** 2) / (8 * L)
    L2 = (4 * L ** 2 + D2 ** 2 - D1 ** 2) / (8 * L)
    h = (2 * _np.sqrt(D1 ** 2 / 4 - L1 ** 2)).real
    overlap = L - 0.5 * (D1 + D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = (coords[cn[:, 0]] + L1[:, None] * unit_vec_P1T)[mask]
    EP2[mask] = (coords[cn[:, 1]] + L2[:, None] * unit_vec_P2T)[mask]
    # Calculate conduit lengths
    Lt = _norm(EP1 - EP2, axis=1)
    L1 = _norm(C1 - EP1, axis=1)
    L2 = _norm(C2 - EP2, axis=1)
    return L1, L2, Lt


def spheres_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    conduit_lengths=None,
):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders

    Parameter
    ---------
    target: OpenPNM Geometry object


    Notes
    -----
    The diffusive shape coefficient is the geometrical part of the pre-factor
    in Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    if conduit_lengths is not None:
        L1 = network[conduit_lengths + ".pore1"][throats]
        L2 = network[conduit_lengths + ".pore2"][throats]
        Lt = network[conduit_lengths + ".throat"][throats]
    else:
        L1, L2, Lt = _calc_conduit_length_spheres_cylinders(
            target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
        )
    # F is Integral(1/A) dx , x : 0 --> L
    if (_np.sum(D1 <= 2 * L1) != 0) or (_np.sum(D2 <= 2 * L2) != 0):
        raise Exception("Some throats are too short, add spherical_pores endpoint model")
    F1 = 2 / (D1 * _pi) * _atanh(2 * L1 / D1)
    F2 = 2 / (D2 * _pi) * _atanh(2 * L2 / D2)
    Ft = Lt / (_pi / 4 * Dt ** 2)
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft
    vals = {"pore1": g1, "throat": gt, "pore2": g2}
    return vals


def spheres_and_cylinders_2D(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    conduit_lengths=None,
):
    r"""
    Compute diffusive shape coefficient for conduits of spheres and cylinders
    in 2D

    Parameter
    ---------
    target: OpenPNM object

    Notes
    -----
    The diffusive shape coefficient is the geometrical part of the pre-factor
    in Fick's law:

    .. math::

        n_A = S_{diffusive} D_{AB} \Delta C_A

    Thus :math:`S_{diffusive}` represents the combined effect of the area and
    length of the *conduit*, which consists of a throat and 1/2 of the pore
    on each end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    if conduit_lengths is not None:
        L1 = network[conduit_lengths + ".pore1"][throats]
        L2 = network[conduit_lengths + ".pore2"][throats]
        Lt = network[conduit_lengths + ".throat"][throats]
    else:
        L1, L2, Lt = _calc_conduit_length_spheres_cylinders_2D(
            target, pore_diameter=pore_diameter, throat_diameter=throat_diameter
        )
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1 = 0.5 * _atanh(2 * L1 / _sqrt(D1 ** 2 - 4 * L1 ** 2))
    F2 = 0.5 * _atanh(2 * L2 / _sqrt(D2 ** 2 - 4 * L2 ** 2))
    Ft = Lt / Dt
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft
    vals = {"pore1": g1, "throat": gt, "pore2": g2}
    return vals


def cones_and_cylinders(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_area="pore.area",
    throat_area="throat.area",
    conduit_lengths="throat.conduit_lengths",
):
    r"""
    Compute diffusive shape coefficient assuming pores are truncated pyramids
    and throats are cylinders.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    A dictionary containing the diffusive_shape_coefficient, which can be
    accessed via the dict keys 'pore1', 'pore2', and 'throat'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network["throat.conns"][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    At = network[throat_area][throats]
    L1 = network[conduit_lengths + ".pore1"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1 = 4 * L1 / (D1 * Dt * _pi)
    F2 = 4 * L2 / (D2 * Dt * _pi)
    Ft = Lt / At
    g1, g2, gt = 1 / F1, 1 / F2, 1 / Ft
    vals = {"pore1": g1, "throat": gt, "pore2": g2}
    return vals


def pyramids_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    return_elements=False,
):
    r"""

    """
    pass


def cubes_and_cuboids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_aspect=[1, 1, 1],
    throat_aspect=[1, 1, 1],
    return_elements=False,
):
    r"""

    """
    pass


def intersecting_cones(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
    return_elements=False,
):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns] / 2).T
    Rt = target[throat_diameter] / 2
    Lt = 0
    L1 = Lt - R1
    L2 = Lt - R2

    alpha1 = (R1 - Rt) / L1
    beta1 = 1 / (1 / (Rt ** 3) - 1 / (R1 ** 3))
    alpha2 = (R2 - Rt) / L2
    beta2 = 1 / (1 / (Rt ** 3) - 1 / (R2 ** 3))
    g1 = (3 * alpha1 * _pi / 8) * beta1
    g2 = (3 * alpha2 * _pi / 8) * beta2
    gt = _pi * Rt ** 4 / (8 * Lt)

    if return_elements:
        g = {"pore1": g1, "throat": gt, "pore2": g2}
    else:
        g = (1 / g1 + 1 / gt + 1 / g2) ** -1
    return g


def intersecting_pyramids(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    midpoint=None,
    return_elements=False,
):
    r"""
    """
    pass


def ncylinders_in_series(
    target,
    pore_diameter="pore.equivalent_diameter",
    throat_diameter="throat.equivalent_diameter",
    throat_length="throat.length",
    n=5,
    return_elements=True,
):
    r"""
    Computes the shape coefficient of pores as N cylinders in series

    Parameters
    ----------
    target : OpenPNM Geometry object
        The object to which this model applies
    pore_diameter : str
        Dictionary key pointing to the pore diameter values
    throat_diameter : str
        Dictionary key pointing to the throat diameter values
    throat_length : str
    """
    project = target.project
    network = project.network
    P12 = network["throat.conns"]
    Pdia1, Pdia2 = network[pore_diameter][P12].T
    Tdia = network[throat_diameter]
    # Ensure throats are never bigger than connected pores
    Tdia = _np.minimum(Tdia, 0.99 * _np.minimum(Pdia1, Pdia2))
    Plen1 = Pdia1 / 2 * (_np.cos(_np.arcsin(Tdia / Pdia1)))
    Plen2 = Pdia2 / 2 * (_np.cos(_np.arcsin(Tdia / Pdia2)))
    Lcyl1 = _np.linspace(0, Plen1, num=n, endpoint=False)
    Lcyl2 = _np.linspace(0, Plen2, num=n, endpoint=False)
    Rcyl1 = Pdia1 / 2 * _np.sin(_np.arccos(Lcyl1 / (Pdia1 / 2)))
    Rcyl2 = Pdia2 / 2 * _np.sin(_np.arccos(Lcyl2 / (Pdia2 / 2)))
    gtemp = (_pi * Rcyl1 ** 2 / (Plen1 / n)).T
    g1 = 1 / _np.sum(1 / gtemp, axis=1)
    gtemp = (_pi * Rcyl2 ** 2 / (Plen2 / n)).T
    g2 = 1 / _np.sum(1 / gtemp, axis=1)
    Tlen = network[throat_length]
    gt = (_pi * (Tdia / 2) ** 2 / (Tlen)).T
    if return_elements:
        result = {"pore1": g1, "throat": gt, "pore2": g2}
    else:
        result = 1 / (1 / g1 + 1 / gt + 1 / g2)
    return result
