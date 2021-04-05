r"""
Pore-scale models for calculating hydraulic conductance of conduits.
"""
import numpy as _np

__all__ = [
    "generic_hydraulic",
    "hagen_poiseuille",
    "hagen_poiseuille_2d",
    "hagen_poiseuille_power_law",
    "valvatne_blunt",
    "classic_hagen_poiseuille"
]


def generic_hydraulic(
        target,
        pore_viscosity='pore.viscosity',
        throat_viscosity='throat.viscosity',
        size_factors='throat.hydraulic_size_factors'
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    size_factors: str
        Dictionary key of the conduit hydraulic size factors' values.

    Returns
    -------
    ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    """
    network = target.project.network
    throats = target.throats(target=network)
    conns = network.conns[throats]
    phase = target.project.find_phase(target)
    F = network[size_factors]
    mu1, mu2 = phase[pore_viscosity][conns].T
    mut = phase[throat_viscosity]

    if isinstance(F, dict):
        g1 = F[f"{size_factors}.pore1"][throats] / mu1
        gt = F[f"{size_factors}.throat"][throats] / mut
        g2 = F[f"{size_factors}.pore2"][throats] / mu2
        gh = 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        gh = F / mut

    return gh


def hagen_poiseuille(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.flow_shape_factors"
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit hydraulic (flow) shape factor values.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area
    can be imposed by passing the proper conduit_shape_factors argument.

    ``conduit_shape_factors`` depends on the physics of the problem, i.e.
    diffusion-like processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + ".pore1"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    # Preallocating g
    g1, g2, gt = _np.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _np.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    Dt = phase[throat_viscosity][throats]
    D1, D2 = phase[pore_viscosity][cn].T
    # Find g for half of pore 1, throat, and half of pore 2
    g1[m1] = A1[m1] ** 2 / (8 * _np.pi * D1 * L1)[m1]
    g2[m2] = A2[m2] ** 2 / (8 * _np.pi * D2 * L2)[m2]
    gt[mt] = At[mt] ** 2 / (8 * _np.pi * Dt * Lt)[mt]
    # Apply shape factors and calculate the final conductance
    return 1 / (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2)


def hagen_poiseuille_2d(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.flow_shape_factors",
):
    r"""
    Calculates the hydraulic conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    conduit_lengths : str
        Dictionary key of the conduit lengths' values.
    conduit_shape_factors : str
        Dictionary key of the conduit flow shape factors' values.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    This function requires that all the necessary phase properties already
    be calculated.

    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

    This function assumes rectangular (2D) throats. Corrections for
    different shapes and variable cross-section area can be imposed by
    passing the proper flow_shape_factor argument.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    # Getting pore/throat diameters
    D1 = network[pore_diameter][cn[:, 0]]
    Dt = network[throat_diameter][throats]
    D2 = network[pore_diameter][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + ".pore1"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Getting viscosity values
    mut = phase[throat_viscosity][throats]
    mu1, mu2 = phase[pore_viscosity][cn].T
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = D1 ** 3 / (12 * mu1 * L1)
    g2 = D2 ** 3 / (12 * mu2 * L2)
    gt = Dt ** 3 / (12 * mut * Lt)

    return 1 / (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2)


def hagen_poiseuille_power_law(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_viscosity_min="pore.viscosity_min",
    throat_viscosity_min="throat.viscosity_min",
    pore_viscosity_max="pore.viscosity_max",
    throat_viscosity_max="throat.viscosity_max",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.flow_shape_factors",
    pore_consistency="pore.consistency",
    throat_consistency="throat.consistency",
    pore_flow_index="pore.flow_index",
    throat_flow_index="throat.flow_index",
    pore_pressure="pore.pressure",
):
    r"""
    Calculates the hydraulic conductance of conduits in network, assuming
    a non Newtonian fluid whose viscosity obeys a power law.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_viscosity_min : str
        Dictionary key of the pore minimum viscosity values.
    throat_viscosity_min : str
        Dictionary key of the throat minimum viscosity values.
    pore_viscosity_max : str
        Dictionary key of the pore maximum viscosity values.
    throat_viscosity_max : str
        Dictionary key of the throat maximum viscosity values.
    conduit_lengths : str
        Dictionary key of the conduit lengths' values.
    conduit_shape_factors : str
        Dictionary key of the conduit flow shape factors' values.
    pore_consistency : str
        Dictionary key of the pore fluid consistency values.
    throat_consistency : str
        Dictionary key of the throat fluid consistency values.
    pore_flow_index : str
        Dictionary key of the pore fluid flow index values.
    throat_flow_index : str
        Dictionary key of the throat fluid flow index values.
    pore_pressure : str
        Dictionary key of the pore pressure values.

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area
    can be imposed by passing the proper conduit_shape_factors argument.

    ``conduit_shape_factors`` depends on the physics of the problem, i.e.
    diffusion-like processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + ".pore1"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    # Preallocating g
    g1, g2, gt = _np.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _np.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    pi = _np.pi

    # Check if pressure field exists
    try:
        phase[pore_pressure]
    except KeyError:
        phase[pore_pressure] = 0
    P = phase[pore_pressure]

    # Interpolate pore phase property values to throats
    try:
        mu_mint = phase[throat_viscosity_min][throats]
    except KeyError:
        mu_mint = phase.interpolate_data(propname=pore_viscosity_min)[throats]
    try:
        mu_maxt = phase[throat_viscosity_max][throats]
    except KeyError:
        mu_maxt = phase.interpolate_data(propname=pore_viscosity_max)[throats]
    try:
        Ct = phase[throat_consistency][throats]
    except KeyError:
        Ct = phase.interpolate_data(propname=pore_consistency)[throats]
    try:
        nt = phase[throat_flow_index][throats]
    except KeyError:
        nt = phase.interpolate_data(propname=pore_flow_index)[throats]
    # Interpolate throat phase property values to pores
    try:
        mu_min1 = phase[pore_viscosity_min][cn[:, 0]]
        mu_min2 = phase[pore_viscosity_min][cn[:, 1]]
    except KeyError:
        mu_min1 = phase.interpolate_data(propname=throat_viscosity_min)[cn[:, 0]]
        mu_min2 = phase.interpolate_data(propname=throat_viscosity_min)[cn[:, 1]]
    try:
        mu_max1 = phase[pore_viscosity_max][cn[:, 0]]
        mu_max2 = phase[pore_viscosity_max][cn[:, 1]]
    except KeyError:
        mu_max1 = phase.interpolate_data(propname=throat_viscosity_max)[cn[:, 0]]
        mu_max2 = phase.interpolate_data(propname=throat_viscosity_max)[cn[:, 1]]
    try:
        C1 = phase[pore_consistency][cn[:, 0]]
        C2 = phase[pore_consistency][cn[:, 1]]
    except KeyError:
        C1 = phase.interpolate_data(propname=throat_consistency)[cn[:, 0]]
        C2 = phase.interpolate_data(propname=throat_consistency)[cn[:, 1]]
    try:
        n1 = phase[pore_flow_index][cn[:, 0]]
        n2 = phase[pore_flow_index][cn[:, 1]]
    except KeyError:
        n1 = phase.interpolate_data(propname=throat_flow_index)[cn[:, 0]]
        n2 = phase.interpolate_data(propname=throat_flow_index)[cn[:, 1]]
    # Interpolate pore pressure values to throats
    Pt = phase.interpolate_data(propname=pore_pressure)[throats]

    # Pressure differences dP
    dP1 = _np.absolute(P[cn[:, 0]] - Pt)
    dP2 = _np.absolute(P[cn[:, 1]] - Pt)
    dPt = _np.absolute(_np.diff(P[cn], axis=1).squeeze())

    dP1 = dP1.clip(min=1e-20)
    dP2 = dP2.clip(min=1e-20)
    dPt = dPt.clip(min=1e-20)

    # Apparent viscosities
    mu1, mu2, mut = _np.zeros((3, len(Lt)))

    mu1[m1] = (dP1 ** (1 - 1 / n1) * C1 ** (1 / n1))[m1] / (
        (4 * n1 / (3 * n1 + 1))[m1]
        * (2 * L1[m1] / ((A1[m1] / pi) ** 0.5)) ** (1 - 1 / n1[m1])
    )

    mu2[m2] = (dP2 ** (1 - 1 / n2) * C2 ** (1 / n2))[m2] / (
        (4 * n2 / (3 * n2 + 1))[m2]
        * (2 * L2[m2] / ((A2[m2] / pi) ** 0.5)) ** (1 - 1 / n2[m2])
    )

    mut[mt] = (dPt ** (1 - 1 / nt) * Ct ** (1 / nt))[mt] / (
        (4 * nt / (3 * nt + 1))[mt]
        * (2 * Lt[mt] / ((At[mt] / pi) ** 0.5)) ** (1 - 1 / nt[mt])
    )

    # Bound the apparent viscosity
    mu1[m1] = _np.minimum(_np.maximum(mu1[m1], mu_min1[m1]), mu_max1[m1])
    mu2[m2] = _np.minimum(_np.maximum(mu2[m2], mu_min2[m2]), mu_max2[m2])
    mut[mt] = _np.minimum(_np.maximum(mut[mt], mu_mint[mt]), mu_maxt[mt])

    phase["throat.viscosity_eff.pore1"] = mu1
    phase["throat.viscosity_eff.pore2"] = mu2
    phase["throat.viscosity_eff.throat"] = mut

    g1[m1] = A1[m1] ** 2 / ((8 * pi * L1) * mu1)[m1]
    g2[m2] = A2[m2] ** 2 / ((8 * pi * L2) * mu2)[m2]
    gt[mt] = At[mt] ** 2 / ((8 * pi * Lt) * mut)[mt]

    # Apply shape factors and calculate the final conductance
    return 1 / (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2)


def valvatne_blunt(
    target,
    pore_viscosity="pore.viscosity",
    throat_viscosity="throat.viscosity",
    pore_shape_factor="pore.shape_factor",
    throat_shape_factor="throat.shape_factor",
    pore_area="pore.area",
    throat_area="throat.area",
    conduit_lengths="throat.conduit_lengths",
):
    r"""
    Calculate the single phase hydraulic conductance of conduits in network,
    where a conduit is ( 1/2 pore - full throat - 1/2 pore ) according to [1].
    Function has been adapted for use with the Statoil imported networks and
    makes use of the shape factor in these networks to apply Hagen-Poiseuille
    flow for conduits of different shape classes: Triangular, Square and
    Circular [2].

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    throat_viscosity : str
        Dictionary key of the throat viscosity values.
    pore_shape_factor : str
        Dictionary key of the pore geometric shape factor values.
    throat_shape_factor : str
        Dictionary key of the throat geometric shape factor values.
    pore_area : str
        Dictionary key of the pore area values. The pore area is
        calculated using following formula:
            pore_area = (pore_radius ** 2) / (4 * pore_shape_factor)
        Where theoratical value of pore_shape_factor in circular tube is
        calculated using following formula:
            pore_shape_factor = pore_area / perimeter **2 = 1/4π
    throat_area : str
        Dictionary key of the throat area values. The throat area is
        calculated using following formula:
            throat_area = (throat_radius ** 2) / (4 * throat_shape_factor)
        Where theoratical value of throat_shape_factor in circular tube is
        calculated using following formula:
            throat_shape_factor = throat_area / perimeter ** 2 = 1/4π
    conduit_lengths : str
        Dictionary key of the throat conduit lengths.

    Returns
    -------
    g : ND-array
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    References
    ----------
    [1] Valvatne, Per H., and Martin J. Blunt. "Predictive pore‐scale
    modeling of two‐phase flow in mixed wet media." Water Resources
    Research 40, no. 7 (2004).

    [2] Patzek, T. W., and D. B. Silin (2001), Shape factor and hydraulic
    conductance in noncircular capillaries I. One-phase creeping flow,
    J. Colloid Interface Sci., 236, 295–304.

    """
    network = target.project.network
    mu_p = target[pore_viscosity]
    try:
        mu_t = target[throat_viscosity]
    except KeyError:
        mu_t = target.interpolate_data(pore_viscosity)
    # Throat Portion
    Gt = network[throat_shape_factor]
    tri = Gt <= _np.sqrt(3) / 36.0
    circ = Gt >= 0.07
    square = ~(tri | circ)
    ntri = _np.sum(tri)
    nsquare = _np.sum(square)
    ncirc = _np.sum(circ)
    kt = _np.ones_like(Gt)
    kt[tri] = 3.0 / 5.0
    kt[square] = 0.5623
    kt[circ] = 0.5
    # Pore Portions
    Gp = network[pore_shape_factor]
    tri = Gp <= _np.sqrt(3) / 36.0
    circ = Gp >= 0.07
    square = ~(tri | circ)
    ntri += _np.sum(tri)
    nsquare += _np.sum(square)
    ncirc += _np.sum(circ)
    kp = _np.ones_like(Gp)
    kp[tri] = 3.0 / 5.0
    kp[square] = 0.5623
    kp[circ] = 0.5
    gp = kp * (network[pore_area] ** 2) * Gp / mu_p
    gt = kt * (network[throat_area] ** 2) * Gt / mu_t
    conns = network["throat.conns"]
    l1 = network[conduit_lengths + ".pore1"]
    lt = network[conduit_lengths + ".throat"]
    l2 = network[conduit_lengths + ".pore2"]
    # Resistors in Series
    value = l1 / gp[conns[:, 0]] + lt / gt + l2 / gp[conns[:, 1]]
    return 1 / value


def classic_hagen_poiseuille(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    throat_length="throat.length",
    pore_viscosity="pore.viscosity",
    shape_factor="throat.shape_factor",
):
    r"""
    Legacy method that calculates the hydraulic conductance of conduits
    assuming pores and throats are cylinders in series.

    This method is based on the Hagen-Poiseuille model.

    Parameters
    ----------
    target : _GenericPhysics
        Physics object with which this model is associated.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    throat_length : str
        Dictionary key of the throat length values.
    pore_viscosity : str
        Dictionary key of the pore viscosity values.
    shape_factor: str
        Dictionary key of the throat shape factor values.

    Notes
    -----
    This function calculates the specified property for the *entire* network
    then extracts the values for the appropriate throats at the end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network["throat.conns"]
    # Get properties in every pore in the network
    phase = target.project.find_phase(target)
    mut = phase.interpolate_data(propname=pore_viscosity)[throats]
    pdia = network[pore_diameter]
    # Get pore lengths
    plen1 = 0.5 * pdia[Ps[:, 0]]
    plen2 = 0.5 * pdia[Ps[:, 1]]
    # Remove any non-positive lengths
    plen1[plen1 <= 1e-12] = 1e-12
    plen2[plen2 <= 1e-12] = 1e-12
    # Find g for half of pore 1
    gp1 = _np.pi * (pdia[Ps[:, 0]]) ** 4 / (128 * plen1 * mut)
    gp1[_np.isnan(gp1)] = _np.inf
    gp1[~(gp1 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf

    # Find g for half of pore 2
    gp2 = _np.pi * (pdia[Ps[:, 1]]) ** 4 / (128 * plen2 * mut)
    gp2[_np.isnan(gp2)] = _np.inf
    gp2[~(gp2 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except KeyError:
        sf = _np.ones(network.num_throats())
    sf[_np.isnan(sf)] = 1.0
    gt = (1 / sf) * _np.pi * (tdia) ** 4 / (128 * tlen * mut)
    gt[~(gt > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    value = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
    return value
