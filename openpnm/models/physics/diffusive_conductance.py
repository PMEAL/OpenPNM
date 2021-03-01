r"""
Pore-scale models for calculating the diffusive conductance of conduits.
"""
import numpy as _np
import scipy.constants as _const

__all__ = [
    "generic_diffusive",
    "ordinary_diffusion",
    "ordinary_diffusion_2d",
    "mixed_diffusion",
    "taylor_aris_diffusion",
    "classic_ordinary_diffusion",
    "multiphase_diffusion"
]


def generic_diffusive(
    target,
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    size_factors="throat.diffusive_size_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    size_factors: str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    """
    network = target.project.network
    throats = target.throats(target=network)
    conns = network.conns[throats]
    phase = target.project.find_phase(target)
    F = network[size_factors]
    DAB1, DAB2 = phase[pore_diffusivity][conns].T
    DABt = phase[throat_diffusivity]

    if isinstance(F, dict):
        g1 = DAB1 * F[f"{size_factors}.pore1"][throats]
        gt = DABt * F[f"{size_factors}.throat"][throats]
        g2 = DAB2 * F[f"{size_factors}.pore2"][throats]
        gd = 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        gd = DABt * F

    return gd


def ordinary_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This method assumes cylindrical throats with constant cross-section
    area, and the pore area corresponds to the cross-sectional area at the
    largest opening of the pore. Corrections for different shapes and
    variable cross-section area can be imposed by passing the proper
    ``conduit_shape_factors`` argument.

    ``conduit_shape_factor`` depends on the physics of the problem,
    i.e. diffusion-like processes and fluid flow need different shape
    factors.

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
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1 * A1) / L1
    g2 = (D2 * A2) / L2
    gt = (Dt * At) / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    # Apply shape factors and calculate the final conductance
    return (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2) ** (-1)


def ordinary_diffusion_2d(
    target,
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in a 2D network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ). The
    conduit consists of 2 flate parallel plates. See the Notes section.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
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

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area
    can be imposed by passing the proper conduit_shape_factors argument.

    ``conduit_shape_factors`` depends on the physics of the problem, i.e.
    diffusion-like processes and fluid flow need different shape factors.

    """
    # Basically, call the 3D version, but pass area with diameter
    return ordinary_diffusion(
        target=target,
        pore_area=pore_diameter,
        throat_area=throat_diameter,
        pore_diffusivity=pore_diffusivity,
        throat_diffusivity=throat_diffusivity,
        conduit_lengths=conduit_lengths,
        conduit_shape_factors=conduit_shape_factors
    )


def mixed_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    pore_temperature="pore.temperature",
    molecular_weight="pore.molecular_weight",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in network with
    Knudsen correction.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ). See
    Notes section for the limitations of this method.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    pore_temperature : str
        Dictionary key of the pore temperature values.
    molecular_weigth : str
        Dictionary key of the pore molecular weight values.
    conduit_lengths : str
        Dictionary key of the conduit lengths' values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function is ONLY suitable for dilute mixtures and NOT those with
    concentrated species.

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
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
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
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Calculating Knudsen diffusivity
    d1, d2 = network[pore_diameter][cn].T
    dt = network[throat_diameter][throats]
    MW1, MW2 = phase[molecular_weight][cn].T
    MWt = phase.interpolate_data(propname=molecular_weight)[throats]
    T1, T2 = phase[pore_temperature][cn].T
    Tt = phase.interpolate_data(propname=pore_temperature)[throats]
    DK1 = d1 / 3 * (8 * _const.R * T1 / _const.pi / MW1) ** 0.5
    DK2 = d2 / 3 * (8 * _const.R * T2 / _const.pi / MW2) ** 0.5
    DKt = dt / 3 * (8 * _const.R * Tt / _const.pi / MWt) ** 0.5
    # Calculate mixed diffusivity
    D1e = (1 / DK1 + 1 / D1) ** (-1)
    D2e = (1 / DK2 + 1 / D2) ** (-1)
    Dte = (1 / DKt + 1 / Dt) ** (-1)
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1e * A1) / L1
    g2 = (D2e * A2) / L2
    gt = (Dte * At) / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    # Apply shape factors and calculate the final conductance
    return (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2) ** (-1)


def taylor_aris_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diffusivity="pore.diffusivity",
    pore_pressure="pore.pressure",
    throat_hydraulic_conductance="throat.hydraulic_conductance",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in network
    considering the Taylor-Aris effect (effect of flow on diffusion).

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore )

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    pore_pressure : str
        Dictionary key of the pore pressure values.
    throat_hydraulic_conductance : str
        Dictionary key of the throat hydraulic_conductance values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties are
    already calculated.

    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area
    can be imposed by passing the proper ``conduit_shape_factors``
    argument.

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
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Fetch properties for calculating Peclet
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]
    Qt = -gh * _np.diff(P[cn], axis=1).squeeze()
    # Find fluid velocity in elements
    u1 = Qt / A1
    u2 = Qt / A2
    ut = Qt / At
    # Find Peclet number in elements
    Pe1 = u1 * ((4 * A1 / _np.pi) ** 0.5) / D1
    Pe2 = u2 * ((4 * A2 / _np.pi) ** 0.5) / D2
    Pet = ut * ((4 * At / _np.pi) ** 0.5) / Dt
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = D1 * (1 + (Pe1 ** 2) / 192) * A1 / L1
    g2 = D2 * (1 + (Pe2 ** 2) / 192) * A2 / L2
    gt = Dt * (1 + (Pet ** 2) / 192) * At / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    # Apply shape factors and calculate the final conductance
    return (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2) ** (-1)


def classic_ordinary_diffusion(
    target,
    pore_molar_density="pore.molar_density",
    pore_diffusivity="pore.diffusivity",
    pore_diameter="pore.diameter",
    throat_diameter="throat.diameter",
    pore_area="pore.area",
    throat_area="throat.area",
    throat_length="throat.length",
    shape_factor="throat.shape_factor"
):
    r"""
    Legacy method that calculates the diffusive conductance of conduits.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore )

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_molar_density : str
        Dictionary key pointing to pore molar density values.
    pore_diffusivity : str
        Dictionary key pointing to pore diffusivity values.
    pore_diameter : str
        Dictionary key pointing to pore diameter values.
    throat_diameter : str
        Dictionary key pointing to throat diameter values.
    pore_area : str
        Dictionary key pointing to pore area values.
    throat_area : str
        Dictionary key pointing to pore area values.
    throat_length : str
        Dictionary key pointing to throat length values.
    shape_factor : str
        Dictionary key pointing to shape factor values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This method requires that all the necessary phase properties already
    be calculated.

    This method calculates the specified property for the *entire* network
    then extracts the values for the appropriate throats at the end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network["throat.conns"]
    # Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Get the properties of every throat
    tdia = network[throat_diameter]
    tarea = _np.pi * (tdia / 2) ** 2
    tlen = network[throat_length]
    # Interpolate pore phase property values to throats
    DABt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    ct = phase.interpolate_data(propname=pore_molar_density)[throats]
    # Get pore lengths
    plen1 = 0.5 * pdia[Ps[:, 0]]
    plen2 = 0.5 * pdia[Ps[:, 1]]
    # Remove any non-positive lengths
    plen1[plen1 <= 1e-12] = 1e-12
    plen2[plen2 <= 1e-12] = 1e-12
    # Find g for half of pore 1
    gp1 = ct * DABt * parea[Ps[:, 0]] / plen1
    gp1[_np.isnan(gp1)] = _np.inf
    gp1[~(gp1 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = ct * DABt * parea[Ps[:, 1]] / plen2
    gp2[_np.isnan(gp2)] = _np.inf
    gp2[~(gp2 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat, remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except KeyError:
        sf = _np.ones(network.num_throats())
    sf[_np.isnan(sf)] = 1.0
    gt = (1 / sf) * ct * DABt * tarea / tlen
    # Set 0 conductance pores (boundaries) to inf
    gt[~(gt > 0)] = _np.inf
    value = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
    return value


def multiphase_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
    partition_coef_global="throat.partition_coef.all"
):
    r"""
    Calculates the diffusive conductance of conduits in network.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array (Nt by 2) containing diffusive conductance values for
        conduits in the geometry attached to the given physics object.

    Notes
    -----
    This method assumes that phase["partition_coef"] contains information on
    binary phase partitioning. See ``MultiPhase`` class documentation for more
    information.

    """
    network = target.network
    throats = target.throats(target=network)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]

    # Getting equivalent areas
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    # Getting conduit lengths
    L1 = network[f"{conduit_lengths}.pore1"][throats]
    Lt = network[f"{conduit_lengths}.throat"][throats]
    L2 = network[f"{conduit_lengths}.pore2"][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]

    # Find g for half of pore 1, throat, and half of pore 2 + apply shape factors
    g1 = (D1 * A1) / L1 * SF1
    g2 = (D2 * A2) / L2 * SF2
    gt = (Dt * At) / Lt * SFt

    # Ensure infinite conductance for elements with zero length
    for gi, Li in zip([g1, gt, g2], [L1, Lt, L2]):
        gi[Li == 0] = _np.inf

    # Apply Henry's partitioning coefficient
    # Note 1: m12 = G21*c1 - G12*c2
    K12 = phase[partition_coef_global][throats]
    G21 = (1/g1 + 0.5/gt + K12 * (1/g2 + 0.5/gt)) ** (-1)
    G12 = K12 * G21

    return _np.vstack((G12, G21)).T
