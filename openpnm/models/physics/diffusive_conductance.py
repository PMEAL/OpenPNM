r"""

.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion
.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion_2D
.. autofunction:: openpnm.models.physics.diffusive_conductance.mixed_diffusion
.. autofunction:: openpnm.models.physics.diffusive_conductance.taylor_aris_diffusion
.. autofunction:: openpnm.models.physics.diffusive_conductance.classic_ordinary_diffusion

"""
import numpy as _np
from numpy import pi as _pi
import scipy as _sp
import scipy.constants as const


def ordinary_diffusion(
    target,
    pore_area='pore.area',
    throat_area='throat.area',
    pore_diffusivity='pore.diffusivity',
    throat_diffusivity='throat.diffusivity',
    conduit_lengths='throat.conduit_lengths',
    conduit_shape_factors='throat.poisson_shape_factors'
):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    * This function assumes cylindrical throats with constant cross-section
    area, and the pore area corresponds to the cross-sectional area at the
    largest opening of the pore. Corrections for different shapes and variable
    cross-section area can be imposed by passing the proper
    ``conduit_shape_factors`` argument.

    * ``conduit_shape_factor`` depends on the physics of the problem,
    i.e. diffusion-like processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1*A1) / L1
    g2 = (D2*A2) / L2
    gt = (Dt*At) / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1==0] = _sp.inf
    g2[L2==0] = _sp.inf
    gt[Lt==0] = _sp.inf
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def ordinary_diffusion_2D(
    target,
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter',
    pore_diffusivity='pore.diffusivity',
    throat_diffusivity='throat.diffusivity',
    conduit_lengths='throat.conduit_lengths',
    conduit_shape_factors='throat.poisson_shape_factors'
):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). The conduit consists of 2
    flate parallel plates. See the notes section.

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

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_diameter][cn[:, 0]]
    At = network[throat_diameter][throats]
    A2 = network[pore_diameter][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1*A1) / L1
    g2 = (D2*A2) / L2
    gt = (Dt*At) / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1==0] = _sp.inf
    g2[L2==0] = _sp.inf
    gt[Lt==0] = _sp.inf
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def mixed_diffusion(
    target,
    pore_area='pore.area',
    throat_area='throat.area',
    pore_diameter='pore.diameter',
    throat_diameter='throat.diameter',
    pore_diffusivity='pore.diffusivity',
    pore_temperature='pore.temperature',
    molecular_weight='pore.molecular_weight',
    throat_diffusivity='throat.diffusivity',
    conduit_lengths='throat.conduit_lengths',
    conduit_shape_factors='throat.poisson_shape_factors',
):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ), assuming Knudsen
    diffusivity. See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (0) This function is ONLY suitable for dilute mixtures and NOT those with
    concentrated species.

    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
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
    DK1 = d1/3 * (8*const.R*T1/const.pi/MW1)**0.5
    DK2 = d2/3 * (8*const.R*T2/const.pi/MW2)**0.5
    DKt = dt/3 * (8*const.R*Tt/const.pi/MWt)**0.5
    # Calculate mixed diffusivity
    D1e = (1/DK1 + 1/D1)**(-1)
    D2e = (1/DK2 + 1/D2)**(-1)
    Dte = (1/DKt + 1/Dt)**(-1)
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1e * A1) / L1
    g2 = (D2e * A2) / L2
    gt = (Dte * At) / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1==0] = _sp.inf
    g2[L2==0] = _sp.inf
    gt[Lt==0] = _sp.inf
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def taylor_aris_diffusion(
    target,
    pore_area='pore.area',
    throat_area='throat.area',
    pore_diffusivity='pore.diffusivity',
    pore_pressure='pore.pressure',
    throat_hydraulic_conductance='throat.hydraulic_conductance',
    throat_diffusivity='throat.diffusivity',
    conduit_lengths='throat.conduit_lengths',
    conduit_shape_factors='throat.poisson_shape_factors'
):
    r"""
    Calculate the diffusive conductance of conduits in network considering the
    Taylor-Aris effect (effect of fluid flow on diffusion), where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    pore_pressure : string
        Dictionary key of the pore pressure values

    throat_hydraulic_conductance : string
        Dictionary key of the throat hydraulic_conductance values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values (with Taylor-Aris effect)
        for conduits in the geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties are
    already calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Fetch properties for calculating Peclet
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]
    Qt = -gh*_sp.diff(P[cn], axis=1).squeeze()
    # Find fluid velocity in elements
    u1 = Qt / A1
    u2 = Qt / A2
    ut = Qt / At
    # Find Peclet number in elements
    Pe1 = u1 * ((4 * A1 / _sp.pi)**0.5) / D1
    Pe2 = u2 * ((4 * A2 / _sp.pi)**0.5) / D2
    Pet = ut * ((4 * At / _sp.pi)**0.5) / Dt
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = D1 * (1 + (Pe1**2)/192) * A1 / L1
    g2 = D2 * (1 + (Pe2**2)/192) * A2 / L2
    gt = Dt * (1 + (Pet**2)/192) * At / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1==0] = _sp.inf
    g2[L2==0] = _sp.inf
    gt[Lt==0] = _sp.inf
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def classic_ordinary_diffusion(
    target,
    pore_molar_density='pore.molar_density',
    pore_diffusivity='pore.diffusivity',
    pore_area='pore.area',
    pore_diameter='pore.diameter',
    throat_area='throat.area',
    throat_length='throat.length',
    throat_diameter='throat.diameter',
    shape_factor='throat.shape_factor',
    **kwargs
):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas
    Parameters
    ----------
    network : OpenPNM Network Object
    phase : OpenPNM Phase Object
        The phase of interest
    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.
    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Get the properties of every throat
    tdia = network[throat_diameter]
    tarea = _sp.pi * (tdia / 2) ** 2
    tlen = network[throat_length]
    # Interpolate pore phase property values to throats
    DABt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    ct = phase.interpolate_data(propname=pore_molar_density)[throats]
    # Get pore lengths
    plen1 = (0.5 * pdia[Ps[:, 0]])
    plen2 = (0.5 * pdia[Ps[:, 1]])
    # Remove any non-positive lengths
    plen1[plen1 <= 1e-12] = 1e-12
    plen2[plen2 <= 1e-12] = 1e-12
    # Find g for half of pore 1
    gp1 = ct * DABt * parea[Ps[:, 0]] / plen1
    gp1[_np.isnan(gp1)] = _sp.inf
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = ct * DABt * parea[Ps[:, 1]] / plen2
    gp2[_np.isnan(gp2)] = _sp.inf
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
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
    gt[~(gt > 0)] = _sp.inf
    value = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
    return value


def multiphase_diffusion(
    target,
    pore_area='pore.area',
    throat_area='throat.area',
    pore_diffusivity='pore.diffusivity',
    throat_diffusivity='throat.diffusivity',
    conduit_lengths='throat.conduit_lengths',
    conduit_shape_factors='throat.poisson_shape_factors'
):
    r"""
    Similar to ordinary_diffusion, except it also accounts for Henry's law
    partitioning.

    Notes
    -----
    This method assumes that phase["partition_coef"] contains information on
    binary phase partitioning. See ``MultiPhase`` class documentation for more
    information.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Find g for half of pore 1, throat, and half of pore 2 + apply shape factors
    g1 = (D1*A1) / L1 * SF1
    g2 = (D2*A2) / L2 * SF2
    gt = (Dt*At) / Lt * SFt
    # Ensure infinite conductance for elements with zero length
    g1[L1==0] = _sp.inf
    g2[L2==0] = _sp.inf
    gt[Lt==0] = _sp.inf
    # Get partition coefficient dictionary key from phase settings
    partition_coef = phase.settings["partition_coef"]
    # Apply Henry's partitioning coefficient
    K12 = phase[partition_coef][throats]
    G12 = K12 * (1.0/g1 + 0.5/gt + K12*(1.0/g2 + 0.5/gt)) ** (-1)
    G21 = 1.0/K12 * (1.0/g2 + 0.5/gt + 1.0/K12*(1.0/g1 + 0.5/gt)) ** (-1)
    return _np.vstack((G12, G21)).T
