r"""

.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille
.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille_power_law
.. autofunction:: openpnm.models.physics.hydraulic_conductance.generic_conductance

"""

import scipy as _sp
import numpy as np


def hagen_poiseuille(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
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

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

    """
    return generic_conductance(target=target, transport_type='flow',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)


def hagen_poiseuille_power_law(
        target,
        pore_area='pore.area',
        throat_area='throat.area',
        pore_viscosity='pore.viscosity',
        throat_viscosity='throat.viscosity',
        conduit_lengths='throat.conduit_lengths',
        conduit_shape_factors='throat.flow_shape_factors',
        pore_consistency='pore.consistency',
        throat_consistency='throat.consistency',
        pore_flow_index='pore.flow_index',
        throat_flow_index='throat.flow_index',
        pore_pressure='pore.pressure'):
    r"""
    Calculate the hydraulic conductance of conduits in network (assuming a non
    Newtonian fluid whose viscosity obeys a power law), where a
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

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    pore_consistency : string
        Dictionary key of the pore fluid consistency values

    throat_consistency : string
        Dictionary key of the throat fluid consistency values

    pore_flow_index : string
        Dictionary key of the pore fluid flow index values

    throat_flow_index : string
        Dictionary key of the throat fluid flow index values

    pore_pressure : string
        Dictionary key of the pore pressure values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

    """
    return generic_conductance(target=target, transport_type='flow_power_law',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors,
                               pore_consistency=pore_consistency,
                               throat_consistency=throat_consistency,
                               pore_flow_index=pore_flow_index,
                               throat_flow_index=throat_flow_index,
                               pore_pressure=pore_pressure)


def valvatne_blunt(target,
                   pore_viscosity='pore.viscosity',
                   throat_viscosity='throat.viscosity',
                   pore_shape_factor='pore.shape_factor',
                   throat_shape_factor='throat.shape_factor',
                   pore_area='pore.area',
                   throat_area='throat.area',
                   conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the single phase hydraulic conductance of conduits in network,
    where a conduit is ( 1/2 pore - full throat - 1/2 pore ) according to [1].
    Function has been adapted for use with the Statoil imported networks and
    makes use of the shape factor in these networks to apply Hagen-Poiseuille
    flow for conduits of different shape classes: Triangular, Square and
    Circular [2].

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    pore_shape_factor : string
        Dictionary key of the pore geometric shape factor values

    throat_shape_factor : string
        Dictionary key of the throat geometric shape factor values

    pore_area : string
        Dictionary key of the pore area values. The pore area is calculated
        using following formula:
            pore_area = (pore_radius ** 2) / (4 * pore_shape_factor)
        Where theoratical value of pore_shape_factor in circular tube is
        calculated using following formula:
            pore_shape_factor = pore_area / perimeter **2 = 1/4π

    throat_area : string
        Dictionary key of the throat area values. The throat area is calculated
        using following formula:
            throat_area = (throat_radius ** 2) / (4 * throat_shape_factor)
        Where theoratical value of throat_shape_factor in circular tube is
        calculated using following formula:
            throat_shape_factor = throat_area / perimeter **2 = 1/4π

    conduit_lengths : string
        Dictionary key of the throat conduit lengths

    Returns
    -------
    g : ND-array
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    References
    ----------
    [1] Valvatne, Per H., and Martin J. Blunt. "Predictive pore‐scale modeling
    of two‐phase flow in mixed wet media." Water Resources Research 40,
    no. 7 (2004).
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
    tri = Gt <= np.sqrt(3)/36.0
    circ = Gt >= 0.07
    square = ~(tri | circ)
    ntri = np.sum(tri)
    nsquare = np.sum(square)
    ncirc = np.sum(circ)
    kt = np.ones_like(Gt)
    kt[tri] = 3.0/5.0
    kt[square] = 0.5623
    kt[circ] = 0.5
    # Pore Portions
    Gp = network[pore_shape_factor]
    tri = Gp <= np.sqrt(3)/36.0
    circ = Gp >= 0.07
    square = ~(tri | circ)
    ntri += np.sum(tri)
    nsquare += np.sum(square)
    ncirc += np.sum(circ)
    kp = np.ones_like(Gp)
    kp[tri] = 3.0/5.0
    kp[square] = 0.5623
    kp[circ] = 0.5
    gp = kp*(network[pore_area]**2)*Gp/mu_p
    gt = kt*(network[throat_area]**2)*Gt/mu_t
    conns = network['throat.conns']
    l1 = network[conduit_lengths + '.pore1']
    lt = network[conduit_lengths + '.throat']
    l2 = network[conduit_lengths + '.pore2']
    # Resistors in Series
    value = (l1/gp[conns[:, 0]] +
             lt/gt +
             l2/gp[conns[:, 1]])
    return 1/value


def classic_hagen_poiseuille(target,
                             pore_diameter='pore.diameter',
                             pore_viscosity='pore.viscosity',
                             throat_length='throat.length',
                             throat_diameter='throat.diameter',
                             shape_factor='throat.shape_factor',
                             **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming cylindrical
    geometry using the Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object

    Notes
    -----
    This function calculates the specified property for the *entire* network
    then extracts the values for the appropriate throats at the end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    phase = target.project.find_phase(target)
    mup = phase[pore_viscosity]
    mut = phase.interpolate_data(propname=pore_viscosity)[throats]
    pdia = network[pore_diameter]
    # Get pore lengths
    plen1 = (0.5*pdia[Ps[:, 0]])
    plen2 = (0.5*pdia[Ps[:, 1]])
    # Remove any non-positive lengths
    plen1[plen1 <= 1e-12] = 1e-12
    plen2[plen2 <= 1e-12] = 1e-12
    # Find g for half of pore 1
    gp1 = _sp.pi*(pdia[Ps[:, 0]])**4/(128*plen1*mut)
    gp1[_sp.isnan(gp1)] = _sp.inf
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf

    # Find g for half of pore 2
    gp2 = _sp.pi*(pdia[Ps[:, 1]])**4/(128*plen2*mut)
    gp2[_sp.isnan(gp2)] = _sp.inf
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except KeyError:
        sf = _sp.ones(network.num_throats())
    sf[_sp.isnan(sf)] = 1.0
    gt = (1/sf)*_sp.pi*(tdia)**4/(128*tlen*mut)
    gt[~(gt > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return value


def generic_conductance(target, transport_type, pore_area, throat_area,
                        pore_diffusivity, throat_diffusivity,
                        conduit_lengths, conduit_shape_factors, **kwargs):
    r"""
    Calculate the generic conductance (could be mass, thermal, electrical,
    ionic, or hydraylic) of conduits in the network, where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    transport_type : string
        The transport type.  Options are:

        *'flow'* - For Newtonian fluids

        *'power_law'* - For power-law fluids

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
    g : ND-array
        Array containing conductance values for conduits in the geometry
        attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper shape factor.

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
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    try:
        Dt = phase[throat_diffusivity][throats]
    except KeyError:
        Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    try:
        D1 = phase[pore_diffusivity][cn[:, 0]]
        D2 = phase[pore_diffusivity][cn[:, 1]]
    except KeyError:
        D1 = phase.interpolate_data(propname=throat_diffusivity)[cn[:, 0]]
        D2 = phase.interpolate_data(propname=throat_diffusivity)[cn[:, 1]]
    # Find g for half of pore 1, throat, and half of pore 2
    pi = _sp.pi
    if transport_type == 'flow':
        g1[m1] = A1[m1]**2 / (8*pi*D1*L1)[m1]
        g2[m2] = A2[m2]**2 / (8*pi*D2*L2)[m2]
        gt[mt] = At[mt]**2 / (8*pi*Dt*Lt)[mt]
    elif transport_type == 'flow_power_law':
        for k, v in kwargs.items():
            if k == 'pore_consistency':
                pore_consistency = v
            elif k == 'throat_consistency':
                throat_consistency = v
            elif k == 'pore_flow_index':
                pore_flow_index = v
            elif k == 'throat_flow_index':
                throat_flow_index = v
            elif k == 'pore_pressure':
                pore_pressure = v

        # Check if pressure field exists
        try:
            phase[pore_pressure]
        except KeyError:
            phase[pore_pressure] = 0
        P = phase[pore_pressure]

        # Interpolate pore phase property values to throats
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
        dP1 = _sp.absolute(P[cn[:, 0]]-Pt)
        dP2 = _sp.absolute(P[cn[:, 1]]-Pt)
        dPt = _sp.absolute(_sp.diff(P[cn], axis=1).squeeze())

        # Apparent viscosities
        mu1 = (dP1**(1-1/n1)[m1] * C1**(1/n1)[m1] / ((4*n1/(3*n1+1)) *
               (2*L1/((A1/pi)**0.5))**(1-1/n1))[m1])

        mu2 = (dP2**(1-1/n2)[m2] * C2**(1/n2)[m2] / ((4*n2/(3*n2+1)) *
               (2*L2/((A2/pi)**0.5))**(1-1/n2))[m2])

        mut = (dPt**(1-1/nt)[mt] * Ct**(1/nt)[mt] / ((4*nt/(3*nt+1)) *
               (2*Lt/((At/pi)**0.5))**(1-1/nt))[mt])

        # Bound the apparent viscosity
        vis_min = 1e-08
        vis_max = 1e+04
        mu1[mu1 < vis_min] = vis_min
        mu1[mu1 > vis_max] = vis_max
        mu2[mu2 < vis_min] = vis_min
        mu2[mu2 > vis_max] = vis_max
        mut[mut < vis_min] = vis_min
        mut[mut > vis_max] = vis_max

        phase['throat.viscosity_eff'] = mut

        g1[m1] = A1[m1]**2 / ((8*pi*L1)[m1]*mu1)
        g2[m2] = A2[m2]**2 / ((8*pi*L2)[m2]*mu2)
        gt[mt] = At[mt]**2 / ((8*pi*Lt)[mt]*mut)
    else:
        raise Exception('Unknown keyword for "transport_type", can only be' +
                        ' "flow" or "flow_power_law"')
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)
