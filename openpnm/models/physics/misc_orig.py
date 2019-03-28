from scipy import pi as pi
import scipy as _sp


def generic_conductance(target, transport_type, pore_diffusivity,
                        throat_diffusivity, pore_area, throat_area,
                        conduit_lengths, conduit_shape_factors, **kwargs):
    r"""
    Calculate the generic conductance (could be mass, thermal, electrical,
    or hydraylic) of conduits in the network, where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values.

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values.

    pore_area : string
        Dictionary key of the pore area values.

    throat_area : string
        Dictionary key of the throat area values.

    shape_factor : string
        Dictionary key of the conduit shape factor values.

    Returns
    -------
    g : ndarray
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
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Find g for half of pore 1, throat, and half of pore 2
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
        dP = _sp.absolute(_sp.diff(P[cn], axis=1).squeeze())

        # Interpolate pore phase property values to throats
        try:
            Ct = phase[throat_consistency][throats]
        except KeyError:
            Ct = phase.interpolate_data(propname=pore_consistency)[throats]
        try:
            nt = phase[throat_flow_index][throats]
        except KeyError:
            nt = phase.interpolate_data(propname=pore_flow_index)[throats]
        # Interpolate pore phase property values to pores
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

        g1[m1] = ((A1**2/(8*pi*L1*C1**(1/n1)))[m1] * (4*n1/(3*n1+1))[m1] *
                  (2*L1/((A1/pi)**0.5))**(1-1/n1)[m1] * (dP)**(1/n1-1)[m1])

        g2[m2] = ((A2**2/(8*pi*L2*C2**(1/n2)))[m2] * (4*n2/(3*n2+1))[m2] *
                  (2*L2/((A2/pi)**0.5))**(1-1/n2)[m2] * (dP)**(1/n2-1)[m2])

        gt[mt] = ((At**2/(8*pi*Lt*Ct**(1/nt)))[mt] * (4*nt/(3*nt+1))[mt] *
                  (2*Lt/((At/pi)**0.5))**(1-1/nt)[mt] * (dP)**(1/nt-1)[mt])
        return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)
    elif transport_type == 'diffusion':
        g1[m1] = (D1*A1)[m1] / L1[m1]
        g2[m2] = (D2*A2)[m2] / L2[m2]
        gt[mt] = (Dt*At)[mt] / Lt[mt]
    elif transport_type == 'taylor_aris_diffusion':
        for k, v in kwargs.items():
            if k == 'pore_pressure':
                pore_pressure = v
            elif k == 'throat_hydraulic_conductance':
                throat_hydraulic_conductance = v
        P = phase[pore_pressure]
        gh = phase[throat_hydraulic_conductance]
        Qt = -gh*_sp.diff(P[cn], axis=1).squeeze()

        u1 = Qt[m1]/A1[m1]
        u2 = Qt[m2]/A2[m2]
        ut = Qt[mt]/At[mt]

        Pe1 = u1 * ((4*A1[m1]/_sp.pi)**0.5) / D1[m1]
        Pe2 = u2 * ((4*A2[m2]/_sp.pi)**0.5) / D2[m2]
        Pet = ut * ((4*At[mt]/_sp.pi)**0.5) / Dt[mt]

        g1[m1] = D1[m1]*(1+(Pe1**2)/192)*A1[m1] / L1[m1]
        g2[m2] = D2[m2]*(1+(Pe2**2)/192)*A2[m2] / L2[m2]
        gt[mt] = Dt[mt]*(1+(Pet**2)/192)*At[mt] / Lt[mt]
    else:
        raise Exception('Unknown keyword for "transport_type", can only be' +
                        ' "flow", "diffusion", "taylor_aris_diffusion"' +
                        ' or "ionic"')
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)
