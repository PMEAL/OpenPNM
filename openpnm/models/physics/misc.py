from scipy import pi as pi


def generic_conductance(target, mechanism, pore_diffusivity,
                        throat_diffusivity, pore_area, throat_area,
                        conduit_lengths, conduit_shape_factors):
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
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    shape_factor : string
        Dictionary key of the conduit shape factor values

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
    L1 = network[conduit_lengths+'.pore1'][throats]
    Lt = network[conduit_lengths+'.throat'][throats]
    L2 = network[conduit_lengths+'.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = network[conduit_shape_factors+'.pore1'][throats]
        SFt = network[conduit_shape_factors+'.throat'][throats]
        SF2 = network[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1
    # Interpolate pore phase property values to throats
    try:
        Dt = phase[throat_diffusivity]
    except KeyError:
        Dt = phase.interpolate_data(propname=pore_diffusivity)
    try:
        Dp = phase[pore_diffusivity]
    except KeyError:
        Dp = phase.interpolate_data(propname=throat_diffusivity)
    # Find g for half of pore 1, throat, and half of pore 2
    if mechanism == 'flow':
        gp1 = A1**2/(8*pi*Dp[cn[:, 0]]*L1)
        gp2 = A2**2/(8*pi*Dp[cn[:, 1]]*L2)
        gt = At**2/(8*pi*Dt*Lt)
    elif mechanism == 'diffusion':
        gp1 = Dp[cn[:, 0]] * A1 / L1
        gp2 = Dp[cn[:, 1]] * A2 / L2
        gt = Dt[throats] * At / Lt
    else:
        raise Exception('Unknown keyword for "mechanism", can only be' +
                        ' "flow" or "diffusion"')
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/gp1/SF1 + 1/gp2/SF2)**(-1)
