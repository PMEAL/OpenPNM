import scipy as _sp


def ordinary_diffusion(target,
                       pore_molar_density='pore.molar_density',
                       throat_molar_density='throat.molar_density',
                       pore_diffusivity='pore.diffusivity',
                       throat_diffusivity='throat.diffusivity',
                       pore_area='pore.area',
                       throat_area='throat.area',
                       pore_diameter='pore.diameter',
                       throat_length='throat.length',
                       shape_factor='throat.shape_factor',
                       throat_conduit_length='throat.conduit_lengths'):
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
    phase = target.project.find_phase(target)
    # Get Nt-by-2 list of pores connected to each throat
    P12 = network['throat.conns']
    parea = network[pore_area][P12]
    tarea = network[throat_area]
    try:  # If conduit lengths are defined, use them
        conduit = target[throat_conduit_length]
        plen1 = conduit[:, 0]
        tlen = conduit[:, 1]
        plen2 = conduit[:, 2]
        plen = _sp.vstack((plen1, plen2)).T
    except KeyError:  # Otherwise fall back to pore-throat-pore size info
        pdia = network[pore_diameter]
        plen1 = (0.5*pdia[P12[:, 0]])
        tlen = network[throat_length]
        plen2 = (0.5*pdia[P12[:, 1]])
        plen = _sp.vstack((plen1, plen2)).T
    # Interpolate pore phase property values to throats
    try:
        ct = phase[throat_molar_density]
    except KeyError:
        ct = phase.interpolate_data(propname=pore_molar_density)
    try:
        cp = phase[pore_molar_density][P12]
    except KeyError:
        cp = phase.interpolate_data(propname=throat_molar_density)
        cp = cp[P12]
    try:
        DABt = phase[throat_diffusivity]
    except KeyError:
        DABt = phase.interpolate_data(propname=pore_diffusivity)
    try:
        DABp = phase[pore_diffusivity][P12]
    except KeyError:
        DABp = phase.interpolate_data(propname=throat_diffusivity)
        DABp = DABp[P12]
    # Remove any non-positive lengths
    plen[plen <= 0] = 1e-12
    # Find g for each half of pore 1 and 2
    gp = (cp*DABp*parea) / plen
    gp[_sp.isnan(gp)] = _sp.inf
    # gp[~(gp > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat, remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except KeyError:
        sf = _sp.ones(network.num_throats())
    sf[_sp.isnan(sf)] = 1.0
    gt = (1/sf)*ct*DABt*tarea/tlen
    # Set 0 conductance pores (boundaries) to inf
    gt[~(gt > 0)] = _sp.inf
    value = (1/gt + 1/gp[:, 0] + 1/gp[:, 1])**(-1)
    value = value[phase.throats(target.name)]
    return value
