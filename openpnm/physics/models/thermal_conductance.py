r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""

import scipy as _sp


def series_resistors(physics, phase, network,
                     thermal_conductivity='pore.thermal_conductivity',
                     pore_diameter='pore.diameter', pore_area='pore.area',
                     throat_area='throat.area', throat_length='throat.length',
                     **kwargs):
    r"""
    Calculate the thermal conductance of void conduits in network ( 1/2 pore - full
    throat - 1/2 pore ) based on size (assuming cylindrical geometry)

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
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    kp = phase[thermal_conductivity]
    kt = phase.interpolate_data(kp)
    # Find g for half of pore 1
    pdia = network[pore_diameter]
    parea = network[pore_area]
    pdia1 = pdia[Ps[:, 0]]
    pdia2 = pdia[Ps[:, 1]]
    # Remove any non-positive lengths
    pdia1[pdia1 <= 0] = 1e-12
    pdia2[pdia2 <= 0] = 1e-12
    gp1 = kt*parea[Ps[:, 0]]/(0.5*pdia1)
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = kt*parea[Ps[:, 1]]/(0.5*pdia2)
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    gt = kt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(physics.name)]
    return value
