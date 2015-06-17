r"""
===============================================================================
Submodule -- electrical_conductance
===============================================================================

"""

import scipy as _sp


def series_resistors(physics, phase, network,
                     pore_conductivity='pore.electrical_conductivity',
                     pore_area='pore.area', pore_diameter='pore.diameter',
                     throat_area='throat.area', throat_length='throat.length',
                     **kwargs):
    r"""
    Calculates the electrical conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object

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
    sigmap = phase[pore_conductivity]
    sigmat = phase.interpolate_data(sigmap)
    # Find g for half of pore 1
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Remove any non-positive lengths
    pdia[pdia <= 0] = 0
    gp1 = sigmap[Ps[:, 0]]*parea[Ps[:, 0]]/(0.5*pdia[Ps[:, 0]])
    # Find g for half of pore 2
    gp2 = sigmap[Ps[:, 1]]*parea[Ps[:, 1]]/(0.5*pdia[Ps[:, 1]])
    # Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 0
    gt = sigmat*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(physics.name)]
    return value
