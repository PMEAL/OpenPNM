r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as sp


def conduit_conductance(physics, phase, network, throat_conductance,
                        throat_occupancy='throat.occupancy',
                        pore_occupancy='pore.occupancy',
                        mode='strict', factor=1e-6, **kwargs):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.

    This method "closes" conduits that are not sufficiently filled with the
    specified phase by multiplying the original conductance by a very small *factor*.

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
        The phase of interest

    occupied_condition : 'occupancy'
        The name of the pore and throat property that dictates whether conduit is
        "closed" or not

    mode : 'strict' or 'medium' or 'loose'
        How agressive the method should be in "closing" conduits.
        'strict' implies that if any pore or throat in the conduit is unoccupied by
         the given phase, the conduit is closed.
        'medium' implies that if either the throat or both pores are unoccupied, the
        conduit is closed
        'loose' will only close the conduit if the throat is unoccupied.

    factor : 1/1e3
        The "closing factor" which becomes multiplied to the original conduit's
        conductance to severely limit transport.

    Notes
    -----
    This function requires that all the necessary phase properties already be
    calculated.

    """
    throats = phase.throats(physics.name)
    if mode == 'loose':
        closed_conduits = -sp.array(phase[throat_occupancy], dtype=bool)
    else:
        throats_closed = -sp.array(phase[throat_occupancy], dtype=bool)
        connected_pores = network.find_connected_pores(throats)
        pores_1 = connected_pores[:, 0]
        pores_2 = connected_pores[:, 1]
        pores_1_closed = -sp.array(phase[pore_occupancy][pores_1], dtype=bool)
        pores_2_closed = -sp.array(phase[pore_occupancy][pores_2], dtype=bool)
        if mode == 'medium':
            closed_conduits = throats_closed | (pores_1_closed & pores_2_closed)

        if mode == 'strict':
            closed_conduits = pores_1_closed | throats_closed | pores_2_closed
    open_conduits = -closed_conduits
    throat_value = phase[throat_conductance]
    value = throat_value*open_conduits + throat_value*closed_conduits*factor
    value = value[phase.throats(physics.name)]
    return value


def late_throat_filling(network, phase, Pc_star, eta):
    pass


def late_pore_filling(physics, phase, network, Pc, Swp_star=0.2, eta=3,
                      wetting_phase=False, pore_occupancy='pore.occupancy',
                      throat_capillary_pressure='throat.capillary_pressure',
                      **kwargs):
    r"""
    Applies a late pore filling model to calculate fractional pore filling as
    a function of applied capillary pressure.

    Parameters
    ----------
    Pc : float
        The capillary pressure in the non-wetting phase (Pc > 0)
    Swp_star : float
        The residual wetting phase in an invaded pore immediately after
        nonwetting phase invasion
    eta : float
        Exponent to control the rate at which wetting phase is displaced
    wetting_phase : boolean
        Indicates whether supplied phase is the wetting or non-wetting phase


    """
    pores = phase.pores(physics.name)
    prop = phase[throat_capillary_pressure]
    neighborTs = network.find_neighbor_throats(pores, flatten=False)
    Pc_star = sp.array([sp.amin(prop[row]) for row in neighborTs])
    Swp = Swp_star*(Pc_star/Pc)**eta
    if wetting_phase:
        values = Swp*phase[pore_occupancy]*(Pc_star < Pc)
    else:
        values = (1-Swp)*(1-phase[pore_occupancy])*(Pc_star < Pc)
    values = values[phase.throats(physics.name)]
    return values
