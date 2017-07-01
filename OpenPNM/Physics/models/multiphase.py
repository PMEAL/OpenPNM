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
    Add a new multiphase conductance property to the conduits of network, where
    a conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.

    This method "closes" conduits that are not sufficiently filled with the
    specified phase by multiplying the original conductance by a very small
    factor.

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
        The phase of interest

    occupied_condition : 'occupancy'
        The name of the pore and throat property that dictates whether conduit
        is "closed" or not

    mode : 'strict' or 'medium' or 'loose'
        How agressive the method should be in "closing" conduits.

        **'strict'** :  If any pore or throat in the conduit is unoccupied by
         the given phase, the conduit is closed.

        **'medium'** : If either the throat or both pores are unoccupied, the
        conduit is closed

        **'loose'** : Only close the conduit if the throat is unoccupied

    factor : float (default is 1e-6)
        The factor which becomes multiplied to the original conduit's
        conductance to severely limit transport, but not set it to zero.

    Notes
    -----
    This function requires that all the necessary phase properties already be
    calculated.

    """
    throats = phase.Ts
    if mode == 'loose':
        closed_conduits = ~sp.array(phase[throat_occupancy], dtype=bool)
    else:
        throats_closed = ~sp.array(phase[throat_occupancy], dtype=bool)
        connected_pores = network.find_connected_pores(throats)
        pores_1 = connected_pores[:, 0]
        pores_2 = connected_pores[:, 1]
        pores_1_closed = ~sp.array(phase[pore_occupancy][pores_1], dtype=bool)
        pores_2_closed = ~sp.array(phase[pore_occupancy][pores_2], dtype=bool)
        if mode == 'medium':
            closed_conduits = throats_closed | (pores_1_closed & pores_2_closed)

        if mode == 'strict':
            closed_conduits = pores_1_closed | throats_closed | pores_2_closed
    open_conduits = ~closed_conduits
    throat_value = phase[throat_conductance]
    value = throat_value*open_conduits + throat_value*closed_conduits*factor
    value = value[phase.throats(physics.name)]
    return value


def late_throat_filling(network, phase, physics,
                        Pc,
                        Swp_star=0.11,
                        eta=3,
                        throat_entry_pressure='throat.capillary_pressure',
                        **kwargs):
    r"""
    Applies a late thraot filling model to calculate fractional throat filling
    as a function of applied capillary pressure.

    Parameters
    ----------
    Pc : float
        The capillary pressure in the non-wetting phase (Pc > 0)

    eta : float
        Exponent to control the rate at which wetting phase is displaced

    Swp_star : float
        The residual wetting phase in an invaded pore immediately after
        nonwetting phase invasion

    throat_entry_pressure : string
        The dictionary key containing throat entry pressures.  The default is
        'throat.capillary_pressure'.

    Returns
    -------
    A Nt-list of containing the fraction of each throat that is filled with
    non-wetting phase.  Note this method does NOT account for whether a throat
    is actually filled or not; this needs to be done using some other external
    criteria such as the 'throat.inv_Pc' array on a *Drainage* algorithm.

    """
    Swp = sp.ones([physics.Nt, ])
    if Pc > 0:
        Swp = Swp_star*(physics[throat_entry_pressure]/Pc)**eta
    values = (1-Swp)*(physics[throat_entry_pressure] <= Pc)
    return values


def late_pore_filling(physics, phase, network,
                      Pc,
                      Swp_star=0.2,
                      eta=3,
                      pc_star='pore.pc_star',
                      throat_entry_pressure='throat.capillary_pressure',
                      **kwargs):
    r"""
    Applies a late pore filling model to calculate fractional pore filling as
    a function of applied capillary pressure.

    Parameters
    ----------
    Pc : float
        The capillary pressure in the non-wetting phase (Pc > 0)

    eta : float
        Exponent to control the rate at which wetting phase is displaced

    Swp_star : float
        The residual wetting phase in an invaded pore immediately after
        nonwetting phase invasion

    pc_star : string
        The dictionary key to find or place the Pc_star array.  Pc_star is
        the minimum pressure at which a pore can be invaded and is found as
        the minimum entery pressure of all the pore's neighboring throats.
        The default is 'pore.Pc_star' and if this array is not found it is
        created.

    throat_entry_pressure : string
        The dictionary key containing throat entry pressures.  The default is
        'throat.capillary_pressure'.

    Returns
    -------
    A Np-list of containing the fraction of each pore that is filled with non-
    wetting phase.  Note this method does NOT account for whether a pore is
    actually filled or not; this needs to be done using some other external
    criteria such as the 'pore.inv_Pc' array on a *Drainage* algorithm.

    """
    # If pc_star has not yet been calculated, do so
    if pc_star not in physics.keys():
        pores = phase.Ps
        prop = phase[throat_entry_pressure]
        neighborTs = network.find_neighbor_throats(pores, flatten=False)
        temp = sp.array([sp.amin(prop[row]) for row in neighborTs])
        physics[pc_star] = temp[physics.Pnet]

    Swp = sp.ones([physics.Np, ])
    if Pc > 0:
        Swp = Swp_star*(physics[pc_star]/Pc)**eta
    else:
        Swp = Swp_star
    values = (1-Swp)*(physics[pc_star] <= Pc)
    return values
