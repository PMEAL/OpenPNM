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
    Swp = sp.ones(physics.Nt,)
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

    Swp = sp.ones(physics.Np,)
    if Pc > 0:
        Swp = Swp_star*(physics[pc_star]/Pc)**eta
    else:
        Swp = Swp_star
    values = (1-Swp)*(physics[pc_star] <= Pc)
    return values

def phase_conductance(physics, phase, network,
           pore_viscosity='pore.viscosity',
           throat_length='throat.length',
           throat_diameter='throat.diameter',
           shape_factor='throat.shape_factor',
           contact_angle='pore.contact_angle',
           throat_occupancy='throat.occupancy',
           pore_surface_tension='pore.surface_tension',
           throat_Pc='throat.capillary_pressure',
           **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming cuboid
    geometry using the Mayer and Stowe approximation where both phases
    may be present at the same time
    
    *** For now asssume that one phase is wetting and one-is non-wetting ***
    *** Mixed wettability may not work as saturation of wetting phase is ***
    *** Solved for in dynamic pressure solver                            ***

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
    # Get properties in every pore in the network
    mup = phase[pore_viscosity]
    mut = phase.interpolate_data(mup)
    sigmap = phase[pore_surface_tension]
    sigmat = phase.interpolate_data(sigmap)
    # Determine the throats where phase is wetting
    theta = sp.mean(phase[contact_angle])
    wetting_phase = theta <= 90
    beta = 0.5
    # Find g for full throat
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except:
        sf = sp.ones(network.num_throats())
    sf[sp.isnan(sf)] = 1.0
    value = sp.zeros(network.num_throats())
    # Geometric throat radius
    r = network[throat_diameter]/2
    # Calculate effective radius single-phase
    r_eff_sp = sp.sqrt(4/sp.pi)*r
    # Calculate effective radius capillary pressure
    r_c = sigmat/sp.absolute(physics[throat_Pc])
    # Calculate effective radius single-phase
    r_eff_mp = 0.5*(sp.sqrt((r**2 - (4-sp.pi)*r_c**2)/sp.pi)+r)
    if wetting_phase:
        # Throats  occupied by wetting phase completely
        Ts = sp.array(phase[throat_occupancy]==1.0, dtype=bool)
        v1 = (sp.pi*r_eff_sp**4)/(8*mut*tlen)
        value[Ts] = v1[Ts]
        # Throats occupied by both phases
        v2 = ((4 -sp.pi)*r_c**4)/(beta*mut*tlen)
        value[~Ts] = v2[~Ts]
    else:
        # Throats occupied by non-wetting phase to some extent
        Ts = sp.array(phase[throat_occupancy]>0.0, dtype=bool)
        v1 = (sp.pi*r_eff_mp**4)/(8*mut*tlen)
        value[Ts] = v1[Ts]
        # Throats occupied by wetting-phase completely have zero conduction 
        # to non-wetting phase

    value = value[phase.throats(physics.name)]
    return value

def lpf_cuboid(physics, phase, network, Pc,
               pore_diameter='pore.diameter',
               **kwargs):
    r"""
    Applies a late pore filling model to calculate fractional pore filling as
    a function of applied capillary pressure.

    Parameters
    ----------
    Pc : float
        The capillary pressure in the non-wetting phase (Pc > 0)

    wetting_phase : boolean
        Indicates whether supplied phase is the wetting or non-wetting phase


    """
    sigma = sp.mean(phase["pore.surface_tension"])
    values = -sp.log(1-(4*sigma/(network[pore_diameter]*Pc)))/6.83
    values = values[phase.pores(physics.name)]
    return values
