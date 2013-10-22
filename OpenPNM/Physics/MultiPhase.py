
"""
module MultiPhase
===============================================================================

"""

import OpenPNM
import scipy as sp

def conduit_filled_state_calculator(network):
    r"""
    
    """
    pores = network.get_connected_pores(network.throat_conditions['numbering'],flatten=0)
    network.throat_conditions['satn_wp_conduits'] = network.pore_conditions['satn_wp'][pores[:,0]]*network.pore_conditions['satn_wp'][pores[:,1]]
    network.throat_conditions['satn_nwp_conduits'] = -network.pore_conditions['satn_wp'][pores[:,0]]*-network.pore_conditions['satn_wp'][pores[:,1]]

def apply_phase_state_to_conduit_conductance(network,fluid_name):
    r"""
    nothing yet
    """
    fluid_wettability = network.phases[fluid_name]['wettability']    
    Conductance = network.throat_conditions['diffusive_conductance'+'_'+fluid_name]
    if fluid_wettability == 'wp':
        Conductance[network.throat_conditions['satn_nwp_conduits']] = 1e-30
    else:
        Conductance[network.throat_conditions['satn_wp_conduits']] = 1e-30
        
    return(Conductance)

def full_pore_filling(network,Pc=0.0,Seq=0):
    r"""
    Determine the filled state of a pore based on given capillary pressure

    Parameters
    ----------
    network : OpenPNM Network Object

    Pc : float, scalar
        The capillary pressure applied to the nonwetting phase

    Notes
    -----
    It is necessary that a capillary pressure curve has been run first, using the OrdinaryPercolation module.

    """
    if Pc:
        network.pore_conditions['satn_wp'] = network.pore_conditions['Pc_invaded']>Pc
    else:
        network.pore_conditions['satn_wp'] = network.pore_conditions['IP_inv_seq']>Seq
        
def late_pore_filling(network,swpi=0.0,eta=1.0,Pc=0.0):
    r"""
    Applies a late pore filling model to determine the fractional saturation of a pore based on the given capillary pressure

    Parameters
    ----------
    network : OpenPNM Network Object

    swpi : float, array_like
        The fraction of each pore still filled by wetting phase upon initial invasion

    eta : float, array_like
        The late pore filling exponent

    Pc : float, scalar
        The capillary pressure applied to the nonwetting phase

    Notes
    -----
    It is necessary that a capillary pressure curve has been run first, using the OrdinaryPercolation module.

    """
    try: swpi = network.pore_conditions['swpi']
    except: pass
    try: eta = network.pore_conditions['eta']
    except: pass

    Pc_star = network.pore_conditions['Pc_invaded']
    swp = swpi*(Pc_star/Pc)**eta*(network.pore_conditions['Pc_invaded']<=Pc)
    swp = swp + (network.pore_conditions['Pc_invaded']>Pc)
    network.pore_conditions['satn_wp'] = swp










