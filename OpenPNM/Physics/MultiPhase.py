
"""
module MultiPhase
===============================================================================

"""
import scipy as sp

def calc_conduit_filling(network,method='strict'):
    r"""

    """
    if method == 'strict':
        #if only EITHER pore is filled an open throat is considered closed
        pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
        network.throat_conditions['conduit_satn_wp'] = network.pore_conditions['satn_wp'][pores[:,0]]*network.pore_conditions['satn_wp'][pores[:,1]]
        network.throat_conditions['conduit_satn_nwp'] = -network.pore_conditions['satn_wp'][pores[:,0]]*-network.pore_conditions['satn_wp'][pores[:,1]]
    elif method == 'moderate':
        #if only ONE pore isfilled an open throat is still considered open
        print 'nothing yet'
    elif method == 'liberal':
        #if BOTH pores are filled an open throat is still considered open
        print 'nothing yet'

def update_satn_from_OP(network,fluid,Pc=0):
    r"""
    ---
    """
    fluid_name = fluid['name']
    network.pore_conditions['satn_'+fluid_name] = network.pore_conditions['Pc_invaded']>Pc

def update_satn_from_IP(network,fluid,Seq=0):
    r"""
    ---
    """
    fluid_name = fluid['name']
    network.pore_conditions['satn_'+fluid_name] = network.pore_conditions['IP_inv_seq']>Seq

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
    network.pore_conditions['satn'+'_'+fluid_name] = swp










