
"""
module MultiPhase
===============================================================================

"""
import OpenPNM
import scipy as sp

def effective_occupancy(network,fluid,method='strict'):
    r"""

    """
    if method == 'strict':
        #if only EITHER pore is filled an open throat is considered closed
        pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
        fluid.throat_conditions['conduit_occupancy'] = fluid.pore_conditions['occupancy'][pores[:,0]]*fluid.pore_conditions['occupancy'][pores[:,1]]
        fluid.throat_conditions['conduit_occupancy'] = -fluid.pore_conditions['occupancy'][pores[:,0]]*-fluid.pore_conditions['occupancy'][pores[:,1]]
    elif method == 'moderate':
        #if only ONE pore isfilled an open throat is still considered open
        print 'nothing yet'
    elif method == 'liberal':
        #if BOTH pores are filled an open throat is still considered open
        print 'nothing yet'

def update_occupancy_OP(fluid,Pc=0):
    r"""
    ---
    """
    #Apply occupancy to given fluid
    try:
        fluid.pore_conditions['occupancy'] = sp.array(fluid.pore_conditions['Pc_invaded']<=Pc,ndmin=1)
        fluid.throat_conditions['occupancy'] = sp.array(fluid.throat_conditions['Pc_invaded']<=Pc,ndmin=1)
    except:
        print ('OP has not been run with this fluid, checking partner fluid')
        try:
            #Apply occupancy to given fluid
            fluid.pore_conditions['occupancy'] = sp.array(~(fluid.partner.pore_conditions['Pc_invaded']<=Pc),ndmin=1)
            fluid.throat_conditions['occupancy'] = sp.array(~(fluid.partner.throat_conditions['Pc_invaded']<=Pc),ndmin=1)
        except:
            raise Exception('It seems that OP has not been run on either fluid')
    #Apply occupancy to partner fluid
    fluid.partner.pore_conditions['occupancy'] = sp.array(~fluid.pore_conditions['occupancy'],ndmin=1)
    fluid.partner.throat_conditions['occupancy'] = sp.array(~fluid.throat_conditions['occupancy'],ndmin=1)

def update_occupancy_IP(network,fluid,Seq=0):
    r"""
    ---
    """
    try: 
        fluid.pore_conditions['occupancy'] = fluid.pore_conditions['IP_inv_seq']<Seq
        fluid.throat_conditions['occupancy'] = fluid.throat_conditions['IP_inv_seq']<Seq
    except: raise Exception('Something bad happened')
    try:
        fluid.partner.pore_conditions['occupancy'] = ~fluid.pore_conditions['occupancy']
        fluid.partner.throat_conditions['occupancy'] = ~fluid.throat_conditions['occupancy']
    except: raise Exception('A partner fluid has not been set so inverse occupancy cannot be set')

def late_pore_filling(network,fluid,swpi=0.0,eta=1.0,Pc=0.0):
    r"""
    Applies a late pore filling model to determine the fractional saturation of a pore based on the given capillary pressure

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object

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

    Pc_star = fluid.pore_conditions['Pc_invaded']
    swp = swpi*(Pc_star/Pc)**eta*(fluid.pore_conditions['Pc_invaded']<=Pc)
    swp = swp + (fluid.pore_conditions['Pc_invaded']>Pc)
    fluid.pore_conditions['volume_fraction'] = swp










