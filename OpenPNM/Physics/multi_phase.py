
"""
module multi_phase
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

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
        print('nothing yet')
    elif method == 'liberal':
        #if BOTH pores are filled an open throat is still considered open
        print('nothing yet')

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

    Pc_star = network.get_pore_data(fluid=fluid.name,prop='Pc_invaded')
    swp = swpi*(Pc_star/Pc)**eta*(Pc_star<=Pc)
    swp = swp + (Pc_star>Pc)
    network.set_throat_data(fluid=fluid.name,prop='volume_fraction',data=swp)










