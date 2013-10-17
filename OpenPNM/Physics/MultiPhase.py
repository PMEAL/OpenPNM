
"""
module MultiPhase
===============================================================================

"""

import OpenPNM
import scipy as sp

def Conduit_Filled_State_Calculator(network):
    r"""
    nothing yet
    """
    if 'Pc_invaded' in network.pore_properties:
        network.pore_properties['swp'] = network.pore_properties['Pc_invaded']>p_val
        network.pore_properties['snwp'] = -network.pore_properties['swp']
    elif 'IP_Pseq' in network.pore_properties:
        network.pore_properties['swp'] = network.pore_properties['IP_Pseq']>p_val
        network.pore_properties['snwp'] = -network.pore_properties['swp']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    network.throat_properties['swp_conduits'] = -(-network.pore_properties['swp'][pores[:,0]]*-network.pore_properties['swp'][pores[:,1]])
    network.throat_properties['snwp_conduits'] = -(-network.pore_properties['snwp'][pores[:,0]]*-network.pore_properties['snwp'][pores[:,1]])

#    if -sp.in1d(neighborPs,self.Pinvaded).all():
#        self._net.throat_properties['UninvadedConduits'] = 1
#            val_name = 'Pc_invaded'
#
#        elif self.Alg=='IP':
#            Alg_var = self.Psequence
#            val_name = 'IP_Pseq'
#
#        for P_val in Alg_var:
#            self._logger.info("Applying Pressure/Sequence = "+str(P_val))
#            Concentration = self._do_one_inner_iteration(P_val)
#            #Store result of Concentration in each step
#            if P_val!=0:
#                Concentration = np.multiply(Concentration,self._net.pore_properties[val_name]>P_val)
#
#            self._net.set_pore_property(name="Concentration_at_"+str(P_val),ndarray=Concentration)
#    return

def Apply_Phase_State_to_Conduit_Conductivity(network):
    r"""
    nothing yet
    """
    C_wet = network.throat_properties['Cdiff']
    C_wet[network.throat_properties['snwp']] = 1e-30
    return(C_wet)

def full_pore_filling(network,Pc=0.0):
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
    network.pore_conditions['satn_wp'] = network.pore_conditions['Pc_invaded']>Pc

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

    .. warning:: If the values of eta and swpi are not already set in the network then the parmeters received by this function will be written to the network. All other calls to this function will use these values regardless of what parameters are received. These parameters must be overwritten or deleted explicity.

    """

    try:    swpi = network.pore_conditions['swpi']
    except: network.pore_conditions['swpi'] = swpi

    try:    eta = network.pore_conditions['eta']
    except: network.pore_conditions['eta'] = eta

    Pc_star = network.pore_conditions['Pc_invaded']
    swp = swpi*(Pc_star/Pc)**eta*(network.pore_conditions['Pc_invaded']<=Pc)
    swp = swp + (network.pore_conditions['Pc_invaded']>Pc)
    network.pore_conditions['satn_wp'] = swp










