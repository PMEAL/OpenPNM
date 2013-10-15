
"""
module FluidFlow
===============================================================================


"""

import OpenPNM
import scipy as sp


def HagenPoiseuille(network,viscosity):
    r"""
    Calculates the hydraulic conductvity of throat assuming cylindrical geometry
    
    Parameters
    ----------
    network : OpenPNM Object
        The network for which the calculations are desired
        
    viscosity : float
        The viscosity of the fluid
    """
    gt = 2.28*(network.throat_properties['diameter']/2)**4/(2*network.throat_properties['length']*viscosity)
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*viscosity)
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*viscosity) 
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return g


