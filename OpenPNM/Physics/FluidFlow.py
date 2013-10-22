
"""
module FluidFlow
===============================================================================


"""

import OpenPNM
import scipy as sp


def HagenPoiseuille(network,fluid):
    r"""
    Calculates the hydraulic conductvity of throat assuming cylindrical geometry
    
    Parameters
    ----------
    network : OpenPNM Object
        The network for which the calculations are desired
        
    viscosity : float
        The viscosity of the fluid
    """
    
    vp = fluid['viscosity']
    vt = network.interpolate_throat_values(vp)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*vp[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*vp[:,1])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = 2.28*(network.throat_properties['diameter']/2)**4/(2*network.throat_properties['length']*vt)
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)

    fluid.update({'hydraulic_conductance': g})    


