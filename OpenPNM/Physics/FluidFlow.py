
"""
module FluidFlow
===============================================================================


"""

import OpenPNM
import scipy as sp


def HydraulicConductance(network,fluid):
    r"""
    Calculates the hydraulic conductvity of throat assuming square geometry using a modified Hagen Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    try:
        mup = fluid.pore_conditions['viscosity']
    except:
        raise Exception('viscosity of the phase has not been specified')
    mut = fluid.interpolate_throat_conditions(network,mup)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*mut)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*mut)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = 2.28*(network.throat_properties['diameter']/2)**4/(2*network.throat_properties['length']*mut)
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    fluid.throat_conditions['hydraulic_conductance'] = g



