
"""
module FluidFlow
===============================================================================


"""

import OpenPNM
import scipy as sp


def HydraulicConductance(network,fluid_name):
    r"""
    Calculates the hydraulic conductvity of throat assuming cylindrical geometry using the Hagen Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid_name : 'string'
    """
    try:
        mup = network.pore_conditions['viscosity'+'_'+fluid_name]
    except:
        raise Exception('viscosity of the '+fluid_name+' phase has not been specified')
    mut = network.interpolate_throat_values(mup)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1  
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*mup)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*mup)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = 2.28*(network.throat_properties['diameter']/2)**4/(2*network.throat_properties['length']*mut)
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions['hydraulic_conductance'+'_'+fluid_name] = g



