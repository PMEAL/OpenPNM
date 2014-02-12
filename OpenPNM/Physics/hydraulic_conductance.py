
"""
module hydraulic_conductance
===============================================================================

"""

import scipy as sp

def constant(physics,network,fluid,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(phase=fluid,prop=propname,data=value)

def na(physics,network,fluid,propname,**params):
    value = -1
    network.set_throat_data(phase=fluid,prop=propname,data=value)

def hagen_poiseuille(physics,network,fluid,propname,**params):
    r"""
    Calculates the hydraulic conductvity of throat assuming square geometry using a modified Hagen Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    mup = network.get_pore_data(phase=fluid,prop='viscosity')
    mut = network.interpolate_throat_data(mup)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.find_connected_pores(network.get_throat_data(prop='numbering'),flatten=0)
    #Find g for half of pore 1
    gp1 = 2.28*(network.get_pore_data(prop='diameter')[pores[:,0]]/2)**4/(network.get_pore_data(prop='diameter')[pores[:,0]]*mut)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.get_pore_data(prop='diameter')[pores[:,1]]/2)**4/(network.get_pore_data(prop='diameter')[pores[:,1]]*mut)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = 2.28*(network.get_throat_data(prop='diameter')/2)**4/(2*network.get_throat_data(prop='length')*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.set_throat_data(phase=fluid,prop=propname,data=value)



