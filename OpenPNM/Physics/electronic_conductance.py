
"""
module electronic_conductance
===============================================================================

"""

import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(physics,network,fluid,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(phase=fluid,prop=propname,data=value)

def na(physics,network,fluid,**params):
    value = -1
    network.set_throat_data(phase=fluid,prop=propname,data=value)

def parallel_resistors(physics,network,fluid,**params):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    sigmap = network.get_pore_data(phase=fluid,prop='electronic_conductivity')
    sigmat = network.interpolate_throat_data(sigmap)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.get_throat_data(prop='numbering'),flatten=0)
    #Find g for half of pore 1
    gp1 = (network.get_pore_data(prop='diameter')[pores[:,0]]/2)**2/(network.get_pore_data(prop='diameter')[pores[:,0]]*sigmat)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = (network.get_pore_data(prop='diameter')[pores[:,1]]/2)**2/(network.get_pore_data(prop='diameter')[pores[:,1]]*sigmat)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = sigmat*2*network.get_throat_data(prop='diameter')/(network.get_throat_data(prop='length'))
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.set_throat_data(phase=fluid,prop=propname,data=value)

