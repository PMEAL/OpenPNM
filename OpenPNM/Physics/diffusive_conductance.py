
"""
module diffusive_conductance
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

def bulk_diffusion(physics,network,fluid,propname,**params):
    r"""
    Calculate the diffusive conductance of conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on the area

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties should already be calculated.

    """
    cp = network.get_pore_data(phase=fluid,prop='molar_density')
    DABp = network.get_pore_data(phase=fluid,prop='diffusivity')
    ct = network.interpolate_throat_data(cp)
    DABt = network.interpolate_throat_data(DABp)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.find_connected_pores(network.get_throat_data(prop='numbering'),flatten=0)
    #Find g for half of pore 1
    gp1 = ct*DABt*network.get_pore_data(prop='diameter')[pores[:,0]]**2/(0.5*network.get_pore_data(prop='diameter')[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*network.get_pore_data(prop='diameter')[pores[:,1]]**2/(0.5*network.get_pore_data(prop='diameter')[pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = ct*DABt*network.get_throat_data(prop='diameter')**2/(network.get_throat_data(prop='length'))
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.set_throat_data(phase=fluid,prop=propname,data=value)

