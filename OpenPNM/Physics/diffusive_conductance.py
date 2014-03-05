
"""
module diffusive_conductance
===============================================================================

"""

import scipy as sp

def constant(physics,
             network,
             geometry,
             fluid,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_throat_data(phase=fluid,prop=propname,data=value,locations=geometry)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    r"""
    """
    value = -1
    network.set_throat_data(phase=fluid,prop=propname,data=value,locations=geometry)

def bulk_diffusion(physics,
                   network,
                   fluid,
                   geometry,
                   propname,
                   diffusivity = 'diffusivity',
                   molar_density = 'molar_density',
                   throat_diameter = 'diameter',
                   throat_length = 'length',
                   pore_diameter = 'diameter',
                   **params):
    r"""
    Calculate the diffusive conductance of conduits in network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties already be 
    calculated.

    """    
    #Get fluid properties
    cp = network.get_pore_data(phase=fluid,prop=molar_density)
    DABp = network.get_pore_data(phase=fluid,prop=diffusivity)
    #Interpolate pore values to throats
    ct = network.interpolate_throat_data(cp)
    DABt = network.interpolate_throat_data(DABp)
    #Get Nt-by-2 list of pores connected to each throat
    tind = network.get_throat_indices()
    pores = network.find_connected_pores(tind,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_pore_data(prop=pore_diameter)
    gp1 = ct*DABt*pdia[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = ct*DABt*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    mask = network.get_throat_indices(geometry)
    network.set_throat_data(phase=fluid,prop=propname,data=value[mask],locations=geometry)

