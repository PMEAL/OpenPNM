r"""
===============================================================================
Submodule -- diffusive_conductance
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
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    r"""
    """
    value = -1
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def conduit_conductance(physics,
                   network,
                   fluid,
                   geometry,
                   conductance,
                   propname = 'conduit_conductance',
                   mode = 'strict',
                   factor = 1/1e3,
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
    #Interpolate pore fluid property values to throats
    ct = fluid.get_data(prop='molar_density',throats='all',mode='interpolate')
    DABt = fluid.get_data(prop='diffusivity',throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    Ts= network.get_throat_indices()
    Ps = network.find_connected_pores(Ts,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_data(prop=pore_diameter,pores='all')
    gp1 = ct*DABt*pdia[Ps[:,0]]**2/(0.5*pdia[Ps[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*pdia[Ps[:,1]]**2/(0.5*pdia[Ps[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_data(prop=throat_diameter,throats='all')
    tlen = network.get_data(prop=throat_length,throats='all')
    if (shape == 'circular'):
        gt = sp.pi*ct*DABt*tdia**2/(tlen*4)
    elif (shape == 'square'):
        gt = ct*DABt*tdia**2/tlen
    else:
        print 'invalid shape chosen.  Either circular or square'
        return
    gt = ct*DABt*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[geometry.throats()]
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

