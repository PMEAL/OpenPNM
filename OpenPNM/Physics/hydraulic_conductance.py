r"""
===============================================================================
Submodule -- hydraulic_conductance
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
    fluid.set_throat_data(prop=propname,data=value,locations=geometry)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    value = -1
    fluid.set_throat_data(prop=propname,data=value,locations=geometry)

def hagen_poiseuille(physics,
                     network,
                     geometry,
                     fluid,
                     propname,
                     viscosity='viscosity',
                     throat_diameter = 'diameter',
                     throat_length = 'length',
                     pore_diameter = 'diameter',
                     **params):
    r"""
    Calculates the hydraulic conductivity of throat assuming square geometry using a modified Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    mut = fluid.get_data(prop='viscosity',throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    tind = network.get_throat_indices()
    pores = network.find_connected_pores(tind,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_pore_data(prop=pore_diameter)
    gp1 = 2.28*(pdia[pores[:,0]]/2)**4/(pdia[pores[:,0]]*mut)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(pdia[pores[:,1]]/2)**4/(pdia[pores[:,1]]*mut)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = 2.28*(tdia/2)**4/(2*tlen*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    mask = network.get_throat_indices(geometry)
    fluid.set_throat_data(prop=propname,data=value[mask],locations=geometry)


