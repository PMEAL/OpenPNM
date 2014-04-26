
"""
module electronic_conductance
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

def parallel_resistors(physics,
                       network,
                       geometry,
                       fluid,
                       propname,
                       electrical_conductivity='electrical_conductivity',
                       throat_diameter = 'diameter',
                       throat_length = 'length',
                       pore_diameter = 'diameter',                       
                       **params):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    sigmap = fluid.get_pore_data(prop=electrical_conductivity)
    sigmat = network.interpolate_throat_data(sigmap)
    #Get Nt-by-2 list of pores connected to each throat
    tind = network.get_throat_indices()
    pores = network.find_connected_pores(tind,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_pore_data(prop=pore_diameter)
    gp1 = sigmat*pdia[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = sigmat*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = sigmat*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    mask = network.get_throat_indices(geometry)
    fluid.set_throat_data(prop=propname,data=value[mask],locations=geometry)

