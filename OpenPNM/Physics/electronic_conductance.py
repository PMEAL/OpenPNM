
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
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    value = -1
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def series_resistors(physics,
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
    sigmat = fluid.get_data(prop=electrical_conductivity,throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    throats = network.get_throat_indices()
    pores = network.find_connected_pores(throats,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_data(prop=pore_diameter,pores='all')
    gp1 = sigmat*pdia[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = sigmat*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_data(prop=throat_diameter,throats='all')
    tlen = network.get_data(prop=throat_length,throats='all')
    gt = sigmat*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[geometry.throats()]
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

