# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:51:40 2013

@author: ARAX
"""


import OpenPNM
import scipy as sp


def ElectronicConductance(network,fluid_name):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid_name : 'string'
    """
    try:
        sigmap = network.pore_properties['electonic_conducty'+'_'+fluid_name]
    except:
        raise Exception('Viscosity of the '+fluid_name+' phase has not been specified')
    sigmat = network.interpolate_throat_values(sigmap)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = 2.28*(network.pore_properties['diameter'][pores[:,0]]/2)**4/(network.pore_properties['diameter'][pores[:,0]]*sigmap[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = 2.28*(network.pore_properties['diameter'][pores[:,1]]/2)**4/(network.pore_properties['diameter'][pores[:,1]]*sigmap[:,1])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = sigmat*2*network.throat_properties['diameter']/network.lattice_spacing
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions['hydraulic_conductance'+'_'+fluid_name] = g


    return gt
def Torey_electrical_conductivity(network):

    k = -4.91*-11*network.pore_properties['temperature']**3 + 1.42e-8*network.pore_properties['temperature']**2\
        -1.46e-6*network.pore_properties['temperature'] + 8.91e-5
    return k