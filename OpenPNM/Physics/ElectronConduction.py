# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:51:40 2013

@author: ARAX
"""


import OpenPNM
import scipy as sp


def ElectronicConductance(network,fluid):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    try:
        sigmap = fluid.pore_conditions['electronic_conductivity']
    except:
        raise Exception('electronic_conductivity of the phase has not been specified')
    sigmat = fluid.interpolate_throat_conditions(network,sigmap)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = (network.pore_properties['diameter'][pores[:,0]]/2)**2/(network.pore_properties['diameter'][pores[:,0]]*sigmat)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = (network.pore_properties['diameter'][pores[:,1]]/2)**2/(network.pore_properties['diameter'][pores[:,1]]*sigmat)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = sigmat*2*network.throat_properties['diameter']/(network.throat_properties['length'])
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    fluid.throat_conditions['electronic_conductance'] = g

