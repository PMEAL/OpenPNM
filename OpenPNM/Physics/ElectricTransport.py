# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:51:40 2013

@author: ARAX
"""


import OpenPNM
import scipy as sp


def ElectricalConductivity(network,sigma):
    r"""
    ---
    """
    gt = sigma*2*network.throat_properties['diameter']/network.lattice_spacing
    
    return gt
def Torey_electrical_conductivity(network):
    
    k = -4.91*-11*network.pore_properties['temperature']**3 + 1.42e-8*network.pore_properties['temperature']**2\
        -1.46e-6*network.pore_properties['temperature'] + 8.91e-5
    return k