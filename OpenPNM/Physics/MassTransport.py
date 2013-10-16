
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def Chapman_Enskog_DiffusionCoefficient(T,P):
    
    
    return
    
def Fuller(T,P,MA=31.99,MB=28.01,VA=16.6,VB=17.9): 
    
    MAB = 2*(1/MA+1/MB)**(-1)    
    D = 0.00143*T**1.75/(P*MAB**0.5*(VA**(1/3)+VB**(1/3))**2)*1e-4    

    return D
    
def FullerScale(network,DABo,To=298,Po=101325): 
    r"""
    
    """
    T1 = network.pore_conditions['temperature']
    P1 = network.pore_conditions['pressure']
    D1 = DABo*(Po/To**1.75)*(T1**1.75/P1)

    return D1

def Conduits_DiffusionConductivity(network):
    r"""
    
    """
    R = 8.314462
    P = network.pore_conditions['pressure']
    T = network.pore_conditions['temperature']
    network.pore_conditions['molar_density'] = P/(R*T)
    network.interpolate_throat_conditions('molar_density')
    ct = network.throat_conditions['molar_density']
    cp = network.pore_conditions['molar_density']
    network.interpolate_throat_conditions('DAB')
    DABt = network.throat_conditions['DAB']
    DABp = network.pore_conditions['DAB']
        
    gt = ct*DABt*network.throat_properties['diameter']**2/network.throat_properties['length']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = cp*DABp*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    #Find g for half of pore 2
    gp2 = cp*DABp*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2) 
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    
    return g
    
