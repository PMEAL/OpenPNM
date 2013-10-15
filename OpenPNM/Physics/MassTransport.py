
"""
module MassTransport
===============================================================================

"""

import OpenPNM
import scipy as sp

def Chapman_Enskog_DiffusionCoefficient(T,P):
    
    
    return
    
def Fuller_DiffusionCoefficient(T,P,MA=31.99,MB=28.01,VA=16.6,VB=17.9): 
    
    MAB = 2*(1/MA+1/MB)**(-1)    
    D = 0.00143*T**1.75/(P*MAB**0.5*(VA**(1/3)+VB**(1/3))**2)*1e-4    

    return D

def Conduits_DiffusionConductivity(network,T,P,DAB_method):
    r"""
    This does not account for temperature gradients
    Binary_Diff_In_a_Half_Pore_Throat_Half_Pore_Conduit
    """
    R = 8.314462    
    C = P/(R*T)
    if DAB_method=='Fuller':
        DAB = OpenPNM.Physics.MassTransport.Fuller_DiffusionCoefficient(T,P)
    elif DAB_method=='Chapman_Enskog':
        DAB = OpenPNM.Physics.MassTransport.Chapman_Enskog_DiffusionCoefficient(T,P)
        
    gt = C*DAB*network.throat_properties['diameter']**2/network.throat_properties['length']
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    gp1 = C*DAB*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp2 = C*DAB*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2) 
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return g
    
